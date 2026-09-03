#!/usr/bin/env python3
"""
Ground-state search for a stealthy point pattern with a hole (void) constraint,
in a 2D periodic box.  Fully vectorized (numpy/BLAS) + parallel (BLAS threads,
cKDTree workers).  Drop-in replacement for CCO.out on the ground path: same CLI
args and stdin token order, with the hole parameters appended.

Objective (minimized over particle positions):
    Phi = Phi_ann + Phi_hole
  Phi_ann  = sum_{k in annuli} |rho_k|^2 / N,   rho_k = sum_j exp(i k.r_j)
  Phi_hole = (lambda_h / M^2) sum_alpha ( max(0, d_alpha - Rmax) / Rmax )^2
    where d_alpha = periodic distance from grid point alpha to its nearest particle.

CLI :  python3 ground_search.py <timelimit_hr> <seed> <beg_idx> <verbosity>
stdin (whitespace tokens, same order the ground slurm already builds):
    dim  M  (K1a deltaa S0)*M  vareps0 sigma phi_fake  threads N Nc
    initmode  savename  mode  [tolerance <v>] [maxsteps <v>] [algorithm <name>] run
    Rmax  lambda_h  <grid sizes ...>          # hole params (appended)
The K1a/deltaa are scaled by a=(phi_fake/V_d)^(1/d) exactly as the slurm builds
them; this script divides by a to recover absolute k (so the slurm is unchanged).

Output: <savename>.h5  (config_i datasets, box + parameters + per-config diagnostics).
"""

import os
import sys
import time
import math

import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
import h5py


class _TimeUp(Exception):
    """Raised from the L-BFGS callback to abort a stage when the wall limit is hit."""
    pass


# ---------------------------------------------------------------- k-space term
def phi_ann_and_grad(r, kvecs, N):
    """Phi_ann = sum_k |rho_k|^2 / N  and its gradient (N,2). Vectorized (BLAS)."""
    kr = kvecs @ r.T                 # (K, N)
    P = np.exp(1j * kr)              # (K, N)
    rho = P.sum(axis=1)              # (K,)
    phi = float(np.vdot(rho, rho).real) / N
    A = np.conj(rho)[:, None] * P    # (K, N)
    grad = -(2.0 / N) * (kvecs.T @ A.imag).T   # (N, 2)
    return phi, grad


def structure_factor_max(r, kvecs, N):
    kr = kvecs @ r.T
    rho = np.exp(1j * kr).sum(axis=1)
    S = (np.abs(rho) ** 2) / N
    return float(S.max()) if S.size else 0.0


# ------------------------------------------------------------- real-space term
def make_grid(L, M):
    c = (np.arange(M) + 0.5) * (L / M)      # cell centers in [0, L)
    gx, gy = np.meshgrid(c, c)
    return np.stack([gx.ravel(), gy.ravel()], axis=1)


def phi_hole_and_grad(r, L, grid_pts, Rmax, lambda_h, M, workers):
    """Phi_hole and its gradient (N,2) using periodic nearest-neighbor (cKDTree)."""
    rw = np.mod(r, L)
    tree = cKDTree(rw, boxsize=L)
    d, idx = tree.query(grid_pts, k=1, workers=workers)     # periodic NN
    u = np.maximum(0.0, d - Rmax)
    phi = (lambda_h / (M * M)) * float(np.sum((u / Rmax) ** 2))
    grad = np.zeros_like(r)
    active = d > Rmax
    if np.any(active):
        da, ida, ga = d[active], idx[active], grid_pts[active]
        disp = rw[ida] - ga
        disp -= L * np.round(disp / L)                      # min-image
        coef = (2.0 * lambda_h / (M * M * Rmax * Rmax)) * (da - Rmax) / da
        np.add.at(grad, ida, coef[:, None] * disp)          # scatter onto j_alpha
    return phi, grad


def hole_diagnostics(r, L, grid_pts, Rmax, M, workers):
    rw = np.mod(r, L)
    tree = cKDTree(rw, boxsize=L)
    d, _ = tree.query(grid_pts, k=1, workers=workers)
    return float(d.max()), float(np.mean(d > Rmax))         # max_d, fraction violating


# ------------------------------------------------------------------- objective
def make_objective(kvecs, N, L, grid_pts, Rmax, lambda_h, M, workers):
    shape = (N, 2)

    def obj(xflat):
        r = xflat.reshape(shape)
        pa, ga = phi_ann_and_grad(r, kvecs, N)
        ph, gh = phi_hole_and_grad(r, L, grid_pts, Rmax, lambda_h, M, workers)
        return pa + ph, (ga + gh).ravel()

    return obj


def report(tag, r, kvecs, N, L, Rmax, lambda_h, M, workers, out=sys.stdout):
    pa, _ = phi_ann_and_grad(r, kvecs, N)
    grid = make_grid(L, M)
    ph, _ = phi_hole_and_grad(r, L, grid, Rmax, lambda_h, M, workers)
    maxS = structure_factor_max(r, kvecs, N)
    max_d, frac = hole_diagnostics(r, L, grid, Rmax, M, workers)
    out.write(f"  [{tag}] M={M}:  Phi_ann={pa:.6e}  Phi_hole={ph:.6e}"
              f"  maxS={maxS:.6e}  max_d={max_d:.6f} (Rmax={Rmax})"
              f"  frac(d>Rmax)={frac:.4f}\n")
    out.flush()
    return dict(phi_ann=pa, phi_hole=ph, maxS=maxS, max_d=max_d, frac=frac)


# ----------------------------------------------------------------------- input
def parse_stdin(tokens):
    it = iter(tokens)
    nxt = lambda: next(it)
    dim = int(nxt())
    M = int(nxt())
    shells = [(float(nxt()), float(nxt()), float(nxt())) for _ in range(M)]  # (K1a,deltaa,S0)
    vareps0 = float(nxt()); sigma = float(nxt()); phi_fake = float(nxt())
    threads = int(nxt()); N = int(nxt()); Nc = int(nxt())
    initmode = nxt(); savename = nxt(); mode = nxt()
    tol, maxsteps, algorithm = 1e-16, 100000, "LBFGS"
    while True:
        t = nxt()
        if t == "run":
            break
        tl = t.lower()
        if tl == "tolerance":   tol = float(nxt())
        elif tl == "maxsteps":  maxsteps = int(nxt())
        elif tl == "algorithm": algorithm = nxt()
    Rmax = float(nxt()); lambda_h = float(nxt())
    grid_sched = [int(t) for t in it] or [64, 128, 256]
    return dict(dim=dim, shells=shells, vareps0=vareps0, sigma=sigma, phi_fake=phi_fake,
                threads=threads, N=N, Nc=Nc, initmode=initmode, savename=savename,
                mode=mode, tol=tol, maxsteps=maxsteps, algorithm=algorithm,
                Rmax=Rmax, lambda_h=lambda_h, grid_sched=grid_sched)


def unit_ball_volume(d):
    return math.pi ** (d / 2.0) / math.gamma(d / 2.0 + 1.0)


def build_kvecs(shells_abs, L, dim):
    """Independent reciprocal vectors k=(2pi/L)*n with k1 < |k| <= k2 for some shell.

    Only ONE of each +-k pair is kept (upper half-plane), i.e. the independent
    collective-coordinate constraints -- matching the C++ GetKs/CCO convention.
    So Phi_ann = sum over these modes, and chi = len(kvecs)/(d(N-1)) is correct.
    """
    if dim != 2:
        raise SystemExit("ground_search.py currently supports dim=2 only.")
    k2max = max(k2 for (_, k2) in shells_abs)
    nmax = int(math.ceil(k2max * L / (2.0 * math.pi))) + 1
    ns = np.arange(-nmax, nmax + 1)
    nx, ny = np.meshgrid(ns, ns)
    nx, ny = nx.ravel(), ny.ravel()
    kv = (2.0 * math.pi / L) * np.stack([nx, ny], axis=1).astype(float)
    km = np.linalg.norm(kv, axis=1)
    in_ann = np.zeros(km.shape, dtype=bool)
    for (k1, k2) in shells_abs:
        in_ann |= (km > k1 + 1e-12) & (km <= k2 + 1e-12)
    half = (nx > 0) | ((nx == 0) & (ny > 0))       # one per +-k pair (independent)
    return kv[in_ann & half]


# ---------------------------------------------------------------------- driver
def main():
    timelimit_hr = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    # argv[3] = beg_idx (unused here), argv[4] = verbosity
    verbosity = int(sys.argv[4]) if len(sys.argv) > 4 else 3

    p = parse_stdin(sys.stdin.read().split())
    dim, N, L = p["dim"], p["N"], p["N"] ** (1.0 / p["dim"])   # unit number density
    workers = max(1, p["threads"])
    t_stop = time.time() + timelimit_hr * 3600 - 300           # 5-min margin

    # recover absolute annuli: divide the a-scaled shells by a=(phi_fake/V_d)^(1/d)
    a = (p["phi_fake"] / unit_ball_volume(dim)) ** (1.0 / dim)
    shells_abs, n_equi = [], 0
    for (K1a, deltaa, S0) in p["shells"]:
        if S0 != 0.0:
            n_equi += 1
            continue                     # equiluminous (S0>0) not handled in v1
        shells_abs.append((K1a / a, (K1a + deltaa) / a))
    if n_equi:
        print(f"WARNING: {n_equi} equiluminous shell(s) (S0>0) ignored; v1 is stealthy-only.")
    if not shells_abs:
        raise SystemExit("No stealthy (S0=0) shells to constrain.")

    kvecs = build_kvecs(shells_abs, L, dim)

    print(f"ground_search: N={N}, L={L:.6f} (rho=1), dim={dim}, threads={workers}", flush=True)
    print(f"  annuli (absolute k): {['[%.4f,%.4f]' % s for s in shells_abs]}", flush=True)
    print(f"  constrained modes: {len(kvecs)}   chi ~ {len(kvecs)/(dim*(N-1)):.4f}", flush=True)
    print(f"  Rmax={p['Rmax']}, lambda_h={p['lambda_h']}, grid schedule={p['grid_sched']}", flush=True)
    print(f"  Nc={p['Nc']}, maxsteps={p['maxsteps']}, ftol(tolerance)={p['tol']}\n", flush=True)
    if p["algorithm"].upper() not in ("LBFGS", "L-BFGS", "L-BFGS-B"):
        print(f"WARNING: algorithm '{p['algorithm']}' not supported; using L-BFGS-B.", flush=True)

    rng = np.random.default_rng(seed)
    configs, diags = [], []
    stop_all = False
    for c in range(p["Nc"]):
        if time.time() > t_stop:
            print(f"Time limit reached; produced {c}/{p['Nc']} configs.", flush=True)
            break
        print(f"--- config {c} ---", flush=True)
        r = rng.random((N, dim)) * L                          # random initial condition
        d_final = None
        for M in p["grid_sched"]:                             # refine 64 -> 128 -> 256
            grid = make_grid(L, M)
            obj = make_objective(kvecs, N, L, grid, p["Rmax"], p["lambda_h"], M, workers)
            state = {"x": r.ravel().copy()}
            def cb(xk, _s=state):                             # abort mid-stage on wall limit
                _s["x"] = np.asarray(xk)
                if time.time() > t_stop:
                    raise _TimeUp()
            try:
                res = minimize(obj, r.ravel(), method="L-BFGS-B", jac=True, callback=cb,
                               options={"maxiter": p["maxsteps"], "maxfun": 10 * p["maxsteps"],
                                        "ftol": p["tol"], "gtol": 1e-14})
                r = res.x.reshape(N, dim)
            except _TimeUp:
                r = state["x"].reshape(N, dim)
                stop_all = True
            d_final = report(f"config {c}", r, kvecs, N, L,
                             p["Rmax"], p["lambda_h"], M, workers)
            if stop_all:
                print("  Time limit reached during L-BFGS; stopping this config.", flush=True)
                break
        configs.append(np.mod(r, L))
        diags.append(d_final)
        if stop_all:
            break

    # ---- HDF5 output ----
    out = p["savename"] + ".h5"
    with h5py.File(out, "w") as f:
        f.attrs["dim"] = dim
        f.attrs["N"] = N
        f.attrs["L"] = L
        f.attrs["n_configs"] = len(configs)
        f.attrs["Rmax"] = p["Rmax"]
        f.attrs["lambda_h"] = p["lambda_h"]
        f.attrs["n_modes"] = len(kvecs)
        f.attrs["annuli_abs"] = np.array(shells_abs, dtype=float)
        for i, (cfg, dg) in enumerate(zip(configs, diags)):
            ds = f.create_dataset(f"config_{i}", data=cfg)     # (N, dim) float64
            for key, val in dg.items():
                ds.attrs[key] = val
    print(f"\nWrote {out}  ({len(configs)} configs)", flush=True)


if __name__ == "__main__":
    main()
