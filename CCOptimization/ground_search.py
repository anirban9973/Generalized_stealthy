#!/usr/bin/env python3
"""
Ground-state search for an annular-stealthy point pattern with a real-space
hole (void) penalty in a 2D periodic square box.

Objective:
    Phi = Phi_ann + Phi_hole

    Phi_ann = sum_{k in annulus} S(k)
            = sum_{k in annulus} |rho_k|^2 / N

    Phi_hole = (lambda_h / Mgrid^2)
               sum_alpha [max(0, d_alpha - Rmax) / Rmax]^2

where d_alpha is the periodic distance from fixed grid probe alpha to its
nearest particle.

For a SINGLE annulus [K1, K2], define
    Keff = K2 - K1

and impose the 2D Zhang-style hole scale
    Rmax = 3*pi / (2*Keff).

IMPORTANT:
This use of Keff = K2-K1 is our operational extension to an annular,
nonhyperuniform stealthy constraint. The Zhang-Stillinger-Torquato result
itself was derived for a stealthy disk 0 < k < K.

The spatial grid is chosen automatically from a requested fractional
resolution epsilon_grid. For grid spacing h=L/Mgrid, require
    h/sqrt(2) <= epsilon_grid * Rmax.
Mgrid is rounded upward to a multiple of 32, with a minimum of 64.

CLI:
    python3 ground_search.py <timelimit_hr> <seed> <beg_idx> <verbosity>

stdin:
    dim  M  (K1a deltaa S0)*M  vareps0 sigma phi_fake  threads N Nc
    initmode savename mode
    [tolerance <v>] [maxsteps <v>] [algorithm <name>] run
    lambda_h epsilon_grid eps_ann eps_hole

Acceptance (BOTH must pass):  Phi_ann < eps_ann  AND  Phi_hole < eps_hole.

Example trailing input:
    run 1.0 0.05 1e-8 1e-8

Output:
    <savename>.h5
"""

import sys
import time
import math

import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
import h5py


class _TimeUp(Exception):
    """Raised from the L-BFGS callback when the wall-time margin is reached."""
    pass


# ---------------------------------------------------------------- k-space term
def phi_ann_and_grad(r, kvecs, N):
    """
    Phi_ann = sum_k |rho_k|^2 / N, with one representative from each +/-k pair.
    Returns (phi, grad) with grad shape (N,2).
    """
    kr = kvecs @ r.T
    P = np.exp(1j * kr)
    rho = P.sum(axis=1)
    phi = float(np.vdot(rho, rho).real) / N

    A = np.conj(rho)[:, None] * P
    grad = -(2.0 / N) * (kvecs.T @ A.imag).T
    return phi, grad


def structure_factor_max(r, kvecs, N):
    if len(kvecs) == 0:
        return 0.0
    kr = kvecs @ r.T
    rho = np.exp(1j * kr).sum(axis=1)
    S = np.abs(rho) ** 2 / N
    return float(S.max())


# ------------------------------------------------------------- real-space term
def make_grid(L, Mgrid):
    """Cell-centered Mgrid x Mgrid probe grid in [0,L)^2."""
    c = (np.arange(Mgrid) + 0.5) * (L / Mgrid)
    gx, gy = np.meshgrid(c, c, indexing="xy")
    return np.stack([gx.ravel(), gy.ravel()], axis=1)


def choose_grid_size(L, Rmax, epsilon_grid, block=32, min_M=64):
    """
    Choose one grid size from
        h/sqrt(2) <= epsilon_grid * Rmax,
    where h=L/Mgrid.

    Round upward to a convenient multiple of `block`.
    """
    if Rmax <= 0.0:
        raise ValueError("Rmax must be positive.")
    if epsilon_grid <= 0.0:
        raise ValueError("epsilon_grid must be positive.")

    raw = math.ceil(L / (math.sqrt(2.0) * epsilon_grid * Rmax))
    rounded = block * math.ceil(raw / block)
    return max(min_M, rounded)


def phi_hole_and_grad(r, L, grid_pts, Rmax, lambda_h, Mgrid, workers):
    """
    Grid-based void penalty:
        Phi_hole = lambda_h/Mgrid^2 *
                   sum_alpha [max(0,d_alpha-Rmax)/Rmax]^2

    d_alpha is obtained by periodic nearest-neighbor search using cKDTree.

    The nearest-particle assignment is held fixed within this evaluation.
    """
    rw = np.mod(r, L)
    tree = cKDTree(rw, boxsize=L)

    d, idx = tree.query(grid_pts, k=1, workers=workers)

    u = np.maximum(0.0, d - Rmax)
    phi = (lambda_h / (Mgrid * Mgrid)) * float(
        np.sum((u / Rmax) ** 2)
    )

    grad = np.zeros_like(r)
    active = d > Rmax

    if np.any(active):
        da = d[active]
        ida = idx[active]
        ga = grid_pts[active]

        disp = rw[ida] - ga
        disp -= L * np.round(disp / L)

        coef = (
            2.0
            * lambda_h
            / (Mgrid * Mgrid * Rmax * Rmax)
            * (da - Rmax)
            / da
        )

        np.add.at(grad, ida, coef[:, None] * disp)

    return phi, grad


def hole_diagnostics(r, L, grid_pts, Rmax, workers):
    rw = np.mod(r, L)
    tree = cKDTree(rw, boxsize=L)
    d, _ = tree.query(grid_pts, k=1, workers=workers)

    return {
        "max_d": float(d.max()),
        "frac": float(np.mean(d > Rmax)),
    }


def radial_Sk(r, N, L, kmax, chunk=40000):
    """
    Azimuthally-averaged structure factor S(|k|) on the reciprocal grid
    k = (2*pi/L) n,  0 < |k| <= kmax,  binned into shells of width dk=2*pi/L.
    Computed once per (saved) config, chunked over k to bound memory.
    Returns (k_centers, S_radial); bins with no modes are NaN.
    """
    dk = 2.0 * math.pi / L
    nmax = int(math.ceil(kmax / dk))
    ns = np.arange(-nmax, nmax + 1)
    nx, ny = np.meshgrid(ns, ns)
    n = np.stack([nx.ravel(), ny.ravel()], axis=1).astype(float)
    km = dk * np.linalg.norm(n, axis=1)
    sel = (km > 1e-12) & (km <= kmax)
    kv = dk * n[sel]
    km = km[sel]

    Ssum = np.zeros(nmax + 1)
    cnt = np.zeros(nmax + 1)
    for i in range(0, len(kv), chunk):
        kb = kv[i:i + chunk]
        rho = np.exp(1j * (kb @ r.T)).sum(axis=1)
        Sk = (np.abs(rho) ** 2) / N
        b = np.clip(np.rint(km[i:i + chunk] / dk).astype(int), 0, nmax)
        np.add.at(Ssum, b, Sk)
        np.add.at(cnt, b, 1.0)

    kc = np.arange(nmax + 1) * dk
    with np.errstate(invalid="ignore"):
        Sr = np.where(cnt > 0, Ssum / np.maximum(cnt, 1.0), np.nan)
    return kc[1:], Sr[1:]           # drop k=0 bin


# ------------------------------------------------------------------- objective
def make_objective(kvecs, N, L, grid_pts, Rmax, lambda_h, Mgrid, workers):
    shape = (N, 2)

    def obj(xflat):
        r = xflat.reshape(shape)

        phi_ann, grad_ann = phi_ann_and_grad(r, kvecs, N)
        phi_hole, grad_hole = phi_hole_and_grad(
            r, L, grid_pts, Rmax, lambda_h, Mgrid, workers
        )

        return phi_ann + phi_hole, (grad_ann + grad_hole).ravel()

    return obj


def report(
    tag,
    r,
    kvecs,
    N,
    L,
    Rmax,
    lambda_h,
    Mgrid,
    grid_pts,
    workers,
    out=sys.stdout,
):
    phi_ann, _ = phi_ann_and_grad(r, kvecs, N)
    phi_hole, _ = phi_hole_and_grad(
        r, L, grid_pts, Rmax, lambda_h, Mgrid, workers
    )

    maxS = structure_factor_max(r, kvecs, N)
    hd = hole_diagnostics(r, L, grid_pts, Rmax, workers)

    h = L / Mgrid
    grid_error_bound = h / math.sqrt(2.0)

    out.write(
        f"  [{tag}] Mgrid={Mgrid}: "
        f"Phi_ann={phi_ann:.6e}  "
        f"Phi_hole={phi_hole:.6e}  "
        f"maxS={maxS:.6e}  "
        f"max_d_grid={hd['max_d']:.6f}  "
        f"Rmax={Rmax:.6f}  "
        f"frac(d>Rmax)={hd['frac']:.6e}  "
        f"h={h:.6f}  "
        f"h/sqrt(2)={grid_error_bound:.6f}\n"
    )
    out.flush()

    return {
        "phi_ann": phi_ann,
        "phi_hole": phi_hole,
        "maxS": maxS,
        "max_d": hd["max_d"],
        "frac": hd["frac"],
        "grid_spacing": h,
        "grid_error_bound": grid_error_bound,
    }


# ----------------------------------------------------------------------- input
def parse_stdin(tokens):
    it = iter(tokens)

    def nxt():
        return next(it)

    dim = int(nxt())
    M = int(nxt())

    shells = [
        (float(nxt()), float(nxt()), float(nxt()))
        for _ in range(M)
    ]

    vareps0 = float(nxt())
    sigma = float(nxt())
    phi_fake = float(nxt())

    threads = int(nxt())
    N = int(nxt())
    Nc = int(nxt())

    initmode = nxt()
    savename = nxt()
    mode = nxt()

    tol = 1e-16
    maxsteps = 100000
    algorithm = "LBFGS"

    while True:
        t = nxt()
        if t == "run":
            break

        tl = t.lower()
        if tl == "tolerance":
            tol = float(nxt())
        elif tl == "maxsteps":
            maxsteps = int(nxt())
        elif tl == "algorithm":
            algorithm = nxt()

    tail = list(it)

    if len(tail) != 5:
        raise SystemExit(
            "After 'run', expected exactly: lambda_h epsilon_grid eps_ann eps_hole n_kmax\n"
            "Example: run 5.0 0.1 1e-12 1e-10 1000\n"
            "  eps_ann  : accept if Phi_ann  < eps_ann   (stealthiness)\n"
            "  eps_hole : accept if Phi_hole < eps_hole  (holes)\n"
            "  n_kmax   : radial S(k) is saved up to kmax = n_kmax * (2*pi/L)\n"
            "Rmax and Mgrid are still derived internally."
        )

    lambda_h = float(tail[0])
    epsilon_grid = float(tail[1])
    eps_ann = float(tail[2])
    eps_hole = float(tail[3])
    n_kmax = int(tail[4])

    return {
        "dim": dim,
        "shells": shells,
        "vareps0": vareps0,
        "sigma": sigma,
        "phi_fake": phi_fake,
        "threads": threads,
        "N": N,
        "Nc": Nc,
        "initmode": initmode,
        "savename": savename,
        "mode": mode,
        "tol": tol,
        "maxsteps": maxsteps,
        "algorithm": algorithm,
        "lambda_h": lambda_h,
        "epsilon_grid": epsilon_grid,
        "eps_ann": eps_ann,
        "eps_hole": eps_hole,
        "n_kmax": n_kmax,
    }


def unit_ball_volume(d):
    return math.pi ** (d / 2.0) / math.gamma(d / 2.0 + 1.0)


def build_kvecs(shells_abs, L, dim):
    """
    Independent reciprocal vectors k=(2*pi/L)n lying in the constrained shell(s).
    Only one of each +/-k pair is kept.
    """
    if dim != 2:
        raise SystemExit("ground_search.py currently supports dim=2 only.")

    k2max = max(k2 for (_, k2) in shells_abs)
    nmax = int(math.ceil(k2max * L / (2.0 * math.pi))) + 1

    ns = np.arange(-nmax, nmax + 1)
    nx, ny = np.meshgrid(ns, ns, indexing="xy")
    nx = nx.ravel()
    ny = ny.ravel()

    kv = (
        2.0 * math.pi / L
        * np.stack([nx, ny], axis=1).astype(float)
    )
    km = np.linalg.norm(kv, axis=1)

    in_ann = np.zeros(km.shape, dtype=bool)
    for k1, k2 in shells_abs:
        in_ann |= (km > k1 + 1e-12) & (km <= k2 + 1e-12)

    half = (nx > 0) | ((nx == 0) & (ny > 0))

    return kv[in_ann & half]


# ---------------------------------------------------------------------- driver
def main():
    timelimit_hr = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    verbosity = int(sys.argv[4]) if len(sys.argv) > 4 else 3

    p = parse_stdin(sys.stdin.read().split())

    dim = p["dim"]
    N = p["N"]

    if dim != 2:
        raise SystemExit(
            "Automatic Rmax = 3*pi/[2*(K2-K1)] is implemented here for d=2 only."
        )

    L = math.sqrt(N)  # unit number density
    workers = max(1, p["threads"])
    t_stop = time.time() + timelimit_hr * 3600.0 - 300.0

    a = (
        p["phi_fake"] / unit_ball_volume(dim)
    ) ** (1.0 / dim)

    shells_abs = []
    n_equi = 0

    for K1a, deltaa, S0 in p["shells"]:
        if S0 != 0.0:
            n_equi += 1
            continue

        K1 = K1a / a
        K2 = (K1a + deltaa) / a
        shells_abs.append((K1, K2))

    if n_equi:
        print(
            f"WARNING: {n_equi} equiluminous shell(s) (S0>0) ignored; "
            "this version is stealthy-only.",
            flush=True,
        )

    if not shells_abs:
        raise SystemExit("No stealthy (S0=0) shell to constrain.")

    if len(shells_abs) != 1:
        raise SystemExit(
            "Automatic Zhang-style Rmax currently requires exactly ONE "
            "stealthy annulus [K1,K2]."
        )

    K1, K2 = shells_abs[0]
    Keff = K2 - K1

    if Keff <= 0.0:
        raise SystemExit("Need K2 > K1.")

    Rmax = 3.0 * math.pi / (2.0 * Keff)

    Mgrid = choose_grid_size(
        L=L,
        Rmax=Rmax,
        epsilon_grid=p["epsilon_grid"],
        block=32,
        min_M=64,
    )

    grid_pts = make_grid(L, Mgrid)
    h = L / Mgrid
    grid_error_bound = h / math.sqrt(2.0)

    kvecs = build_kvecs(shells_abs, L, dim)
    chi = len(kvecs) / (dim * (N - 1))

    print(
        f"ground_search: N={N}, L={L:.6f} (rho=1), "
        f"dim={dim}, threads={workers}",
        flush=True,
    )
    print(
        f"  annulus: K1={K1:.8f}, K2={K2:.8f}, "
        f"Keff=K2-K1={Keff:.8f}",
        flush=True,
    )
    print(
        f"  constrained independent modes: {len(kvecs)}   "
        f"chi={chi:.8f}",
        flush=True,
    )
    print(
        f"  automatic Rmax = 3*pi/(2*Keff) = {Rmax:.8f}",
        flush=True,
    )
    print(f"  lambda_h={p['lambda_h']}", flush=True)
    print(
        f"  epsilon_grid={p['epsilon_grid']}  "
        f"Mgrid={Mgrid}  h={h:.8f}  "
        f"h/sqrt(2)={grid_error_bound:.8f}  "
        f"(target <= epsilon_grid*Rmax={p['epsilon_grid'] * Rmax:.8f})",
        flush=True,
    )
    print(
        f"  Nc={p['Nc']}, maxsteps={p['maxsteps']}, "
        f"ftol={p['tol']}\n",
        flush=True,
    )

    if p["algorithm"].upper() not in ("LBFGS", "L-BFGS", "L-BFGS-B"):
        print(
            f"WARNING: algorithm '{p['algorithm']}' not supported; "
            "using L-BFGS-B.",
            flush=True,
        )

    kmax = p["n_kmax"] * (2.0 * math.pi / L)
    print(f"  radial S(k) saved up to kmax={kmax:.4f} "
          f"(n_kmax={p['n_kmax']} shells of 2*pi/L)\n", flush=True)

    rng = np.random.default_rng(seed)

    configs = []
    diags = []
    sk_list = []
    k_radial = None
    stop_all = False

    for c in range(p["Nc"]):
        if time.time() > t_stop:
            print(
                f"Time limit reached; produced {c}/{p['Nc']} configs.",
                flush=True,
            )
            break

        print(f"--- config {c} ---", flush=True)
        r = rng.random((N, dim)) * L

        obj = make_objective(
            kvecs=kvecs,
            N=N,
            L=L,
            grid_pts=grid_pts,
            Rmax=Rmax,
            lambda_h=p["lambda_h"],
            Mgrid=Mgrid,
            workers=workers,
        )

        state = {"x": r.ravel().copy()}

        def callback(xk, _state=state):
            _state["x"] = np.asarray(xk).copy()
            if time.time() > t_stop:
                raise _TimeUp()

        try:
            res = minimize(
                obj,
                r.ravel(),
                method="L-BFGS-B",
                jac=True,
                callback=callback,
                options={
                    "maxiter": p["maxsteps"],
                    "maxfun": 10 * p["maxsteps"],
                    "ftol": p["tol"],
                    "gtol": 1e-14,
                },
            )

            r = res.x.reshape(N, dim)

            if verbosity >= 2:
                print(
                    f"  optimizer success={res.success} "
                    f"status={res.status} nit={res.nit} "
                    f"nfev={res.nfev} message={res.message}",
                    flush=True,
                )

        except _TimeUp:
            r = state["x"].reshape(N, dim)
            stop_all = True
            print(
                "  Time limit reached during L-BFGS; saving latest iterate.",
                flush=True,
            )

        d_final = report(
            tag=f"config {c}",
            r=r,
            kvecs=kvecs,
            N=N,
            L=L,
            Rmax=Rmax,
            lambda_h=p["lambda_h"],
            Mgrid=Mgrid,
            grid_pts=grid_pts,
            workers=workers,
        )

        # Acceptance: BOTH thresholds must pass.
        pass_ann = d_final["phi_ann"] < p["eps_ann"]
        pass_hole = d_final["phi_hole"] < p["eps_hole"]
        accepted = pass_ann and pass_hole
        d_final["pass_ann"] = int(pass_ann)
        d_final["pass_hole"] = int(pass_hole)
        d_final["accepted"] = int(accepted)
        print(
            f"  ACCEPT: STEALTHY {'PASS' if pass_ann else 'FAIL'} "
            f"(Phi_ann={d_final['phi_ann']:.3e} < {p['eps_ann']:.1e})  |  "
            f"HOLE {'PASS' if pass_hole else 'FAIL'} "
            f"(Phi_hole={d_final['phi_hole']:.3e} < {p['eps_hole']:.1e})  "
            f"=>  {'ACCEPTED' if accepted else 'REJECTED'}",
            flush=True,
        )

        # radial S(k) of the final config (saved whether accepted or not)
        kc, Sr = radial_Sk(r, N, L, kmax)
        if k_radial is None:
            k_radial = kc
        sk_list.append(Sr)

        configs.append(np.mod(r, L))
        diags.append(d_final)

        if stop_all:
            break

    out = p["savename"] + ".h5"

    with h5py.File(out, "w") as f:
        f.attrs["dim"] = dim
        f.attrs["N"] = N
        f.attrs["L"] = L
        f.attrs["rho"] = N / (L * L)
        f.attrs["n_configs"] = len(configs)

        f.attrs["K1"] = K1
        f.attrs["K2"] = K2
        f.attrs["Keff"] = Keff
        f.attrs["n_modes"] = len(kvecs)
        f.attrs["chi"] = chi
        f.attrs["annuli_abs"] = np.array(shells_abs, dtype=float)

        f.attrs["Rmax"] = Rmax
        f.attrs["Rmax_definition"] = "3*pi/[2*(K2-K1)]"

        f.attrs["lambda_h"] = p["lambda_h"]
        f.attrs["epsilon_grid"] = p["epsilon_grid"]
        f.attrs["eps_ann"] = p["eps_ann"]
        f.attrs["eps_hole"] = p["eps_hole"]
        f.attrs["n_kmax"] = p["n_kmax"]
        f.attrs["kmax"] = kmax

        # radial S(k): one shared |k| axis + one S(k) per config
        if k_radial is not None:
            f.create_dataset("k_radial", data=k_radial)
            for i, Sr in enumerate(sk_list):
                f.create_dataset(f"Sk_radial_{i}", data=Sr)
        f.attrs["Mgrid"] = Mgrid
        f.attrs["grid_spacing"] = h
        f.attrs["grid_error_bound"] = grid_error_bound

        for i, (cfg, dg) in enumerate(zip(configs, diags)):
            ds = f.create_dataset(f"config_{i}", data=cfg)
            for key, val in dg.items():
                ds.attrs[key] = val

    print(f"\nWrote {out}  ({len(configs)} configs)", flush=True)


if __name__ == "__main__":
    main()
