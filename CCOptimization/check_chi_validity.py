"""
Scan crystalmd SLURM .out logs and report which chi runs are valid (stable MD)
vs blown-up, writing a summary to chi_validity.txt.

A run's MD is "valid" if its sampling-stage kinetic energy stays near the
equipartition value and never explodes:
    Ek_eq = 0.5 * d * N * T_E        (d = 2)   e.g. 625 * 1e-6 = 6.25e-4
Per task we look at the sampling lines ("... 2:i/nsamp ... E_k=..."):
    BLEWUP : max Ek over sampling  > BLOWUP_FACTOR * Ek_eq   (integrator exploded)
    HOT    : not blown, but final Ek > HOT_FACTOR * Ek_eq    (never settled to T_E)
    OK     : final Ek within HOT_FACTOR of Ek_eq             (clean, at temperature)

Usage:
    python3 check_chi_validity.py [root_dir]      (default: current directory)
Finds all files whose name contains 'CCO_crystalmd_' anywhere under root_dir.
"""

import sys
import os
import glob
import re

# ---- thresholds (tunable) ----
BLOWUP_FACTOR = 1000.0   # max sampling Ek above this * Ek_eq  -> exploded
HOT_FACTOR    = 3.0      # final Ek above this * Ek_eq         -> not equilibrated
DIM           = 2

root = sys.argv[1] if len(sys.argv) > 1 else "."

re_chi   = re.compile(r"chi \(target\)\s*:\s*([-\d.eE+]+)")
re_task  = re.compile(r"Array task ID\s*:\s*(\d+)")
re_total = re.compile(r"total\s*=\s*(\d+)")
re_TE    = re.compile(r"T_E\s*:\s*([-\d.eE+]+)")
re_samp  = re.compile(r"\b2:(\d+)/\d+.*?E_relax=([-\d.eE+]+).*?E_k=([-\d.eE+]+)")

files = sorted(glob.glob(os.path.join(root, "**", "*CCO_crystalmd_*"), recursive=True))
files = [f for f in files if os.path.isfile(f)]
if not files:
    raise SystemExit(f"No CCO_crystalmd_* logs found under {root!r}.")

rows = []   # (chi, task, n_samp, max_Ek, final_Ek, Ek_eq, status)
for path in files:
    chi = task = total = T_E = None
    eks = []
    final_ep = None
    with open(path, errors="replace") as fh:
        for line in fh:
            if chi is None:
                m = re_chi.search(line);   chi   = float(m.group(1)) if m else None
            if task is None:
                m = re_task.search(line);  task  = int(m.group(1))   if m else None
            if total is None:
                m = re_total.search(line); total = int(m.group(1))   if m else None
            if T_E is None:
                m = re_TE.search(line);    T_E   = float(m.group(1)) if m else None
            m = re_samp.search(line)
            if m:
                eks.append(float(m.group(3)))
                final_ep = float(m.group(2))

    if chi is None or total is None or T_E is None:
        print(f"  skip (missing header): {path}")
        continue

    Ek_eq = 0.5 * DIM * total * T_E
    if not eks:
        status, max_ek, final_ek = "NO_SAMPLES", float("nan"), float("nan")
    else:
        max_ek, final_ek = max(eks), eks[-1]
        if max_ek > BLOWUP_FACTOR * Ek_eq:
            status = "BLEWUP"
        elif final_ek > HOT_FACTOR * Ek_eq:
            status = "HOT"
        else:
            status = "OK"
    rows.append((chi, task if task is not None else -1,
                 len(eks), max_ek, final_ek, Ek_eq, status))

rows.sort(key=lambda r: (r[0], r[1]))

# ---- summary per chi ----
chis = sorted(set(r[0] for r in rows))
summary = []
for chi in chis:
    sub = [r for r in rows if r[0] == chi]
    n = len(sub)
    n_ok  = sum(1 for r in sub if r[6] == "OK")
    n_hot = sum(1 for r in sub if r[6] == "HOT")
    n_bad = sum(1 for r in sub if r[6] in ("BLEWUP", "NO_SAMPLES"))
    verdict = "VALID" if n_ok == n else ("PARTIAL" if n_ok > 0 else "INVALID")
    summary.append((chi, n, n_ok, n_hot, n_bad, verdict))

# ---- write + print ----
out = "chi_validity.txt"
with open(out, "w") as f:
    f.write("# Crystal MD validity check\n")
    f.write(f"# Ek_eq = 0.5*d*N*T_E (d={DIM}); BLEWUP if max Ek > {BLOWUP_FACTOR:g}*Ek_eq; "
            f"HOT if final Ek > {HOT_FACTOR:g}*Ek_eq\n#\n")
    f.write(f"# {'chi':>6} {'task':>4} {'n_samp':>7} {'max_Ek':>12} {'final_Ek':>12} {'Ek_eq':>12}  status\n")
    for chi, task, ns, mx, fn, eq, st in rows:
        f.write(f"  {chi:>6.4g} {task:>4d} {ns:>7d} {mx:>12.4g} {fn:>12.4g} {eq:>12.4g}  {st}\n")
    f.write("#\n# per-chi summary:\n")
    f.write(f"# {'chi':>6} {'n_tasks':>7} {'n_OK':>5} {'n_HOT':>5} {'n_bad':>5}  verdict\n")
    for chi, n, n_ok, n_hot, n_bad, verdict in summary:
        f.write(f"  {chi:>6.4g} {n:>7d} {n_ok:>5d} {n_hot:>5d} {n_bad:>5d}  {verdict}\n")

print(f"Scanned {len(rows)} log(s). Summary:")
print(f"  {'chi':>6} {'n_tasks':>7} {'n_OK':>5} {'n_HOT':>5} {'n_bad':>5}  verdict")
for chi, n, n_ok, n_hot, n_bad, verdict in summary:
    print(f"  {chi:>6.4g} {n:>7d} {n_ok:>5d} {n_hot:>5d} {n_bad:>5d}  {verdict}")
print(f"\nWrote {out}")
