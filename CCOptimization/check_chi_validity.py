"""
Scan crystalmd SLURM .out logs and report which chi runs are valid (stable MD)
vs blown-up, writing a summary to chi_validity.txt.

Classification (per array task), using the sampling-stage kinetic energy E_k
compared to the equipartition value Ek_eq = 0.5*d*N*T_E  (d=2), e.g. 625*1e-6:
    OK     : MD stayed stable AND the final E_k sits at the target temperature
             (final E_k <= HOT_FACTOR * Ek_eq).  -> usable.
    HOT    : MD did NOT explode, but the final E_k is still > HOT_FACTOR * Ek_eq,
             i.e. it never cooled/equilibrated to T_E.  -> configs too hot, not usable as-is.
    BLEWUP : the integrator went unstable and E_k exploded
             (max sampling E_k > BLOWUP_FACTOR * Ek_eq).  -> garbage.
    NO_SAMPLES : the log has no sampling lines (job crashed/ended before sampling).
    "bad" in the summary = BLEWUP + NO_SAMPLES.

Only *.out files are read (SLURM .err files are ignored).

Usage:
    python3 check_chi_validity.py [root_dir]      (default: current directory)
Finds every *CCO_crystalmd_*.out under root_dir (e.g. chi_*/logs/).
"""

import sys
import os
import glob
import re

# ---- thresholds (tunable) ----
BLOWUP_FACTOR = 1000.0   # max sampling E_k above this * Ek_eq  -> exploded (BLEWUP)
HOT_FACTOR    = 3.0      # final E_k above this * Ek_eq         -> not equilibrated (HOT)
DIM           = 2

root = sys.argv[1] if len(sys.argv) > 1 else "."

re_chi   = re.compile(r"chi \(target\)\s*:\s*([-\d.eE+]+)")
re_task  = re.compile(r"Array task ID\s*:\s*(\d+)")
re_total = re.compile(r"total\s*=\s*(\d+)")
re_TE    = re.compile(r"T_E\s*:\s*([-\d.eE+]+)")
re_samp  = re.compile(r"\b2:(\d+)/\d+.*?E_relax=([-\d.eE+]+).*?E_k=([-\d.eE+]+)")

# only .out files (ignore .err)
files = sorted(glob.glob(os.path.join(root, "**", "*CCO_crystalmd_*.out"), recursive=True))
files = [f for f in files if os.path.isfile(f)]
if not files:
    raise SystemExit(f"No *CCO_crystalmd_*.out logs found under {root!r}.")

rows = []   # (chi, task, n_samp, max_Ek, final_Ek, Ek_eq, status)
for path in files:
    chi = task = total = T_E = None
    eks = []
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

def fmt_tasks(ts):
    return ",".join(str(t) for t in ts) if ts else "-"

# ---- summary per chi ----
chis = sorted(set(r[0] for r in rows))
summary = []
for chi in chis:
    sub = [r for r in rows if r[0] == chi]
    n = len(sub)
    ok_t  = sorted(r[1] for r in sub if r[6] == "OK")
    hot_t = sorted(r[1] for r in sub if r[6] == "HOT")
    bad_t = sorted(r[1] for r in sub if r[6] in ("BLEWUP", "NO_SAMPLES"))
    verdict = "VALID" if len(ok_t) == n else ("PARTIAL" if len(ok_t) > 0 else "INVALID")
    summary.append((chi, n, ok_t, hot_t, bad_t, verdict))

legend = (
    "# Legend:\n"
    "#   OK     = stable MD, final E_k at target temperature (usable)\n"
    "#   HOT    = did NOT explode but never cooled to T_E (final E_k > "
    f"{HOT_FACTOR:g}x Ek_eq); too hot, not usable as-is\n"
    "#   BLEWUP = integrator exploded (max E_k > "
    f"{BLOWUP_FACTOR:g}x Ek_eq); garbage\n"
    "#   bad    = BLEWUP + NO_SAMPLES\n"
)

# ---- write ----
out = "chi_validity.txt"
with open(out, "w") as f:
    f.write("# Crystal MD validity check\n")
    f.write(f"# Ek_eq = 0.5*d*N*T_E (d={DIM})\n")
    f.write(legend + "#\n")
    f.write(f"# {'chi':>6} {'task':>4} {'n_samp':>7} {'max_Ek':>12} {'final_Ek':>12} {'Ek_eq':>12}  status\n")
    for chi, task, ns, mx, fn, eq, st in rows:
        f.write(f"  {chi:>6.4g} {task:>4d} {ns:>7d} {mx:>12.4g} {fn:>12.4g} {eq:>12.4g}  {st}\n")
    f.write("#\n# per-chi summary (array task numbers listed by status):\n")
    f.write(f"# {'chi':>6} {'n':>3} {'OK':>3} {'HOT':>3} {'bad':>3}  {'verdict':<8} "
            f"{'OK_tasks':<28} {'HOT_tasks':<14} bad_tasks\n")
    for chi, n, ok_t, hot_t, bad_t, verdict in summary:
        f.write(f"  {chi:>6.4g} {n:>3d} {len(ok_t):>3d} {len(hot_t):>3d} {len(bad_t):>3d}  "
                f"{verdict:<8} {fmt_tasks(ok_t):<28} {fmt_tasks(hot_t):<14} {fmt_tasks(bad_t)}\n")

# ---- print ----
print(legend)
print("Summary (array task numbers by status):")
print(f"  {'chi':>6} {'n':>3} {'OK':>3} {'HOT':>3} {'bad':>3}  {'verdict':<8}")
for chi, n, ok_t, hot_t, bad_t, verdict in summary:
    print(f"  {chi:>6.4g} {n:>3d} {len(ok_t):>3d} {len(hot_t):>3d} {len(bad_t):>3d}  {verdict}")
    print(f"        OK  : {fmt_tasks(ok_t)}")
    print(f"        HOT : {fmt_tasks(hot_t)}")
    print(f"        bad : {fmt_tasks(bad_t)}")
print(f"\nWrote {out}")
