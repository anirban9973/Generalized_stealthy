"""
Overlay two configurations from a configs*.h5 file in the rhombic box, each in a
different color with large points, to compare them by eye.

Auto-loads the first configs*.h5 in the current directory and overlays configs
i and j (set below). Saves a PNG.
"""

import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt

# ---- which two configurations to overlay ----
i, j = 0, 1
# ---------------------------------------------

files = sorted(glob.glob("configs*.h5"))
if not files:
    raise SystemExit("No configs*.h5 file found in the current directory.")
h5path = files[0]
print(f"Loading {h5path}")

with h5py.File(h5path, "r") as f:
    basis = np.array(f.attrs["basis"])   # (2, 2): rows are box vectors a1, a2
    N     = int(f.attrs["N"])
    chi   = float(f.attrs["chi"]) if "chi" in f.attrs else None
    n     = int(f.attrs["n_configs"])
    if i >= n or j >= n:
        raise SystemExit(f"File has only {n} config(s); cannot load {i} and {j}.")
    pos_i = f[f"config_{i}"][:]           # (N, 2)
    pos_j = f[f"config_{j}"][:]           # (N, 2)

a1, a2 = basis[0], basis[1]
box = np.array([[0, 0], a1, a1 + a2, a2, [0, 0]])

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(box[:, 0], box[:, 1], "k-", linewidth=1.0)
ax.scatter(pos_i[:, 0], pos_i[:, 1], s=60, c="crimson",   alpha=0.6,
           linewidths=0, label=f"config_{i}")
ax.scatter(pos_j[:, 0], pos_j[:, 1], s=60, c="royalblue", alpha=0.6,
           linewidths=0, label=f"config_{j}")
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
title = f"config_{i} vs config_{j}   (N={N}" + (f", chi={chi:g}" if chi is not None else "") + ")"
ax.set_title(title)
ax.legend(loc="upper right")

plt.tight_layout()
outname = h5path.replace(".h5", f"_overlay_{i}_{j}.png")
plt.savefig(outname, dpi=150)
print(f"Saved {outname}")
plt.show()
