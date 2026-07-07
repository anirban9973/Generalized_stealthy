"""
Quick look at the first configuration in a configs*.h5 file, drawn inside the
rhombic simulation box.

Auto-loads the first file in the current directory matching configs*.h5
(e.g. configs_N2500_chi0.55.h5), plots config_0, and saves a PNG.
"""

import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt

# find a configs*.h5 file
files = sorted(glob.glob("configs*.h5"))
if not files:
    raise SystemExit("No configs*.h5 file found in the current directory.")
h5path = files[0]
print(f"Loading {h5path}")

with h5py.File(h5path, "r") as f:
    basis = np.array(f.attrs["basis"])   # (2, 2): rows are box vectors a1, a2
    N     = int(f.attrs["N"])
    chi   = float(f.attrs["chi"]) if "chi" in f.attrs else None
    pos   = f["config_0"][:]             # (N, 2) Cartesian coordinates

a1, a2 = basis[0], basis[1]

# rhombic box outline: O -> a1 -> a1+a2 -> a2 -> O
box = np.array([[0, 0], a1, a1 + a2, a2, [0, 0]])

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(box[:, 0], box[:, 1], "k-", linewidth=1.0)
ax.scatter(pos[:, 0], pos[:, 1], s=6, c="steelblue", linewidths=0)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
title = f"config_0   (N={N}" + (f", chi={chi:g}" if chi is not None else "") + ")"
ax.set_title(title)

plt.tight_layout()
outname = h5path.replace(".h5", "_config0.png")
plt.savefig(outname, dpi=150)
print(f"Saved {outname}")
plt.show()
