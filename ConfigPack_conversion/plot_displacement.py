"""
Visualize how two configurations differ, per particle, as a displacement map.

Particle i is the same lattice site in every config (all start from the same
perfect lattice), so we draw an arrow from its position in config i to its
position in config j. Displacements use the minimum-image convention; the mean
(rigid translation) is removed by default so the interesting structure (phonons,
defects) stands out, and arrows are magnified so tiny shifts are visible.

Auto-loads the first configs*.h5 in the current directory.
"""

import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt

# ---- settings ----
i, j        = 0, 1     # which two configs to compare
magnify     = 20.0     # arrow length multiplier (tiny displacements -> visible)
remove_mean = True     # subtract the global translation before plotting
# ------------------

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
        raise SystemExit(f"File has only {n} config(s); cannot compare {i} and {j}.")
    pos_i = f[f"config_{i}"][:]           # (N, 2)
    pos_j = f[f"config_{j}"][:]           # (N, 2)

a1, a2 = basis[0], basis[1]
A = np.column_stack([a1, a2])            # columns are box vectors; r = A @ rel
Ainv = np.linalg.inv(A)

# minimum-image displacement: wrap the relative-coordinate difference to [-0.5, 0.5]
drel = (Ainv @ (pos_j - pos_i).T).T      # (N, 2) in relative coords
drel -= np.round(drel)
disp = (A @ drel.T).T                     # back to Cartesian (N, 2)

if remove_mean:
    disp = disp - disp.mean(axis=0)

mag = np.linalg.norm(disp, axis=1)
print(f"displacement (min-image{' , mean removed' if remove_mean else ''}): "
      f"mean={mag.mean():.4g}, max={mag.max():.4g}")

# rhombic box outline
box = np.array([[0, 0], a1, a1 + a2, a2, [0, 0]])

fig, ax = plt.subplots(figsize=(6.5, 6))
ax.plot(box[:, 0], box[:, 1], "k-", linewidth=1.0)
q = ax.quiver(pos_i[:, 0], pos_i[:, 1],
              disp[:, 0] * magnify, disp[:, 1] * magnify,
              mag, cmap="viridis",
              angles="xy", scale_units="xy", scale=1, width=0.004)
fig.colorbar(q, ax=ax, label="|displacement|")
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
title = (f"config_{i} -> config_{j}  (arrows x{magnify:g}"
         + (", mean removed" if remove_mean else "") + f")\nN={N}"
         + (f", chi={chi:g}" if chi is not None else ""))
ax.set_title(title)

plt.tight_layout()
outname = h5path.replace(".h5", f"_disp_{i}_{j}.png")
plt.savefig(outname, dpi=150)
print(f"Saved {outname}")
plt.show()
