"""
Compute the stealthiness parameter chi for M annular shells.

    chi = sum_{n=1}^{M}  (k_{n,2}^d - k_{n,1}^d)
                         --------------------------------
                         2 Gamma(1 + d/2) * 2^d * pi^{d/2} * d

Evaluated at unit number density (rho = 1).
Chi counts constrained k-modes and is independent of S0.
"""

import numpy as np
from scipy.special import gamma
import math


def compute_chi(d, shells):
    """
    Parameters
    ----------
    d      : int
        Spatial dimension.
    shells : list of (k1, k2) tuples
        Inner and outer radii of each annular shell in absolute k units
        (at unit number density). S0 value is irrelevant.

    Returns
    -------
    chi : float
    """
    denom = 2 * gamma(1 + d / 2) * (2**d) * (np.pi ** (d / 2)) * d
    return sum((k2**d - k1**d) / denom for k1, k2 in shells)


if __name__ == "__main__":
    # ---------------------------------------------------------------
    # Edit these parameters to match your simulation setup
    # ---------------------------------------------------------------
    d = 2
    delta = 2.481
    N = 4000
    L = np.sqrt(N)
    p = 25
    k1 = p * (2 * math.pi / L)
    k2 = k1 + delta

    # One (k1_n, k2_n) tuple per shell — same values as K1s and K1s+deltas
    # in ex_ground.sh at unit number density.
    # Example: two shells, both stealthy (S0 irrelevant here)
    shells = [
        (0.0, k1),   # shell 0: 
        (k1, k2),   # shell 1:
    ]
    # ---------------------------------------------------------------

    chi = compute_chi(d, shells)

    print(f"d       = {d}")
    print(f"M       = {len(shells)}")
    for n, (k1, k2) in enumerate(shells):
        print(f"  shell {n}: k1={k1}, k2={k2}")
    print(f"chi     = {chi:.6f}")
