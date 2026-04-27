# Generalized Stealthy Ground States — Branch `generalized_stealthy`

## Overview

This branch extends the stealthy collective-coordinate optimization to enforce S(**k**) = 0
over a **union of M annular shells** in k-space:

$$S(\mathbf{k}) = 0 \quad \text{for } |\mathbf{k}| \in \bigcup_{n=1}^{M} [k_1^{(n)},\, k_1^{(n)} + \delta^{(n)}]$$

Previously, only a single shell [K1, K2] was supported. The generalization allows
arbitrary non-overlapping excluded regions anywhere in k-space.


## Parameterization

Each shell n is defined by two parameters (both in unit-number-density units):

| Parameter    | Meaning                                      |
|--------------|----------------------------------------------|
| `K1s[n]`     | lower bound k₁⁽ⁿ⁾ of shell n                |
| `deltas[n]`  | width δ⁽ⁿ⁾ of shell n; upper bound is k₂⁽ⁿ⁾ = k₁⁽ⁿ⁾ + δ⁽ⁿ⁾ |

The former `chi` parameter (which implicitly set K2 via K2 = 2π(2dχ/Vd)^(1/d)) is
**no longer used**. Shell widths are now set directly via `deltas`.

The `phi_fake` encoding is preserved: the script multiplies K1 and delta by
`a = (phi_fake / Vd)^(1/d)` before passing to the CLI, and the code divides by `a`
to recover absolute k-space values.


## Changes to `ex_ground.sh`

### Removed
```bash
K1=0.0   # single lower bound
chi=0.3  # stealthiness parameter
```

### Added
```bash
M=2                   # number of stealthy shells
K1s=(0.0  0.5)        # lower bounds k1^(n), one per shell (unit number density)
deltas=(0.3  0.2)     # widths delta^(n), one per shell (unit number density)
```

### Changed: pre-processing block
The script now loops over all M shells to build the encoded shell string passed to the CLI:

```bash
shell_str="${M}"
for (( n=0; n<M; n++ )); do
    K1a=`echo "${K1s[$n]}*${a}"     | bc -l`
    deltaa=`echo "${deltas[$n]}*${a}" | bc -l`
    shell_str="${shell_str} ${K1a} ${deltaa}"
done
```

The executable is then called as:
```bash
${exe} ${timelimit} ${seed} ${beg_idx} ${Verbosity} \
    <<< "${d} ${shell_str} $S0 $vareps0 ${sigma} ${phi_fake} \
         ${threads} ${N} ${Nc} random ${fname}_GS \
         ground $eps0 $max_eval $algorithm run"
```

### Changed: soft-core repulsion default
`vareps0` is set to `0.0` (no soft-core repulsion) by default in this branch.


## Changes to `main.cpp`

### Variable declarations (line ~164)

**Removed:**
```cpp
double K1, K2;
```

**Added:**
```cpp
std::vector<std::pair<double,double>> shells; // (K1a, deltaa) per shell
```

### Input reading (line ~171)

**Before:**
```cpp
ofile << "K1a = "; ifile >> K1;
ofile << "K2a = "; ifile >> K2;
```

**After:**
```cpp
size_t M = 1;
ofile << "M (number of stealthy shells) = "; ifile >> M;
for (size_t n = 0; n < M; n++) {
    double K1n, deltan;
    ofile << "K1a[" << n << "] = ";    ifile >> K1n;
    ofile << "delta_a[" << n << "] = "; ifile >> deltan;
    shells.emplace_back(K1n, deltan);
}
```

### New CLI input format

```
dim  M  K1a_0 deltaa_0  K1a_1 deltaa_1  ...  S0 val sigma phi  threads N Nc  initial savename  mode ...
```

where `K1a_n = K1^(n) * a` and `deltaa_n = delta^(n) * a` are the encoded values.

### Potential setup — unit conversion (line ~234)

**Before:**
```cpp
double k1 = K1 / a, k2 = K2 / a;
```

**After:**
```cpp
std::vector<std::pair<double,double>> shells_abs; // (k1, k2) in absolute units
double k2_max = 0;
for (size_t n = 0; n < M; n++) {
    double k1n = shells[n].first  / a;
    double k2n = k1n + shells[n].second / a;   // k2 = k1 + delta
    shells_abs.emplace_back(k1n, k2n);
    k2_max = std::max(k2_max, k2n);
}
```

### Constraint-adding loop (line ~265)

**Before** — single shell filter:
```cpp
std::vector<GeometryVector> ks_temp = GetKs(Config, k2, k2, 1);
for (auto k = ks_temp.begin(); k != ks_temp.end(); k++) {
    if (k->Modulus2() > K1_modulus) {
        potential->CCPotential->AddConstraint(*k, vals);
    }
}
```

**After** — union of M shells:
```cpp
std::vector<GeometryVector> ks_temp = GetKs(Config, k2_max, k2_max, 1);
for (auto k = ks_temp.begin(); k != ks_temp.end(); k++) {
    double km2 = k->Modulus2();
    for (auto& s : shells_abs) {
        if (km2 > s.first * s.first && km2 <= s.second * s.second) {
            potential->CCPotential->AddConstraint(*k, vals);
            chi++;
            break;  // each k belongs to at most one shell
        }
    }
}
```

The same pattern is applied to both the `S0 == 0` (stealthy) and `S0 > 0` (equiluminous) branches.

The `vtilde` normalization reference is changed from `k2` to `k2_max` for the `overlap`
and `power-law` options, since there is no longer a single outer cutoff.


## Files Not Changed

- `Potential/CollectiveCoordinatePotential.h` — `ShiftedCCPotential` and its `Energy()`,
  `Force()`, `AddConstraint()` are unchanged; they operate on whatever constraints are
  in the vector regardless of how they were assembled.
- `CCOptimization/RepulsiveCCPotential.h` — unchanged.
- All other modules (`PairStat`, `cores`, `Potential`, etc.) — unchanged.


## Backward Compatibility

A single standard stealthy shell (original behavior) is recovered by setting:
```bash
M=1
K1s=(0.0)
deltas=(<your_delta>)
```
where `delta = K2 - K1` from the old parameterization.
