# Geometric Sieve

A structural framework for prime distribution based on square layers, plus a
representation function that reveals a striking asymptotic constant.

**Core result:** every integer n lies in a layer `L_a = [a², (a+1)²)` with
geometric center `C(a) = a(a+1)`. The condition `gcd(|b|, C(a)) > 1`, where
`b = C(a) - n`, deterministically flags composite numbers with **zero false
positives**, reaching ~73.5% efficiency by 10⁶.

Separately, for primes `p`, define `R(p)` as the number of primes `a < p`
such that `b = 2a - p` is 1 or prime. Across 664,578 primes up to 10⁷:

R(p) · (ln p)² / p → 0.73 (std. dev. ≈ 0.009)

---

## Contents

- [Overview](#overview)
- [Main Results](#main-results)
- [Exploratory Extensions](#exploratory-extensions)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Usage](#usage)
- [Data](#data)
- [Status and Limitations](#status-and-limitations)
- [License](#license)
- [Citation](#citation)

---

## Overview

This repository documents two related but distinct layers of work:

1. A **proved geometric sieve** and a **numerically well-supported
   representation constant**, both described in
   [`paper_geometric_sieve.pdf`](./paper_geometric_sieve.pdf).
2. A set of **exploratory numerical extensions** (mod-4 prime classes,
   square-band symmetry, a unified comparison score) that are still at the
   observational stage and are not yet backed by proofs.

The two layers are kept separate on purpose. The first is publication-ready;
the second is active research documented here for transparency and
reproducibility.

---

## Main Results

These results are proved or numerically established in the paper.

| Result | Statement | Status |
|---|---|---|
| **Geometric Sieve** | `gcd(\|b\|, C(a)) > 1` ⟹ n composite, with no false positives | Proved (Theorem 1) |
| **Sieve efficiency** | ~73.5% of composites caught by 10⁶ | Numerically verified |
| **Representation constant** | `R(p)·(ln p)²/p → C ≈ 0.73`, std. dev. 0.009 over 664,578 primes to 10⁷ | Numerically established, unproved asymptotic |
| **Convergence** | Constant drops from 0.774 (p < 10⁴) to 0.725 (p > 10⁶) | Observed trend |

See `paper_geometric_sieve.pdf` for full definitions, proofs, and the
heuristic derivation of the 0.73 constant.

---

## Exploratory Extensions

The following observations come from later, independent scans that use a
different geometric frame (square **bands** between consecutive squares,
rather than the layer/center system above). They are included for
completeness but should be read as **preliminary numerical patterns**, not
theorems:

| Observation | Description | Status |
|---|---|---|
| **Isolation** | Primes show fewer `R`-type representations than composites (~1.69 vs. ~2.79 mean `sum_norm`) | Numerical, multi-scale |
| **Band symmetry** | Primes show near left-right symmetry within square bands | Numerical, multi-scale |
| **4n+1 vs 4n-1 offset** | The two residue classes show a small, consistently signed gap in representation ratio inside square bands | Numerical, stable across N = 25k–250k, small effect size |
| **Unified score** | Aggregated gap between 4n+1 and 4n-1 ≈ 0.048 across tested dimensions | Descriptive statistic, not a proven constant |

These extensions are not yet integrated into the main paper's claims and
should not be cited as established results. They are tracked here as a
future-work line.

---

## Repository Structure

```
├── README.md
├── requirements.txt
├── LICENSE
├── paper_geometric_sieve.pdf
├── code/
│   ├── verify_sieve.py
│   ├── verify_Rp.py
│   ├── singular_series.py
│   ├── odd_vs_prime_scan.py
│   ├── odd_vs_prime_band.py
│   └── mod4_unified_analysis.py
├── data/
│   ├── R_p_results_10M.csv
│   └── scale_summary.csv
└── figures/
    └── (generated plots)
```

---

## Requirements

```bash
pip install -r requirements.txt
```

`requirements.txt`:

```
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
```

---

## Usage

### 1. Verify the geometric sieve

```bash
python code/verify_sieve.py
```

Reproduces Theorem 1 and the efficiency table (73.5% at 10⁶) with zero false
positives.

### 2. Verify the representation constant R(p)

```bash
python code/verify_Rp.py
```

Computes `R(p)` for primes up to a chosen bound and reports the convergence
of `R(p)·(ln p)²/p` toward ≈ 0.73.

### 3. Reproduce the exploratory scan

```bash
python code/odd_vs_prime_scan.py
```

Computes `R₋(n)` and `R₊(n)` for odd `n` up to `N_MAX`, saving
`odd_vs_prime_scan.csv`. Exploratory only — see [Status](#status-and-limitations).

### 4. Band-based analysis

```bash
python code/odd_vs_prime_band.py
```

Analyzes primes within square bands `[m², (m+1)²]` and computes symmetry and
ratio profiles.

### 5. Unified mod-4 comparison

```bash
python code/mod4_unified_analysis.py
```

Computes the unified score and dimension summaries across `t_bin`,
`dist_bin`, and `side_of_band`.

### 6. Singular series check

```bash
python code/singular_series.py
```

Compares the empirical 0.73 constant against the Hardy–Littlewood singular
series and the twin prime constant C₂.

---

## Data

- `R_p_results_10M.csv` — `R(p)` values for 664,578 primes up to 10⁷.
- `scale_summary.csv` — multi-scale comparison (25k to 250k) of exploratory
  ratio profiles.

---

## Status and Limitations

- The geometric sieve (Theorem 1) is proved. No false positives is a
  guarantee, not an empirical estimate.
- The 0.73 representation constant is a numerically supported conjecture. No
  proof of convergence is claimed.
- The exploratory extensions (mod-4 offset, band symmetry, unified score) are
  observational, based on finite-scale scans (N up to 250,000–10,000,000
  depending on the script), and have not been tested against a rigorous null
  model. Effect sizes there are small (unified score ≈ 0.048) and should be
  treated as a research direction, not a finding.

---

## License

MIT. See [`LICENSE`](./LICENSE) for details.

---

## Citation

If you use this work, please cite:

```
[To be updated after publication / preprint DOI assignment]
```
