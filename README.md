markdown
# Geometric Sieve

A numerical study of prime distributions through the geometric framework
**n = 2a ± b**, where *a* and *b* are primes or 1.

---

## Overview

This repository explores the representation of odd integers as:
n = 2a + b or n = 2a - b

text

where **a** is prime and **b** is either 1 or prime.

The study documents observed differences in representation density, symmetry,
and ratio between:

- **Prime vs. composite** odd integers
- **Primes of the form 4n+1 vs. 4n-1**

across scales from 25,000 to 10,000,000.

---

## Key Observations

| Observation | Description |
|-------------|-------------|
| **Isolation** | Primes have fewer representations than composites (~1.69 vs. ~2.79 mean `sum_norm`) |
| **Symmetry** | Primes show near-perfect left-right symmetry within square bands |
| **Topological charge** | 4n+1 and 4n-1 primes exhibit opposite attraction to square boundaries |
| **Unified Score** | The gap between 4n+1 and 4n-1 converges to **U ≈ 0.04825** across all dimensions |
| **Standing wave** | The gap oscillates with a period of ~4 band units |

---

## Repository Structure
├── README.md
├── requirements.txt
├── LICENSE
├── code/
│ ├── verify_sieve.py
│ ├── verify_Rp.py
│ ├── singular_series.py
│ ├── odd_vs_prime_scan.py
│ ├── odd_vs_prime_band.py
│ └── mod4_unified_analysis.py
├── data/
│ ├── R_p_results_10M.csv
│ └── scale_summary.csv
└── figures/
└── (generated plots)

text

---

## Requirements

Install the required packages:

```bash
pip install -r requirements.txt
Usage
1. Reproduce the main scan
bash
python code/odd_vs_prime_scan.py
Computes R₋(n) and R₊(n) for all odd n up to N_MAX and saves
odd_vs_prime_scan.csv.

2. Band-based analysis
bash
python code/odd_vs_prime_band.py
Analyzes primes within square bands [m², (m+1)²] and computes
symmetry and ratio profiles.

3. Unified mod-4 comparison
bash
python code/mod4_unified_analysis.py
Computes the Unified Score U and dimension summaries across
t_bin, dist_bin, and side_of_band.

4. Singular series verification
bash
python code/singular_series.py
Compares the empirical constant with the Hardy–Littlewood
singular series and the twin prime constant C₂.

Data
R_p_results_10M.csv: R(p) values for 664,578 primes up to 10⁷

scale_summary.csv: Multi-scale comparison (25K to 250K) of ratio profiles

License
MIT. See LICENSE for details.

Citation
If you use this work, please cite:

text
[To be updated after publication]
text

---

## ✅ الخطوة 2: requirements.txt
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0

text

---

## ✅ الخطوة 3: LICENSE

رخصة MIT القياسية:
MIT License

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
