
markdown
# The Geometric Sieve

A structural framework for analyzing prime number distribution based on square layers.

## Author

**Reddouane BERRAMDANE**
- Email: reddoma@gmail.com
- Date: August 2026

## Overview

This repository contains the complete implementation and data for the Geometric Sieve framework.

### Key Results

| Result | Status | Description |
|:---|:---|:---|
| **Theorem 1** | ✅ Proven | gcd(\|b\|, C(a)) > 1 ⇒ composite (zero false positives) |
| **Conjecture 1** | ⚠️ Empirical | R(p) · (ln p)² / p → 0.73 (664,578 primes tested) |
| **Golden K Arcs** | 📊 Observed | K ∈ {2, 3, 5, 11, 17, 41} produce longest prime chains |

## Quick Start

```bash
pip install -r requirements.txt
python code/verify_sieve.py
python code/verify_Rp.py
python code/geometric_sieve_full.py
Repository Structure
text
geometric-sieve/
├── README.md
├── LICENSE
├── requirements.txt
├── code/
│   ├── verify_sieve.py
│   ├── verify_Rp.py
│   ├── geometric_sieve_full.py
│   ├── generate_Rp_data.py
│   └── scan_10M_results.py
├── data/
│   ├── R_p_results_10M.csv
│   ├── scan_10M_results.csv
│   └── heegner_scan_top30.csv
└── Rp_cpp/
    └── main.cpp
Citation
text
BERRAMDANE, Reddouane. The Geometric Sieve: A Structural Framework 
for Prime Distribution. August 2026. GitHub: [repository-url]
License
MIT License - see LICENSE file

text

---

## 📄 ملف 4: `paper_geometric_sieve.pdf`

هذا الملف موجود لديك من Claude. تأكد من تحديث:

- **Author**: Reddouane BERRAMDANE
- **Email**: reddoma@gmail.com
- **GitHub URL**: [(https://github.com/reddoma742/geometric-sieve.git)
