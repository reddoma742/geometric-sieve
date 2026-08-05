"""
mod4_unified_analysis.py
-------------------------
Unified comparison of 4n+1 vs 4n-1 primes across three dimensions:
  - t_bin (normalised band position)
  - dist_bin (distance to nearest square)
  - side_of_band (left or right half of band)

Computes the Unified Score U: the mean absolute gap between the
two prime classes across all contexts.
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Parameters
# ============================================================
LIMIT = 200_000
N_MIN = 3
N_MAX = 50_000
N_T_BINS = 20
N_DIST_BINS = 20

# ============================================================
# 1. Sieve
# ============================================================
def sieve(limit: int):
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    primes = np.nonzero(is_prime)[0]
    return is_prime, primes

print("Sieving...")
is_prime, primes = sieve(LIMIT)

def is_prime_or_one(n: int) -> bool:
    return n == 1 or (0 <= n < len(is_prime) and is_prime[n])

# ============================================================
# 2. Representation counts
# ============================================================
def count_minus_any(n: int) -> int:
    cnt = 0
    for a in primes:
        if a >= n: break
        b = 2*a - n
        if b >= 1 and is_prime_or_one(int(b)): cnt += 1
    return cnt

def count_plus_any(n: int) -> int:
    cnt = 0
    for a in primes:
        if 2*a >= n: break
        b = n - 2*a
        if b >= 1 and is_prime_or_one(int(b)): cnt += 1
    return cnt

# ============================================================
# 3. Main scan
# ============================================================
print(f"Scanning primes from {N_MIN} to {N_MAX}...")
rows = []

for n in range(N_MIN, N_MAX + 1, 2):
    if not is_prime[n]:
        continue

    r_minus = count_minus_any(n)
    r_plus = count_plus_any(n)
    scale = n / (math.log(max(n, 3))**2)

    m = int(math.isqrt(n))
    lower_sq = m * m
    upper_sq = (m + 1) * (m + 1)
    band_width = upper_sq - lower_sq
    band_t = (n - lower_sq) / band_width if band_width > 0 else 0.0
    t_bin = int(band_t * N_T_BINS)
    if t_bin >= N_T_BINS:
        t_bin = N_T_BINS - 1

    dist_lower = n - lower_sq
    dist_upper = upper_sq - n
    nearest_sq_dist = min(dist_lower, dist_upper)
    side = 'left_of_band' if band_t <= 0.5 else 'right_of_band'
    near_low_flag = dist_lower <= dist_upper

    rows.append({
        'n': n,
        'prime_class': 'prime_4k1' if n % 4 == 1 else 'prime_4k3',
        'R_minus': r_minus,
        'R_plus': r_plus,
        'R_sum': r_minus + r_plus,
        'R_diff': r_plus - r_minus,
        'R_ratio': (r_plus / r_minus) if r_minus > 0 else np.nan,
        'sum_norm': (r_minus + r_plus) / scale,
        'diff_norm': (r_plus - r_minus) / scale,
        'band_t': band_t,
        't_bin': t_bin,
        'nearest_square_dist': nearest_sq_dist,
        'side_of_band': side,
        'near_low_flag': near_low_flag,
    })

df = pd.DataFrame(rows)
print(f"Done. {len(df)} primes.")

# ============================================================
# 4. Unified comparison
# ============================================================
# Compare prime_4k1 vs prime_4k3 across t_bin, side_of_band, near_low_flag
group_cols = ['t_bin', 'side_of_band', 'near_low_flag']

unified_rows = []
for (t, side, near_low), grp in df.groupby(group_cols):
    k1 = grp[grp['prime_class'] == 'prime_4k1']
    k3 = grp[grp['prime_class'] == 'prime_4k3']
    if len(k1) < 5 or len(k3) < 5:
        continue
    delta_sum = k1['sum_norm'].mean() - k3['sum_norm'].mean()
    delta_diff = k1['diff_norm'].mean() - k3['diff_norm'].mean()
    delta_ratio = k1['R_ratio'].dropna().mean() - k3['R_ratio'].dropna().mean()
    unified_rows.append({
        't_bin': t,
        'side_of_band': side,
        'near_low_flag': near_low,
        'count_4k1': len(k1),
        'count_4k3': len(k3),
        'delta_sum_norm': delta_sum,
        'delta_diff_norm': delta_diff,
        'delta_ratio': delta_ratio,
    })

unified_df = pd.DataFrame(unified_rows)

# ============================================================
# 5. Unified Score
# ============================================================
# Mean absolute delta_ratio across all contexts
unified_score = unified_df['delta_ratio'].abs().mean()
print(f"\n{'='*60}")
print(f"UNIFIED SCORE")
print(f"{'='*60}")
print(f"U = {unified_score:.15f}")

# ============================================================
# 6. Dimension summaries
# ============================================================
print(f"\n{'='*60}")
print(f"DIMENSION SUMMARY")
print(f"{'='*60}")

for dim_name, dim_col in [('t_bin', 't_bin'),
                           ('side_of_band', 'side_of_band'),
                           ('near_low_flag', 'near_low_flag')]:
    grp = unified_df.groupby(dim_col).agg(
        delta_ratio=('delta_ratio', 'mean'),
        delta_sum_norm=('delta_sum_norm', 'mean'),
        count_4k1=('count_4k1', 'sum'),
        count_4k3=('count_4k3', 'sum')
    ).reset_index()
    grp['dimension'] = dim_name
    print(f"\n{dim_name}:")
    for _, row in grp.iterrows():
        print(f"  {row[dim_col]}: delta_ratio={row['delta_ratio']:.6f}, "
              f"delta_sum_norm={row['delta_sum_norm']:.6f}")

# ============================================================
# 7. Save
# ============================================================
unified_df.to_csv('mod4_unified_compare.csv', index=False)
print(f"\nSaved to mod4_unified_compare.csv")
print(f"Unified Score U = {unified_score:.15f}")
