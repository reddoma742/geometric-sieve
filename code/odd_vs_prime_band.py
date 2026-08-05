"""
odd_vs_prime_band.py
--------------------
Analyzes odd numbers within square bands [m^2, (m+1)^2].

For each odd number, computes:
  - band_t: normalised position inside the band (0 = left, 1 = right)
  - nearest_square_dist: distance to nearest square boundary
  - side: left_of_band or right_of_band

Then compares prime vs. composite behaviour within bands.
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
N_BINS = 20

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
# 3. Band computation
# ============================================================
def get_band_info(n):
    """Return (m, band_t, lower_sq, upper_sq, nearest_sq_dist, side)"""
    m = int(math.isqrt(n))
    lower_sq = m * m
    upper_sq = (m + 1) * (m + 1)
    band_width = upper_sq - lower_sq
    band_t = (n - lower_sq) / band_width if band_width > 0 else 0.0
    dist_lower = n - lower_sq
    dist_upper = upper_sq - n
    nearest_sq_dist = min(dist_lower, dist_upper)
    side = 'left_of_band' if band_t <= 0.5 else 'right_of_band'
    return m, band_t, lower_sq, upper_sq, nearest_sq_dist, side

print(f"Scanning odd n from {N_MIN} to {N_MAX}...")
rows = []

for n in range(N_MIN, N_MAX + 1, 2):
    r_minus = count_minus_any(n)
    r_plus = count_plus_any(n)
    scale = n / (math.log(max(n, 3))**2)

    m, band_t, lower_sq, upper_sq, nearest_sq_dist, side = get_band_info(n)
    center_dist = abs(band_t - 0.5)

    rows.append({
        'n': n,
        'is_prime_n': bool(is_prime[n]),
        'prime_type': 'prime_4k1' if (is_prime[n] and n % 4 == 1) else
                      ('prime_4k3' if is_prime[n] else 'composite'),
        'R_minus': r_minus,
        'R_plus': r_plus,
        'R_sum': r_minus + r_plus,
        'R_diff': r_plus - r_minus,
        'R_ratio': (r_plus / r_minus) if r_minus > 0 else np.nan,
        'minus_norm': r_minus / scale,
        'plus_norm': r_plus / scale,
        'sum_norm': (r_minus + r_plus) / scale,
        'diff_norm': (r_plus - r_minus) / scale,
        'band_t': band_t,
        'center_dist': center_dist,
        'nearest_square_dist': nearest_sq_dist,
        'side_of_band': side,
        'm': m,
        'lower_sq': lower_sq,
        'upper_sq': upper_sq,
    })

df = pd.DataFrame(rows)
print(f"Done. {len(df)} odd numbers.")

# ============================================================
# 4. Group summaries
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY BY GROUP")
print(f"{'='*70}")

groups = {
    'all_odd': df,
    'odd_prime': df[df['is_prime_n']],
    'odd_composite': df[(~df['is_prime_n']) & (df['n'] % 2 == 1)],
}

# Add band_middle and band_edge
mid_mask = (df['band_t'] >= 0.25) & (df['band_t'] <= 0.75)
edge_mask = ~mid_mask
groups['band_middle'] = df[mid_mask]
groups['band_edge'] = df[edge_mask]

# Add near_square and far_square
near_sq = df['nearest_square_dist'] < 5
far_sq = df['nearest_square_dist'] >= 5
groups['near_square'] = df[near_sq]
groups['far_square'] = df[far_sq]

for name, g in groups.items():
    print(f"\n{name} (n={len(g)}):")
    print(f"  mean_sum_norm  = {g['sum_norm'].mean():.4f}")
    print(f"  mean_ratio     = {g['R_ratio'].dropna().mean():.4f}")
    print(f"  mean_band_t    = {g['band_t'].mean():.4f}")

# ============================================================
# 5. Symmetry analysis
# ============================================================
print(f"\n{'='*70}")
print("SYMMETRY (left vs right within bands)")
print(f"{'='*70}")

for name in ['all_odd', 'odd_prime', 'odd_composite']:
    if name not in groups:
        continue
    g = groups[name]
    left = g[g['side_of_band'] == 'left_of_band']
    right = g[g['side_of_band'] == 'right_of_band']
    if len(left) > 0 and len(right) > 0:
        sym_diff = abs(left['R_sum'].mean() - right['R_sum'].mean())
        signed_gap = right['R_sum'].mean() - left['R_sum'].mean()
        print(f"  {name}: sym_diff={sym_diff:.6f}, signed_gap={signed_gap:.6f}")

# ============================================================
# 6. Save
# ============================================================
output_file = "odd_vs_prime_band_scan.csv"
df.to_csv(output_file, index=False)
print(f"\nSaved to {output_file}")
