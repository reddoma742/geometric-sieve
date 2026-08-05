"""
mod4_unified_analysis.py
-------------------------
Unified comparison of 4n+1 vs 4n-1 primes across three dimensions:
  - t_bin (normalised band position)
  - dist_bin (distance to nearest square)
  - side_of_band (left or right half of band)

Computes the Unified Score U: the mean absolute gap between the
two prime classes across all contexts.

Output:
    - mod4_unified_compare.csv              (full comparison data)
    - mod4_unified_dimension_summary.csv    (dimension summary)
    - figures/mod4_unified_analysis.png      (plot)
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ============================================================
# Parameters
# ============================================================
LIMIT = 200_000
N_MIN = 3
N_MAX = 50_000
N_T_BINS = 20

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
print(f"Found {len(primes)} primes up to {LIMIT}.")

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
# 3. Main scan (primes only)
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
    dist_bin = min(int(nearest_sq_dist), 30)

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
        'dist_bin': dist_bin,
        'side_of_band': side,
        'near_low_flag': near_low_flag,
    })

df = pd.DataFrame(rows)
print(f"Done. {len(df)} primes scanned.")
print(f"  prime_4k1: {len(df[df['prime_class'] == 'prime_4k1'])}")
print(f"  prime_4k3: {len(df[df['prime_class'] == 'prime_4k3'])}")

# ============================================================
# 4. Unified comparison across all contexts
# ============================================================
group_cols = ['t_bin', 'dist_bin', 'side_of_band', 'near_low_flag']

unified_rows = []
for keys, grp in df.groupby(group_cols):
    t, d, side, near_low = keys
    k1 = grp[grp['prime_class'] == 'prime_4k1']
    k3 = grp[grp['prime_class'] == 'prime_4k3']
    if len(k1) < 3 or len(k3) < 3:
        continue
    delta_sum = k1['sum_norm'].mean() - k3['sum_norm'].mean()
    delta_diff = k1['diff_norm'].mean() - k3['diff_norm'].mean()
    delta_ratio = k1['R_ratio'].dropna().mean() - k3['R_ratio'].dropna().mean()
    unified_rows.append({
        't_bin': t,
        'dist_bin': d,
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
# 5. Unified Score U
# ============================================================
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

dim_rows = []
for dim_name, dim_col in [('t_bin', 't_bin'),
                           ('dist_bin', 'dist_bin'),
                           ('side_of_band', 'side_of_band'),
                           ('near_low_flag', 'near_low_flag')]:
    grp = unified_df.groupby(dim_col).agg(
        delta_ratio=('delta_ratio', 'mean'),
        delta_sum_norm=('delta_sum_norm', 'mean'),
        count_4k1=('count_4k1', 'sum'),
        count_4k3=('count_4k3', 'sum')
    ).reset_index()
    grp['dimension'] = dim_name
    dim_rows.append(grp)
    print(f"\n{dim_name}:")
    for _, row in grp.iterrows():
        print(f"  {row[dim_col]}: delta_ratio={row['delta_ratio']:.6f}, "
              f"delta_sum_norm={row['delta_sum_norm']:.6f}")

dim_summary = pd.concat(dim_rows, ignore_index=True)
dim_summary.to_csv('mod4_unified_dimension_summary.csv', index=False)
print(f"\nDimension summary saved to mod4_unified_dimension_summary.csv")

# ============================================================
# 7. PLOTS — SAVE TO FILE
# ============================================================
os.makedirs("figures", exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: delta_ratio by t_bin
ax = axes[0, 0]
t_bin_grp = unified_df.groupby('t_bin')['delta_ratio'].mean()
t_bins_sorted = sorted(t_bin_grp.index)
ax.bar(t_bins_sorted, [t_bin_grp[t] for t in t_bins_sorted],
       color='steelblue', edgecolor='navy')
ax.axhline(y=unified_score, color='red', linestyle='--', linewidth=2,
           label=f'U = {unified_score:.4f}')
ax.set_xlabel('t_bin')
ax.set_ylabel('Mean delta_ratio')
ax.set_title('delta_ratio by t_bin')
ax.legend()
ax.grid(alpha=0.3, axis='y')

# Plot 2: delta_ratio by dist_bin (top 15)
ax = axes[0, 1]
dist_grp = unified_df.groupby('dist_bin')['delta_ratio'].agg(['mean', 'count'])
dist_grp = dist_grp[dist_grp['count'] >= 3].head(15)
ax.bar(dist_grp.index.astype(str), dist_grp['mean'],
       color='darkgreen', edgecolor='darkgreen')
ax.axhline(y=unified_score, color='red', linestyle='--', linewidth=2,
           label=f'U = {unified_score:.4f}')
ax.set_xlabel('dist_bin')
ax.set_ylabel('Mean delta_ratio')
ax.set_title('delta_ratio by dist_bin')
ax.legend()
ax.grid(alpha=0.3, axis='y')
ax.tick_params(axis='x', rotation=45)

# Plot 3: delta_ratio by side_of_band
ax = axes[1, 0]
side_grp = unified_df.groupby('side_of_band').agg(
    mean_delta=('delta_ratio', 'mean'),
    count=('delta_ratio', 'count')
)
colors_side = ['blue', 'darkorange']
for i, (side, row) in enumerate(side_grp.iterrows()):
    ax.bar(side, row['mean_delta'], color=colors_side[i % len(colors_side)],
           edgecolor='black', alpha=0.8)
ax.axhline(y=unified_score, color='red', linestyle='--', linewidth=2,
           label=f'U = {unified_score:.4f}')
ax.set_ylabel('Mean delta_ratio')
ax.set_title('delta_ratio by side_of_band')
ax.legend()
ax.grid(alpha=0.3, axis='y')

# Plot 4: scatter of delta_ratio vs delta_sum_norm
ax = axes[1, 1]
scatter = ax.scatter(unified_df['delta_sum_norm'], unified_df['delta_ratio'],
                     c=unified_df['t_bin'], cmap='viridis', alpha=0.7, edgecolors='black', linewidth=0.3)
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('delta_sum_norm')
ax.set_ylabel('delta_ratio')
ax.set_title('delta_ratio vs delta_sum_norm (colored by t_bin)')
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('t_bin')
ax.grid(alpha=0.3)

plt.suptitle('Mod-4 Unified Analysis — Geometric Sieve', fontsize=14, fontweight='bold')
plt.tight_layout()

output_path = os.path.join("figures", "mod4_unified_analysis.png")
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"\nFigure saved to: {output_path}")

# ============================================================
# 8. Save data
# ============================================================
unified_df.to_csv('mod4_unified_compare.csv', index=False)
print(f"Data saved to mod4_unified_compare.csv")
print(f"\nUnified Score U = {unified_score:.15f}")
print(f"\nDone.")
