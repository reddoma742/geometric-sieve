"""
odd_vs_prime_band.py
--------------------
Analyzes odd numbers within square bands [m^2, (m+1)^2].

For each odd number, computes:
  - band_t: normalised position inside the band (0 = left, 1 = right)
  - nearest_square_dist: distance to nearest square boundary
  - side_of_band: left or right half

Compares prime vs. composite behaviour within bands, and 4n+1 vs 4n-1.

Output:
    - odd_vs_prime_band_scan.csv       (full data)
    - odd_vs_prime_band_summary.csv    (group summary)
    - odd_vs_prime_band_symmetry.csv   (symmetry analysis)
    - figures/odd_vs_prime_band_scan.png (plot)
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

    # Classify prime type
    if is_prime[n]:
        if n % 4 == 1:
            prime_class = 'prime_4k1'
        else:
            prime_class = 'prime_4k3'
    else:
        prime_class = 'composite'

    # dist_bin (distance from nearest square, binned)
    dist_bin = min(int(nearest_sq_dist), 30)

    rows.append({
        'n': n,
        'is_prime_n': bool(is_prime[n]),
        'prime_class': prime_class,
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
        'dist_bin': dist_bin,
        'side_of_band': side,
        'm': m,
        'lower_sq': lower_sq,
        'upper_sq': upper_sq,
    })

    if n % 5001 == 3:
        print(f"  n = {n} ...")

df = pd.DataFrame(rows)
print(f"Done. {len(df)} odd numbers scanned.")

# ============================================================
# 4. Group summaries
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY BY GROUP")
print(f"{'='*70}")

group_defs = {
    'all_odd': df,
    'prime_4k1': df[df['prime_class'] == 'prime_4k1'],
    'prime_4k3': df[df['prime_class'] == 'prime_4k3'],
    'composite': df[df['prime_class'] == 'composite'],
}

# Add band_middle and band_edge
mid_mask = (df['band_t'] >= 0.25) & (df['band_t'] <= 0.75)
edge_mask = ~mid_mask
group_defs['band_middle'] = df[mid_mask]
group_defs['band_edge'] = df[edge_mask]

# Add near_square and far_square
near_sq = df['nearest_square_dist'] < 5
far_sq = df['nearest_square_dist'] >= 5
group_defs['near_square'] = df[near_sq]
group_defs['far_square'] = df[far_sq]

summary_rows = []
for name, g in group_defs.items():
    valid_ratio = g['R_ratio'].replace([np.inf, -np.inf], np.nan).dropna()
    s = {
        'group': name,
        'count': len(g),
        'mean_minus_norm': g['minus_norm'].mean(),
        'mean_plus_norm': g['plus_norm'].mean(),
        'mean_sum_norm': g['sum_norm'].mean(),
        'mean_diff_norm': g['diff_norm'].mean(),
        'mean_ratio': valid_ratio.mean() if len(valid_ratio) else np.nan,
        'mean_band_t': g['band_t'].mean(),
        'mean_center_dist': g['center_dist'].mean(),
        'mean_nearest_square_dist': g['nearest_square_dist'].mean(),
    }
    summary_rows.append(s)
    print(f"\n{name} (n={len(g)}):")
    print(f"  mean_sum_norm  = {g['sum_norm'].mean():.4f}")
    print(f"  mean_ratio     = {valid_ratio.mean():.4f}" if len(valid_ratio) else "  mean_ratio     = N/A")
    print(f"  mean_band_t    = {g['band_t'].mean():.4f}")

summary = pd.DataFrame(summary_rows)
summary.to_csv('odd_vs_prime_band_summary.csv', index=False)
print(f"\nSummary saved to odd_vs_prime_band_summary.csv")

# ============================================================
# 5. Symmetry analysis (left vs right within bands)
# ============================================================
print(f"\n{'='*70}")
print("SYMMETRY (left vs right within bands)")
print(f"{'='*70}")

sym_rows = []
for name in ['all_odd', 'prime_4k1', 'prime_4k3', 'composite']:
    if name not in group_defs:
        continue
    g = group_defs[name]
    left = g[g['side_of_band'] == 'left_of_band']
    right = g[g['side_of_band'] == 'right_of_band']
    if len(left) > 0 and len(right) > 0:
        sym_diff = abs(left['R_sum'].mean() - right['R_sum'].mean())
        signed_gap = right['R_sum'].mean() - left['R_sum'].mean()
        sym_rows.append({
            'group': name,
            'sym_diff': sym_diff,
            'signed_gap': signed_gap,
        })
        print(f"  {name}: sym_diff={sym_diff:.6f}, signed_gap={signed_gap:.6f}")

sym_df = pd.DataFrame(sym_rows)
sym_df.to_csv('odd_vs_prime_band_symmetry.csv', index=False)

# ============================================================
# 6. PLOTS — SAVE TO FILE
# ============================================================
os.makedirs("figures", exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Plot 1: sum_norm vs n, colored by class
ax = axes[0, 0]
colors = {'prime_4k1': 'blue', 'prime_4k3': 'green', 'composite': 'red'}
alphas = {'prime_4k1': 0.5, 'prime_4k3': 0.5, 'composite': 0.1}
sizes = {'prime_4k1': 4, 'prime_4k3': 4, 'composite': 1}
for cls in ['prime_4k1', 'prime_4k3', 'composite']:
    sub = df[df['prime_class'] == cls]
    ax.plot(sub['n'], sub['sum_norm'], '.', alpha=alphas[cls],
            markersize=sizes[cls], color=colors[cls], label=cls)
ax.set_xscale('log')
ax.set_xlabel('n')
ax.set_ylabel('sum_norm')
ax.set_title('sum_norm vs n')
ax.legend(markerscale=3)
ax.grid(alpha=0.3)

# Plot 2: Binned means
ax = axes[0, 1]
bins = np.logspace(np.log10(N_MIN), np.log10(N_MAX), 20)
for cls, col in [('prime_4k1', 'blue'), ('prime_4k3', 'green'), ('composite', 'red')]:
    sub = df[df['prime_class'] == cls]
    ctrs, means = [], []
    for i in range(len(bins)-1):
        msk = (sub['n'] >= bins[i]) & (sub['n'] < bins[i+1])
        if msk.sum() > 5:
            ctrs.append(np.sqrt(bins[i]*bins[i+1]))
            means.append(sub.loc[msk, 'sum_norm'].mean())
    ax.plot(ctrs, means, 'o-', color=col, label=cls, markersize=4)
ax.set_xscale('log')
ax.set_xlabel('n')
ax.set_ylabel('Binned mean sum_norm')
ax.set_title('Binned sum_norm')
ax.legend()
ax.grid(alpha=0.3)

# Plot 3: sum_norm vs band_t (scatter)
ax = axes[0, 2]
for cls, col in [('prime_4k1', 'blue'), ('prime_4k3', 'green')]:
    sub = df[df['prime_class'] == cls]
    ax.plot(sub['band_t'], sub['sum_norm'], '.', alpha=0.4,
            markersize=3, color=col, label=cls)
ax.set_xlabel('band_t')
ax.set_ylabel('sum_norm')
ax.set_title('sum_norm vs band_t (primes only)')
ax.legend()
ax.grid(alpha=0.3)

# Plot 4: ratio vs band_t (scatter)
ax = axes[1, 0]
for cls, col in [('prime_4k1', 'blue'), ('prime_4k3', 'green')]:
    sub = df[df['prime_class'] == cls]
    valid = sub['R_ratio'].replace([np.inf, -np.inf], np.nan).dropna()
    valid_band_t = sub.loc[valid.index, 'band_t']
    ax.plot(valid_band_t, valid, '.', alpha=0.4, markersize=3,
            color=col, label=cls)
ax.set_xlabel('band_t')
ax.set_ylabel('R_plus / R_minus')
ax.set_title('Ratio vs band_t (primes only)')
ax.legend()
ax.grid(alpha=0.3)

# Plot 5: Histogram — sum_norm for primes
ax = axes[1, 1]
for cls, col in [('prime_4k1', 'blue'), ('prime_4k3', 'green')]:
    sub = df[df['prime_class'] == cls]
    ax.hist(sub['sum_norm'], bins=30, alpha=0.5, color=col, label=cls)
ax.set_xlabel('sum_norm')
ax.set_ylabel('Frequency')
ax.set_title('Histogram: sum_norm (4n+1 vs 4n-1)')
ax.legend()
ax.grid(alpha=0.3)

# Plot 6: Histogram — R_ratio for primes
ax = axes[1, 2]
for cls, col in [('prime_4k1', 'blue'), ('prime_4k3', 'green')]:
    sub = df[df['prime_class'] == cls]
    valid = sub['R_ratio'].replace([np.inf, -np.inf], np.nan).dropna()
    ax.hist(valid, bins=30, alpha=0.5, color=col, label=cls)
ax.set_xlabel('R_plus / R_minus')
ax.set_ylabel('Frequency')
ax.set_title('Histogram: Ratio (4n+1 vs 4n-1)')
ax.legend()
ax.grid(alpha=0.3)

plt.suptitle('Odd vs Prime Band Analysis — Geometric Sieve', fontsize=14, fontweight='bold')
plt.tight_layout()

output_path = os.path.join("figures", "odd_vs_prime_band_scan.png")
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"\nFigure saved to: {output_path}")

# ============================================================
# 7. Save data
# ============================================================
df.to_csv('odd_vs_prime_band_scan.csv', index=False)
print(f"Data saved to odd_vs_prime_band_scan.csv")
print(f"\nDone.")
