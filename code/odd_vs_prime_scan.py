"""
odd_vs_prime_scan.py
--------------------
Compares odd primes vs. odd composites using the representation
functions R_minus(n) and R_plus(n), where:

    R_minus(n) = #{a in P, a < n : 2a - n in {1} U P}
    R_plus(n)  = #{a in P, 2a < n : n - 2a in {1} U P}

Normalised sums are compared across groups to reveal structural
differences between primes and composites.

Output:
    - odd_vs_prime_scan.csv       (full data)
    - odd_vs_prime_summary.csv    (group summary)
    - figures/odd_vs_prime_scan.png (plot)
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ============================================================
# Parameters
# ============================================================
LIMIT = 200_000          # sieve limit (must be >= 2 * N_MAX)
N_MIN = 3
N_MAX = 50_000

# ============================================================
# 1. Sieve of Eratosthenes
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

# ============================================================
# 2. Helper functions
# ============================================================
def is_prime_or_one(n: int) -> bool:
    return n == 1 or (0 <= n < len(is_prime) and is_prime[n])

def count_minus_any(n: int) -> int:
    cnt = 0
    for a in primes:
        if a >= n:
            break
        b = 2 * a - n
        if b >= 1 and is_prime_or_one(int(b)):
            cnt += 1
    return cnt

def count_plus_any(n: int) -> int:
    cnt = 0
    for a in primes:
        if 2 * a >= n:
            break
        b = n - 2 * a
        if b >= 1 and is_prime_or_one(int(b)):
            cnt += 1
    return cnt

# ============================================================
# 3. Main computation
# ============================================================
print(f"Scanning odd n from {N_MIN} to {N_MAX}...")
rows = []

for n in range(N_MIN, N_MAX + 1, 2):
    r_minus = count_minus_any(n)
    r_plus = count_plus_any(n)
    scale = n / (math.log(max(n, 3))**2)

    rows.append({
        'n': n,
        'is_prime_n': bool(is_prime[n]),
        'R_minus': r_minus,
        'R_plus': r_plus,
        'R_sum': r_minus + r_plus,
        'R_diff': r_plus - r_minus,
        'R_ratio': (r_plus / r_minus) if r_minus > 0 else np.nan,
        'minus_norm': r_minus / scale,
        'plus_norm': r_plus / scale,
        'sum_norm': (r_minus + r_plus) / scale,
    })

    if n % 5001 == 3:
        print(f"  n = {n} ...")

df = pd.DataFrame(rows)
print(f"Done. {len(df)} odd numbers scanned.")

# ============================================================
# 4. Split into groups
# ============================================================
prime_df = df[df['is_prime_n']].copy()
composite_df = df[(~df['is_prime_n']) & (df['n'] % 2 == 1)].copy()

print(f"  Primes    : {len(prime_df)}")
print(f"  Composites: {len(composite_df)}")

# ============================================================
# 5. Summary statistics
# ============================================================
print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")

summary_rows = []
for label, sub in [('all_odd', df), ('odd_prime', prime_df), ('odd_composite', composite_df)]:
    valid_ratio = sub['R_ratio'].replace([np.inf, -np.inf], np.nan).dropna()
    s = {
        'group': label,
        'count': len(sub),
        'mean_minus_norm': sub['minus_norm'].mean(),
        'mean_plus_norm': sub['plus_norm'].mean(),
        'mean_sum_norm': sub['sum_norm'].mean(),
        'mean_diff_norm': (sub['R_plus'] - sub['R_minus']).mean() / (sub['n'] / (np.log(sub['n'])**2)).mean() if len(sub) > 0 else np.nan,
        'mean_ratio': valid_ratio.mean() if len(valid_ratio) else np.nan,
    }
    summary_rows.append(s)
    print(f"\n{label} (n={len(sub)}):")
    print(f"  mean sum_norm   = {sub['sum_norm'].mean():.4f}")
    print(f"  mean minus_norm = {sub['minus_norm'].mean():.4f}")
    print(f"  mean plus_norm  = {sub['plus_norm'].mean():.4f}")
    print(f"  mean R_ratio    = {valid_ratio.mean():.4f}" if len(valid_ratio) else "  mean R_ratio    = N/A")

summary = pd.DataFrame(summary_rows)
summary.to_csv('odd_vs_prime_summary.csv', index=False)
print(f"\nSummary saved to odd_vs_prime_summary.csv")

# ============================================================
# 6. Binned means for plotting
# ============================================================
bins = np.logspace(np.log10(df['n'].min()), np.log10(df['n'].max()), 25)
centers = []
all_means = []
prime_means = []
composite_means = []

for i in range(len(bins) - 1):
    mask_all = (df['n'] >= bins[i]) & (df['n'] < bins[i + 1])
    mask_prime = (prime_df['n'] >= bins[i]) & (prime_df['n'] < bins[i + 1])
    mask_comp = (composite_df['n'] >= bins[i]) & (composite_df['n'] < bins[i + 1])
    
    if mask_all.sum() > 10:
        centers.append((bins[i] * bins[i + 1]) ** 0.5)
        all_means.append(df.loc[mask_all, 'sum_norm'].mean())
        prime_means.append(prime_df.loc[mask_prime, 'sum_norm'].mean() if mask_prime.sum() > 0 else np.nan)
        composite_means.append(composite_df.loc[mask_comp, 'sum_norm'].mean() if mask_comp.sum() > 0 else np.nan)

centers = np.array(centers)
all_means = np.array(all_means)
prime_means = np.array(prime_means)
composite_means = np.array(composite_means)

# ============================================================
# 7. PLOTS — SAVE TO FILE
# ============================================================
os.makedirs("figures", exist_ok=True)

plt.figure(figsize=(14, 10))

# Plot 1: Scatter plot
plt.subplot(2, 2, 1)
plt.plot(df['n'], df['sum_norm'], '.', alpha=0.08, markersize=2,
         color='gray', label='All odd n')
plt.plot(prime_df['n'], prime_df['sum_norm'], '.', alpha=0.15, markersize=2,
         color='blue', label='Odd primes')
plt.plot(composite_df['n'], composite_df['sum_norm'], '.', alpha=0.12, markersize=2,
         color='red', label='Odd composites')
plt.xscale('log')
plt.xlabel('n')
plt.ylabel('sum_norm')
plt.title('R_sum / (n / (log n)^2)')
plt.legend()
plt.grid(alpha=0.3)

# Plot 2: Binned means
plt.subplot(2, 2, 2)
plt.plot(centers, all_means, color='black', linewidth=2, label='All odd')
plt.plot(centers, prime_means, color='blue', linewidth=2, label='Odd primes')
plt.plot(centers, composite_means, color='red', linewidth=2, label='Odd composites')
plt.xscale('log')
plt.xlabel('n')
plt.ylabel('Binned mean sum_norm')
plt.title('Binned comparison')
plt.legend()
plt.grid(alpha=0.3)

# Plot 3: Histogram — sum_norm
plt.subplot(2, 2, 3)
plt.hist(prime_df['sum_norm'], bins=40, alpha=0.7, color='blue', label='Odd primes')
plt.hist(composite_df['sum_norm'], bins=40, alpha=0.5, color='red', label='Odd composites')
plt.xlabel('sum_norm')
plt.ylabel('Frequency')
plt.title('Histogram comparison: sum_norm')
plt.legend()
plt.grid(alpha=0.3)

# Plot 4: Histogram — R_ratio
plt.subplot(2, 2, 4)
valid_prime_ratio = prime_df['R_ratio'].replace([np.inf, -np.inf], np.nan).dropna()
valid_comp_ratio = composite_df['R_ratio'].replace([np.inf, -np.inf], np.nan).dropna()
plt.hist(valid_prime_ratio, bins=40, alpha=0.7, color='blue', label='Odd primes')
plt.hist(valid_comp_ratio, bins=40, alpha=0.5, color='red', label='Odd composites')
plt.xlabel('R_plus / R_minus')
plt.ylabel('Frequency')
plt.title('Ratio comparison')
plt.legend()
plt.grid(alpha=0.3)

plt.suptitle('Odd vs Prime Scan — Geometric Sieve', fontsize=14, fontweight='bold')
plt.tight_layout()

# SAVE
output_path = os.path.join("figures", "odd_vs_prime_scan.png")
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"\nFigure saved to: {output_path}")

# ============================================================
# 8. Save data
# ============================================================
df.to_csv('odd_vs_prime_scan.csv', index=False)
print(f"Data saved to odd_vs_prime_scan.csv")
print(f"\nDone.")
