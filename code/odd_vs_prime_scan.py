"""
odd_vs_prime_scan.py
--------------------
Compares odd primes vs. odd composites using the representation
functions R_minus(n) and R_plus(n), where:

    R_minus(n) = #{a in P, a < n : 2a - n in {1} U P}
    R_plus(n)  = #{a in P, 2a < n : n - 2a in {1} U P}

Normalised sums are compared across groups to reveal structural
differences between primes and composites.
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    """Count representations n = 2a - b (a prime, b prime or 1)"""
    cnt = 0
    for a in primes:
        if a >= n:
            break
        b = 2 * a - n
        if b >= 1 and is_prime_or_one(int(b)):
            cnt += 1
    return cnt

def count_plus_any(n: int) -> int:
    """Count representations n = 2a + b (a prime, b prime or 1)"""
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

print(f"  Primes   : {len(prime_df)}")
print(f"  Composites: {len(composite_df)}")

# ============================================================
# 5. Summary statistics
# ============================================================
print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")

for label, sub in [('All odd', df), ('Odd primes', prime_df), ('Odd composites', composite_df)]:
    print(f"\n{label} (n={len(sub)}):")
    print(f"  mean sum_norm  = {sub['sum_norm'].mean():.4f}")
    print(f"  mean minus_norm = {sub['minus_norm'].mean():.4f}")
    print(f"  mean plus_norm  = {sub['plus_norm'].mean():.4f}")
    print(f"  mean R_ratio    = {sub['R_ratio'].dropna().mean():.4f}")

# ============================================================
# 6. Save data
# ============================================================
output_file = "odd_vs_prime_scan.csv"
df.to_csv(output_file, index=False)
print(f"\nSaved to {output_file}")
