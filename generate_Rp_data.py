"""
generate_Rp_data.py
Generate R(p) data for primes up to a specified limit.
Outputs CSV file with full statistics.

Usage: python generate_Rp_data.py [limit]
Default limit: 1,000,000
"""

import math
import numpy as np
import csv
import sys
import time


def sieve_of_eratosthenes(limit):
    """Generate primes up to limit."""
    is_p = np.ones(limit + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if is_p[i]:
            is_p[i * i::i] = False
    return is_p


def compute_R(p, primes_list, is_prime_arr):
    """Compute R(p)."""
    cnt = 0
    for a in primes_list:
        if a >= p:
            break
        b = 2 * a - p
        if b >= 1 and (b == 1 or (b < len(is_prime_arr) and is_prime_arr[b])):
            cnt += 1
    return cnt


def generate_data(limit=1_000_000):
    """Generate full R(p) dataset."""
    print(f"⏳ Generating primes up to {limit:,}...")
    t0 = time.time()
    
    is_prime = sieve_of_eratosthenes(limit)
    primes = np.nonzero(is_prime)[0].tolist()
    
    t1 = time.time()
    print(f"   Generated {len(primes):,} primes in {t1-t0:.2f}s")
    print(f"⏳ Computing R(p) for each prime...")
    
    results = []
    
    for idx, p in enumerate(primes):
        if p <= 2:
            r = 0
        else:
            r = compute_R(p, primes, is_prime)
        
        a = math.isqrt(p)
        C = a * (a + 1)
        b = C - p
        c = p - a * a
        K = a + c
        
        constant = r * (math.log(p)) ** 2 / p if p > 2 else 0
        gcd_val = math.gcd(abs(b), C) if b != 0 else C
        
        results.append({
            'p': p,
            'a': a,
            'a_squared': a * a,
            'c': c,
            'center': C,
            'b': b,
            'R_p': r,
            'R_p_abs': abs(b),
            'K': K,
            'constant': constant,
            'gcd_b_center': gcd_val,
            'is_escape': gcd_val == 1,
            'four_K_minus_1': 4 * K - 1,
            'layer_size': 2 * a + 1
        })
        
        if (idx + 1) % 100000 == 0:
            elapsed = time.time() - t1
            print(f"   Progress: {idx+1:,}/{len(primes):,} ({elapsed:.1f}s)")
    
    t2 = time.time()
    print(f"   Completed in {t2-t1:.1f}s")
    
    return results


def save_csv(results, filename='R_p_results.csv'):
    """Save results to CSV."""
    fieldnames = list(results[0].keys())
    
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"✅ Saved {len(results):,} rows to {filename}")


def print_summary(results):
    """Print summary statistics."""
    constants = [r['constant'] for r in results if r['constant'] > 0]
    rp_values = [r['R_p'] for r in results]
    
    print()
    print("=" * 70)
    print("📊 Summary Statistics")
    print("=" * 70)
    print(f"   Total primes: {len(results):,}")
    print(f"   Mean R(p): {np.mean(rp_values):.4f}")
    print(f"   Max R(p): {max(rp_values)}")
    print(f"   Mean constant: {np.mean(constants):.6f}")
    print(f"   Std constant: {np.std(constants):.6f}")
    print(f"   Median constant: {np.median(constants):.6f}")
    
    # By ranges
    ranges = [(1, 1000), (1000, 10000), (10000, 100000), (100000, 1000000), (1000000, 10000000)]
    for lo, hi in ranges:
        vals = [r['constant'] for r in results if lo <= r['p'] < hi and r['constant'] > 0]
        if vals:
            print(f"   Range [{lo:>8,}, {hi:>8,}): mean={np.mean(vals):.6f}, std={np.std(vals):.6f}, n={len(vals)}")


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
    results = generate_data(limit)
    save_csv(results, f'R_p_results_{limit}.csv')
    print_summary(results)