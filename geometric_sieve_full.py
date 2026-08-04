"""
geometric_sieve_full.py
Complete Geometric Sieve system with all components.

Implements:
1. Layer system: a = floor(sqrt(n)), C(a) = a(a+1), b = C(a) - n
2. Geometric Sieve (Theorem 1): gcd(|b|, C(a)) > 1 ⇒ composite
3. R(p) computation
4. Golden K values and Heegner arcs
5. Escape position analysis

Usage: python geometric_sieve_full.py
"""

import math
import numpy as np
from collections import defaultdict, Counter
import csv


class GeometricSieve:
    """Complete Geometric Sieve system."""
    
    def __init__(self, limit=1_000_000):
        self.limit = limit
        self.primes_set = set()
        self.primes_list = []
        self._generate_primes()
    
    def _generate_primes(self):
        """Generate primes using Sieve of Eratosthenes."""
        is_prime = np.ones(self.limit + 1, dtype=bool)
        is_prime[:2] = False
        for i in range(2, int(self.limit ** 0.5) + 1):
            if is_prime[i]:
                is_prime[i * i::i] = False
        self.primes_list = np.nonzero(is_prime)[0].tolist()
        self.primes_set = set(self.primes_list)
    
    def get_layer(self, n):
        """Get layer a = floor(sqrt(n))."""
        return math.isqrt(n)
    
    def get_center(self, a):
        """Get geometric center C(a) = a(a+1)."""
        return a * (a + 1)
    
    def get_b(self, n, a=None):
        """Get coordinate b = C(a) - n."""
        if a is None:
            a = self.get_layer(n)
        return self.get_center(a) - n
    
    def get_c(self, n, a=None):
        """Get position c = n - a^2."""
        if a is None:
            a = self.get_layer(n)
        return n - a * a
    
    def is_structural_composite(self, n):
        """Theorem 1: gcd(|b|, C(a)) > 1 ⇒ composite."""
        if n <= 3:
            return n == 0 or n == 1 or n == 4  # 4 is composite
        
        a = self.get_layer(n)
        C = self.get_center(a)
        b = C - n
        d = math.gcd(abs(b), C)
        return d > 1
    
    def is_escape_position(self, n):
        """Check if n is an escape position."""
        if n <= 3:
            return n in {2, 3}
        return not self.is_structural_composite(n)
    
    def sieve_efficiency(self, limit=None):
        """Measure Geometric Sieve efficiency."""
        if limit is None:
            limit = self.limit
        
        total_composites = 0
        caught = 0
        false_positives = 0
        
        for n in range(4, limit + 1):
            flagged = self.is_structural_composite(n)
            actually = n not in self.primes_set
            
            if actually:
                total_composites += 1
                if flagged:
                    caught += 1
            else:
                if flagged:
                    false_positives += 1
        
        return total_composites, caught, false_positives
    
    def compute_R(self, p):
        """Compute R(p) = number of representations p = 2a - b."""
        if p <= 2:
            return 0
        
        cnt = 0
        for a in self.primes_list:
            if a >= p:
                break
            b = 2 * a - p
            if b >= 1 and (b == 1 or b in self.primes_set):
                cnt += 1
        return cnt
    
    def get_golden_K_values(self):
        """Return known golden K values (Heegner-related)."""
        return {2, 3, 5, 11, 17, 41}
    
    def scan_K_arcs(self, max_K=200, max_a=1000):
        """Scan for prime-generating arcs n = a^2 - a + K."""
        golden_K = self.get_golden_K_values()
        results = []
        
        for K in range(1, max_K + 1):
            prime_count = 0
            first_failure = None
            last_prime = None
            
            for a in range(1, max_a + 1):
                n = a * a - a + K
                if n > self.limit:
                    break
                if n in self.primes_set:
                    prime_count += 1
                    last_prime = n
                else:
                    if first_failure is None:
                        first_failure = n
                    if prime_count > 40:
                        break
                    if a > prime_count + 5:
                        break
            
            results.append({
                'K': K,
                'arc_length': prime_count,
                'last_prime': last_prime,
                'first_failure': first_failure,
                'is_golden': K in golden_K,
                'four_K_minus_1': 4 * K - 1,
                'four_K_minus_1_prime': (4 * K - 1) in self.primes_set
            })
        
        return results
    
    def analyze_escape_positions(self, max_a=50):
        """Analyze escape positions in each layer."""
        results = []
        
        for a in range(1, max_a + 1):
            a_sq = a * a
            next_sq = (a + 1) * (a + 1)
            C = self.get_center(a)
            
            escape_count = 0
            prime_in_escape = 0
            composite_in_escape = 0
            structural_count = 0
            
            for n in range(a_sq, min(next_sq, self.limit + 1)):
                if n <= 1:
                    continue
                
                is_escape = self.is_escape_position(n)
                is_prime = n in self.primes_set
                is_structural = self.is_structural_composite(n)
                
                if is_structural:
                    structural_count += 1
                elif is_escape:
                    escape_count += 1
                    if is_prime:
                        prime_in_escape += 1
                    else:
                        composite_in_escape += 1
            
            results.append({
                'a': a,
                'layer_size': 2 * a + 1,
                'center': C,
                'structural_composites': structural_count,
                'escape_positions': escape_count,
                'primes_in_escape': prime_in_escape,
                'external_composites': composite_in_escape
            })
        
        return results


def main():
    print("=" * 70)
    print("🔬 Geometric Sieve - Complete System")
    print("=" * 70)
    print()
    
    # Initialize
    limit = 100_000
    print(f"⏳ Initializing system up to {limit:,}...")
    gs = GeometricSieve(limit)
    print(f"   Generated {len(gs.primes_list):,} primes")
    print()
    
    # Test 1: Sieve Efficiency
    print("=" * 70)
    print("📊 Test 1: Geometric Sieve Efficiency")
    print("=" * 70)
    
    for test_limit in [1_000, 10_000, 100_000]:
        tc, c, fp = gs.sieve_efficiency(test_limit)
        eff = 100 * c / tc if tc > 0 else 0
        print(f"   Limit {test_limit:>8,}: {eff:>6.2f}% efficiency, {fp} false positives")
    print()
    
    # Test 2: R(p) Constant
    print("=" * 70)
    print("📊 Test 2: R(p) Constant (Conjecture 1)")
    print("=" * 70)
    
    sample_primes = [p for p in gs.primes_list if 1000 <= p <= 50000][:1000]
    rp_values = []
    
    for p in sample_primes:
        r = gs.compute_R(p)
        rp_values.append(r * (math.log(p)) ** 2 / p)
    
    mean_rp = np.mean(rp_values)
    std_rp = np.std(rp_values)
    print(f"   Sample: {len(sample_primes)} primes in [1000, 50000)")
    print(f"   Mean R(p)*(ln p)^2/p: {mean_rp:.6f}")
    print(f"   Std: {std_rp:.6f}")
    print()
    
    # Test 3: Golden K values
    print("=" * 70)
    print("📊 Test 3: Heegner-type Golden K Arcs")
    print("=" * 70)
    
    arcs = gs.scan_K_arcs(max_K=50, max_a=500)
    arcs_sorted = sorted(arcs, key=lambda x: x['arc_length'], reverse=True)
    
    print(f"   {'K':>4} {'Arc Length':>12} {'Last Prime':>10} {'4K-1':>6} {'4K-1 Prime':>12} {'Golden':>8}")
    print("   " + "-" * 60)
    for arc in arcs_sorted[:10]:
        marker = " ⭐" if arc['is_golden'] else ""
        print(f"   {arc['K']:>4} {arc['arc_length']:>12} {arc['last_prime'] or 'N/A':>10} "
              f"{arc['four_K_minus_1']:>6} {str(arc['four_K_minus_1_prime']):>12} {marker:>8}")
    print()
    
    # Test 4: Escape positions
    print("=" * 70)
    print("📊 Test 4: Escape Position Analysis (First 10 layers)")
    print("=" * 70)
    
    escape_data = gs.analyze_escape_positions(10)
    print(f"   {'a':>3} {'Size':>5} {'Structural':>10} {'Escape':>8} {'Primes':>8} {'External':>9}")
    print("   " + "-" * 50)
    for row in escape_data:
        print(f"   {row['a']:>3} {row['layer_size']:>5} {row['structural_composites']:>10} "
              f"{row['escape_positions']:>8} {row['primes_in_escape']:>8} {row['external_composites']:>9}")
    
    print()
    print("=" * 70)
    print("✅ All tests complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()