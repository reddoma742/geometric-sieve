"""
singular_series.py
------------------
Compares the empirical 0.73 constant against the Hardy-Littlewood
singular series and the twin prime constant C2.

The representation function R(p) counts primes a < p such that
b = 2a - p is 1 or prime. This is a classic two-linear-forms problem:
    L1(a) = a,   L2(a) = 2a - p

The singular series for this pair converges to 2*C2, where C2 is the
twin prime constant (≈ 0.66016). The crude normalization
R(p)*(ln p)^2 / p is expected to tend toward C2 (≈ 0.66016), not 0.73,
with the observed 0.73 being a finite-range effective value.
"""

import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

# ============================================================
# 1. Load data
# ============================================================
filename = "R_p_results_10M.csv"

# Try to find the file
if not os.path.exists(filename):
    # Try in data/ folder
    if os.path.exists(os.path.join("data", filename)):
        filename = os.path.join("data", filename)
    else:
        print(f"File {filename} not found in current directory or data/ folder.")
        print("Please run verify_Rp.py first, or place the file in the data/ folder.")
        exit(1)

df = pd.read_csv(filename)
p = df['p'].values.astype(int)
R = df['R_p'].values.astype(int)

print(f"Loaded {len(p)} primes from {filename}.")

# ============================================================
# 2. Twin prime constant
# ============================================================
C2 = 0.660161815846869573927812110014555778432623360284733413319
two_C2 = 2 * C2

print(f"\nTwin prime constant C2      = {C2:.15f}")
print(f"2 * C2                      = {two_C2:.15f}")

# ============================================================
# 3. Helper: is (p+1)/2 prime? (safe prime check)
# ============================================================
def is_prime_naive(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

# ============================================================
# 4. Compute R*(p) and I*(p) for a sample
# ============================================================
# Take a sample for speed (every N primes + the last 200)
step = max(1, len(p) // 1000)
indices = np.arange(0, len(p), step)
indices = np.unique(np.concatenate([indices, np.arange(max(0, len(p)-200), len(p))]))

p_samp = p[indices]
R_samp = R[indices]

R_star = []
I_star = []

print(f"\nComputing I*(p) for {len(p_samp)} sample primes...")

for i, p_val in enumerate(p_samp):
    # delta(p) = 1 if (p+1)/2 is prime (safe prime)
    a1 = (p_val + 1) // 2
    delta = 1 if is_prime_naive(a1) else 0
    R_star.append(R_samp[i] - delta)

    # I*(p) = integral from (p+3)/2 to p of dt / (ln t * ln(2t-p))
    lower = (p_val + 3) / 2.0
    upper = p_val - 1.0
    if lower >= upper:
        I_star.append(0.0)
        continue

    def integrand(t):
        return 1.0 / (np.log(t) * np.log(2*t - p_val))

    res, _ = quad(integrand, lower, upper, limit=200, epsabs=1e-10, epsrel=1e-10)
    I_star.append(res)

    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{len(p_samp)} done")

R_star = np.array(R_star)
I_star = np.array(I_star)

# ============================================================
# 5. Compute Q*(p) = R*(p) / I*(p)
# ============================================================
mask = I_star > 0
Q_star = np.full_like(I_star, np.nan)
Q_star[mask] = R_star[mask] / I_star[mask]

valid = mask & ~np.isnan(Q_star)
Q_valid = Q_star[valid]
p_valid = p_samp[valid]

# ============================================================
# 6. Compute C_tilde(p) = R*(p)*(ln p)^2 / p
# ============================================================
log_p = np.log(p_valid)
C_tilde = R_star[valid] * (log_p**2) / p_valid

# ============================================================
# 7. Statistics
# ============================================================
mean_Q_all = np.mean(Q_valid)
mean_Q_last100 = np.mean(Q_valid[-100:]) if len(Q_valid) >= 100 else mean_Q_all
mean_C_all = np.mean(C_tilde)
mean_C_last100 = np.mean(C_tilde[-100:]) if len(C_tilde) >= 100 else mean_C_all

print(f"\n{'='*60}")
print(f"RESULTS")
print(f"{'='*60}")
print(f"Mean Q*(p) overall           : {mean_Q_all:.6f}")
print(f"Mean Q*(p) last 100 points   : {mean_Q_last100:.6f}")
print(f"Expected (2*C2)              : {two_C2:.6f}")
print(f"Deviation (last 100)         : {abs(mean_Q_last100 - two_C2)/two_C2*100:.2f}%")
print(f"")
print(f"Mean C_tilde(p) overall      : {mean_C_all:.6f}")
print(f"Mean C_tilde(p) last 100     : {mean_C_last100:.6f}")
print(f"Twin prime constant C2       : {C2:.6f}")
print(f"Observed 0.73 value          : 0.730000")

# ============================================================
# 8. Convergence trend
# ============================================================
bins = [1e3, 1e4, 1e5, 5e5, 1e6, 5e6, 1e7]
print(f"\n{'='*60}")
print(f"CONVERGENCE OF C_tilde(p) BY RANGE")
print(f"{'='*60}")
for i in range(len(bins)-1):
    mask_bin = (p_valid >= bins[i]) & (p_valid < bins[i+1])
    if mask_bin.sum() > 5:
        mean_bin = np.mean(C_tilde[mask_bin])
        print(f"  p in [{bins[i]:.0e}, {bins[i+1]:.0e}): "
              f"mean = {mean_bin:.6f}  (n = {mask_bin.sum()})")

# ============================================================
# 9. PLOTS — SAVE TO FILE
# ============================================================
# Create figures directory if it doesn't exist
os.makedirs("figures", exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# --- Plot 1: Q*(p) vs p ---
ax = axes[0, 0]
ax.plot(p_valid, Q_valid, '.', alpha=0.4, markersize=2, color='blue')
ax.axhline(y=two_C2, color='red', linestyle='--', linewidth=2,
           label=f'2·C₂ = {two_C2:.6f}')
ax.set_xscale('log')
ax.set_xlabel('p')
ax.set_ylabel('Q*(p) = R*(p) / I*(p)')
ax.set_title('Q*(p) vs p (Heuristic: Q*(p) → 2·C₂)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 2: Binned Q*(p) ---
ax = axes[0, 1]
log_bins = np.logspace(np.log10(p_valid.min()), np.log10(p_valid.max()), 25)
bin_centers = []
bin_means = []
for i in range(len(log_bins)-1):
    bm = (p_valid >= log_bins[i]) & (p_valid < log_bins[i+1])
    if bm.sum() > 5:
        bin_centers.append(np.sqrt(log_bins[i] * log_bins[i+1]))
        bin_means.append(np.mean(Q_valid[bm]))
ax.plot(bin_centers, bin_means, 'o-', color='darkblue', markersize=4,
        label='Binned mean Q*(p)')
ax.axhline(y=two_C2, color='red', linestyle='--', linewidth=2,
           label=f'2·C₂ = {two_C2:.6f}')
ax.set_xscale('log')
ax.set_xlabel('p (bin center)')
ax.set_ylabel('Mean Q*(p)')
ax.set_title('Binned Q*(p) → 2·C₂')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 3: C_tilde(p) vs p ---
ax = axes[1, 0]
ax.plot(p_valid, C_tilde, '.', alpha=0.4, markersize=2, color='green')
ax.axhline(y=C2, color='red', linestyle='--', linewidth=2,
           label=f'C₂ = {C2:.6f}')
ax.axhline(y=0.73, color='orange', linestyle=':', linewidth=2,
           label='0.73 (observed)')
ax.set_xscale('log')
ax.set_xlabel('p')
ax.set_ylabel('C_tilde(p) = R*(p)·(ln p)² / p')
ax.set_title('C_tilde(p) vs p (Expected limit: C₂)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 4: Convergence of C_tilde(p) by range ---
ax = axes[1, 1]
range_labels = []
range_means = []
for i in range(len(bins)-1):
    mask_bin = (p_valid >= bins[i]) & (p_valid < bins[i+1])
    if mask_bin.sum() > 5:
        range_labels.append(f'{bins[i]:.0e}\nto\n{bins[i+1]:.0e}')
        range_means.append(np.mean(C_tilde[mask_bin]))
x_pos = range(len(range_means))
ax.bar(x_pos, range_means, color='steelblue', edgecolor='navy')
ax.axhline(y=C2, color='red', linestyle='--', linewidth=2,
           label=f'C₂ = {C2:.6f}')
ax.set_xticks(x_pos)
ax.set_xticklabels(range_labels, fontsize=8)
ax.set_ylabel('Mean C_tilde(p)')
ax.set_title('Convergence of C_tilde(p) by range')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

plt.suptitle('Singular Series Analysis — Geometric Sieve', fontsize=14, fontweight='bold')
plt.tight_layout()

# SAVE THE FIGURE
output_path = os.path.join("figures", "singular_series_analysis.png")
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"\nFigure saved to: {output_path}")
print(f"\nDone.")
