#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Real-line Bendall lens: log-domain, band-limited prime-kernel estimator
Computes ψ_hat(x) and π_hat(x) from primes ≤ Pmax using even kernels in u=log x.

Usage (example):
  python prime_kernel_pi_lens_logdomain.py --x 1000 10000 100000 1000000 \
      --Pmax 2000000 --kernel fejer --Delta_u 0.04 --lam 0.05 --out RL2b_fejer.csv

Kernel options:
  --kernel fejer | gauss
Bandwidth:
  --Delta_u (e.g., 0.03..0.10). Think “window half-width” in u=log x.
Regularization:
  --lam: small Tikhonov-like softening applied to weights (keeps numerics tame)

Outputs CSV with columns: x, psi_hat, pi_hat
"""

import argparse, math, csv
from math import log, sqrt, pi, erf, sin

# ---------- primes & Möbius ----------

def sieve(n: int):
    bs=[True]*(n+1); bs[:2]=[False,False]
    for p in range(2, int(n**0.5)+1):
        if bs[p]:
            step=p; start=p*p
            bs[start:n+1:step]=[False]*(((n-start)//step)+1)
    return [i for i,b in enumerate(bs) if b]

def mobius(n: int) -> int:
    # tiny utility for R(x); brute is fine for small n in benchmark
    m=n; p=2; cnt=0
    while p*p<=m:
        if m%p==0:
            m//=p; cnt+=1
            if m%p==0: return 0
        p += 1 if p==2 else 2
    if m>1: cnt+=1
    return -1 if (cnt%2) else 1

# ---------- special functions ----------

def li(x: float) -> float:
    # use mpmath for numerical li (offset logarithmic integral)
    import mpmath as mp
    mp.mp.dps = 50
    return float(mp.li(x))

def R_of_x(x: float, K: int = 20) -> float:
    s = 0.0
    for k in range(1, K+1):
        mu = mobius(k)
        if mu == 0: continue
        s += mu/k * li(x**(1.0/k))
    return s

# ---------- kernels in u = log x ----------

def fejer_kernel(u: float, Delta_u: float) -> float:
    # Continuous Fejér window in frequency-like coordinate:
    # K(u) ∝ (sin(π u / (2Δ)) / (π u / (2Δ)))^2, normalized to integrate ≈ 1.
    if u == 0.0:
        return 1.0
    z = pi * u / (2.0 * Delta_u)
    s = sin(z) / z
    return s * s  # Fejér ~ sinc^2

def gauss_kernel(u: float, Delta_u: float) -> float:
    # Gaussian with std related to Delta_u: choose σ = Δ_u / sqrt(2 ln 2) so FWHM ≈ 2Δ_u.
    sigma = Delta_u / sqrt(2.0*math.log(2.0))
    return math.exp(-0.5*(u/sigma)**2)

def even_kernel(u: float, Delta_u: float, kind: str) -> float:
    return fejer_kernel(u, Delta_u) if kind=="fejer" else gauss_kernel(u, Delta_u)

# ---------- lens estimators on ℝ⁺ ----------

def psi_hat(x: float, primes, Delta_u: float, kernel_kind: str, lam: float) -> float:
    """
    ψ̂(x) ≈ sum_{p ≤ Pmax} (log p) * K( log x - log p ; Δ_u ), even in u.
    lam enters as a soft global scale (1 / (1+lam)) to tame oversharp windows.
    """
    u = log(x)
    acc = 0.0
    for p in primes:
        if p > x: break  # optional cheap cut (can include all if desired)
        up = log(p)
        acc += log(p) * even_kernel(u - up, Delta_u, kernel_kind)
    return acc / (1.0 + lam)

def pi_hat_from_psi(x: float, primes, Delta_u: float, kernel_kind: str, lam: float) -> float:
    """
    Direct π̂(x) estimator using prime counting lens:
      π̂(x) ≈ sum_{p ≤ Pmax} K( log x - log p ; Δ_u )
    
    Each prime p contributes its kernel weight K(log x - log p) to the count.
    This gives a smoothed version of the prime counting function π(x).
    
    Note: This is NOT derived from ψ̂ by partial summation, but is a direct
    lens estimator for π(x) using the same kernel smoothing approach.
    """
    u = log(x)
    acc = 0.0
    for p in primes:
        if p > x: break  # only count primes ≤ x for π(x)
        up = log(p)
        # Direct kernel weighting: each prime contributes ~1 when close to x
        acc += even_kernel(u - up, Delta_u, kernel_kind)
    return acc / (1.0 + lam)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", nargs="+", type=float, required=True, help="evaluation points (e.g., 1e3 1e4 1e5 1e6)")
    ap.add_argument("--Pmax", type=int, default=2_000_000)
    ap.add_argument("--kernel", choices=["fejer","gauss"], default="fejer")
    ap.add_argument("--Delta_u", type=float, default=0.05)
    ap.add_argument("--lam", type=float, default=0.05)
    ap.add_argument("--out", type=str, default="RL2b_out.csv")
    args = ap.parse_args()

    xmax = int(max(args.x))
    Pmax = max(args.Pmax, xmax)
    primes = sieve(Pmax)

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x","psi_hat","pi_hat","pi_exact","Li","R"])
        for x in args.x:
            psi = psi_hat(x, primes, args.Delta_u, args.kernel, args.lam)
            pi_est = pi_hat_from_psi(x, primes, args.Delta_u, args.kernel, args.lam)
            # baselines for context
            # exact pi(x):
            import bisect
            pi_exact = bisect.bisect_right(primes, int(x))
            lix = li(x)
            Rx  = R_of_x(x)
            w.writerow([int(x), psi, pi_est, pi_exact, lix, Rx])

if __name__ == "__main__":
    main()
