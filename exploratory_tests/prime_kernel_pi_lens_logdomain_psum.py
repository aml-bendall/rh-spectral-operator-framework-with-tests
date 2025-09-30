#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Real-line Bendall lens (corrected):
 - Work in u = log x
 - Even, band-limited kernels (Fejer / Gaussian) with PROPER NORMALIZATION so integral K(u) du = 1
 - Include prime powers in psi-hat (Lambda(n) = log p for n = p^k)
 - pi-hat via discrete partial summation: pi_hat(x) ~= sum_{n<=X} [Lambda(n)/log n] * K_norm( log x - log n )

Usage example:
  python prime_kernel_pi_lens_logdomain_psum.py --x 1e3 1e4 1e5 1e6     --Nmax 2000000 --kernel fejer --Delta_u 0.05 --lam 0.05 --out RL2b_corrected_psum.csv
"""
import argparse, math, csv, bisect

# ---------- primes and prime powers ----------
def sieve(n: int):
    bs=[True]*(n+1); bs[:2]=[False,False]
    for p in range(2, int(n**0.5)+1):
        if bs[p]:
            step=p; start=p*p
            bs[start:n+1:step]=[False]*(((n-start)//step)+1)
    return [i for i,b in enumerate(bs) if b]

def prime_powers(Nmax, primes):
    """Generate all prime powers n = p^k <= Nmax with weight Lambda(n) = log p.
    Returns a list of (n, logn, lambdan) sorted by n."""
    out = []
    for p in primes:
        v = p
        lp = math.log(p)
        while v <= Nmax:
            out.append((v, math.log(v), lp))
            v *= p
    out.sort(key=lambda t: t[0])
    return out

# ---------- special functions ----------
def li(x: float) -> float:
    import mpmath as mp
    mp.mp.dps = 50
    return float(mp.li(x))

def mobius(n: int) -> int:
    m=n; p=2; cnt=0
    while p*p<=m:
        if m%p==0:
            m//=p; cnt+=1
            if m%p==0: return 0
        p += 1 if p==2 else 2
    if m>1: cnt+=1
    return -1 if (cnt%2) else 1

def R_of_x(x: float, K: int = 20) -> float:
    s = 0.0
    for k in range(1, K+1):
        mu = mobius(k)
        if mu == 0: continue
        s += mu/k * li(x**(1.0/k))
    return s

# ---------- normalized even kernels in u = log x ----------
# Ensure integral K_norm(u) du = 1

def fejer_norm(u: float, Delta_u: float) -> float:
    # Base Fejer: K(u) = sinc^2( pi u / (2 Delta_u) ), integral over u is Delta_u
    # Normalize by 1/Delta_u
    if u == 0.0:
        return 1.0 / Delta_u
    z = math.pi * u / (2.0 * Delta_u)
    s = math.sin(z) / z
    return (s*s) / Delta_u

def gauss_norm(u: float, Delta_u: float) -> float:
    # Choose sigma = Delta_u / sqrt(2 ln 2): FWHM ~ 2*Delta_u
    sigma = Delta_u / math.sqrt(2.0*math.log(2.0))
    # Normalized Gaussian: (1/(sigma sqrt(2*pi))) * exp(-u^2/(2 sigma^2))
    return math.exp(-0.5*(u/sigma)**2) / (sigma * math.sqrt(2.0*math.pi))

def kernel_norm(u: float, Delta_u: float, kind: str) -> float:
    return fejer_norm(u, Delta_u) if kind=="fejer" else gauss_norm(u, Delta_u)

# ---------- exact pi(x) by sieve ----------
def pi_exact_of_x(x: int, primes):
    import bisect
    return bisect.bisect_right(primes, x)

# ---------- lens estimators ----------
def psi_hat(x: float, pp_list, Delta_u: float, kernel_kind: str, lam: float) -> float:
    """psi_hat(x) = sum_{n <= Nmax} Lambda(n) * K_norm( log x - log n ; Delta_u ) / (1+lam)"""
    u = math.log(x)
    acc = 0.0
    for n, ln, lam_n in pp_list:
        if n > x: break
        acc += lam_n * kernel_norm(u - ln, Delta_u, kernel_kind)
    return acc / (1.0 + lam)

def pi_hat_psum(x: float, pp_list, Delta_u: float, kernel_kind: str, lam: float) -> float:
    """pi_hat(x) via discrete partial summation:
       pi_hat(x) ~= sum_{n <= x} [Lambda(n)/log n] * K_norm( log x - log n ; Delta_u ) / (1+lam)"""
    u = math.log(x)
    acc = 0.0
    for n, ln, lam_n in pp_list:
        if n > x: break
        acc += (lam_n / ln) * kernel_norm(u - ln, Delta_u, kernel_kind)
    return acc / (1.0 + lam)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", nargs="+", type=float, required=True, help="evaluation points (e.g., 1e3 1e4 1e5 1e6)")
    ap.add_argument("--Nmax", type=int, default=2000000, help="max n for prime powers")
    ap.add_argument("--kernel", choices=["fejer","gauss"], default="fejer")
    ap.add_argument("--Delta_u", type=float, default=0.05)
    ap.add_argument("--lam", type=float, default=0.05)
    ap.add_argument("--out", type=str, default="RL2b_corrected_psum.csv")
    args = ap.parse_args()

    xmax = int(max(args.x))
    Nmax = max(args.Nmax, xmax)
    primes = sieve(Nmax)
    pp_list = prime_powers(Nmax, primes)

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x","psi_hat","pi_hat","pi_exact","Li","R"])
        for x in args.x:
            psi = psi_hat(x, pp_list, args.Delta_u, args.kernel, args.lam)
            piH = pi_hat_psum(x, pp_list, args.Delta_u, args.kernel, args.lam)
            pix = pi_exact_of_x(int(x), primes)
            lix = li(x)
            Rx  = R_of_x(x)
            w.writerow([int(x), psi, piH, pix, lix, Rx])

if __name__ == "__main__":
    main()
