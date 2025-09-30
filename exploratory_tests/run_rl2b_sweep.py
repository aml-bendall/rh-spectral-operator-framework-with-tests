#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RL2b sweep (standalone):
  - Real-line Bendall lens in u=log x with normalized even kernels (Fejer/Gaussian)
  - Includes prime powers (Lambda(n)=log p for n=p^k)
  - Discrete partial-summation bridge to pi(x)
  - Sweeps kernels in {fejer, gauss} and Delta_u in {0.03, 0.05}
  - Evaluates at x in {1e3, 1e4, 1e5, 1e6} (modifiable)
  - Outputs a CSV and prints a compact comparison table

Usage:
  python3 run_rl2b_sweep.py --Nmax 2000000 --lam 0.05 --out RL2b_sweep.csv
"""
import argparse, math, csv, bisect
from math import log, sqrt, pi, sin

def sieve(n: int):
    bs=[True]*(n+1); bs[:2]=[False,False]
    for p in range(2, int(n**0.5)+1):
        if bs[p]:
            step=p; start=p*p
            bs[start:n+1:step]=[False]*(((n-start)//step)+1)
    return [i for i,b in enumerate(bs) if b]

def prime_powers(Nmax, primes):
    out = []
    for p in primes:
        v = p
        lp = math.log(p)
        while v <= Nmax:
            out.append((v, math.log(v), lp))
            v *= p
    out.sort(key=lambda t: t[0])
    return out

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

# Normalized kernels: integral in u equals 1
def fejer_norm(u: float, Delta_u: float) -> float:
    if u == 0.0:
        return 1.0 / Delta_u
    z = pi * u / (2.0 * Delta_u)
    s = sin(z) / z
    return (s*s) / Delta_u

def gauss_norm(u: float, Delta_u: float) -> float:
    sigma = Delta_u / math.sqrt(2.0*math.log(2.0))  # FWHM ~ 2*Delta_u
    return math.exp(-0.5*(u/sigma)**2) / (sigma * math.sqrt(2.0*math.pi))

def K_norm(u: float, Delta_u: float, kind: str) -> float:
    return fejer_norm(u, Delta_u) if kind == "fejer" else gauss_norm(u, Delta_u)

def pi_exact_of_x(x: int, primes):
    return bisect.bisect_right(primes, x)

def psi_hat(x: float, pp_list, Delta_u: float, kernel_kind: str, lam: float) -> float:
    u = math.log(x)
    acc = 0.0
    for n, ln, lam_n in pp_list:
        if n > x: break
        acc += lam_n * K_norm(u - ln, Delta_u, kernel_kind)
    return acc / (1.0 + lam)

def pi_hat_psum(x: float, pp_list, Delta_u: float, kernel_kind: str, lam: float) -> float:
    u = math.log(x)
    acc = 0.0
    for n, ln, lam_n in pp_list:
        if n > x: break
        acc += (lam_n / ln) * K_norm(u - ln, Delta_u, kernel_kind)
    return acc / (1.0 + lam)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Nmax", type=int, default=2000000)
    ap.add_argument("--lam", type=float, default=0.05)
    ap.add_argument("--out", type=str, default="RL2b_sweep.csv")
    args = ap.parse_args()

    xs = [10**3, 10**4, 10**5, 10**6]
    grids = [("fejer", 0.03), ("fejer", 0.05), ("gauss", 0.03), ("gauss", 0.05)]

    primes = sieve(args.Nmax)
    pp_list = prime_powers(args.Nmax, primes)

    rows = []
    for kernel, Du in grids:
        for x in xs:
            psi = psi_hat(x, pp_list, Du, kernel, args.lam)
            piH = pi_hat_psum(x, pp_list, Du, kernel, args.lam)
            pix = pi_exact_of_x(x, primes)
            lix = li(x)
            Rx  = R_of_x(x)
            rows.append({
                "kernel": kernel, "Delta_u": Du, "lam": args.lam,
                "x": x, "psi_hat": psi, "pi_hat": piH,
                "pi_exact": pix, "Li": lix, "R": Rx,
                "pi_err": piH - pix, "Li_err": lix - pix, "R_err": Rx - pix
            })

    # save CSV
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows: w.writerow(r)

    # print compact tables per grid
    for kernel, Du in grids:
        print(f"\n=== {kernel}, Delta_u={Du}, lam={args.lam} ===")
        print("x      |  pi_hat   |  pi(x)  |  err    ||   R(x)   | R-err   |   Li(x)  | Li-err")
        for r in rows:
            if r["kernel"]==kernel and r["Delta_u"]==Du:
                x, piH, pix, Rx, lix = r["x"], r["pi_hat"], r["pi_exact"], r["R"], r["Li"]
                print(f"{x:<7}| {piH:8.3f}| {pix:7.0f}| {piH-pix:+8.3f} || {Rx:8.3f}| {Rx-pix:+8.3f} | {lix:8.3f}| {lix-pix:+8.3f}")

if __name__ == "__main__":
    main()
