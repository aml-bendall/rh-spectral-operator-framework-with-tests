#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RL2b sweep with calibration:
 - Real-line Bendall lens in u=log x (Fejer/Gaussian), normalized kernels
 - Includes prime powers and discrete partial summation to pi(x)
 - Sweeps kernels in {fejer, gauss} and Delta_u in {0.03, 0.05}
 - Calibrates a global scale (or affine) factor on a small-x set, then evaluates all x

Usage:
  python3 run_rl2b_sweep_calibrated.py --Nmax 2000000 --lam 0.05 --out RL2b_sweep_cal.csv       --fit scale --calx 1000 10000
  (use --fit affine to fit y ≈ c * yhat + b)
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

def fit_scale(yhat_list, y_list):
    num = sum(yh*y for yh,y in zip(yhat_list, y_list))
    den = sum(yh*yh for yh in yhat_list)
    return (num/den) if den>0 else 1.0

def fit_affine(yhat_list, y_list):
    n = len(yhat_list)
    Sx  = sum(yhat_list)
    Sy  = sum(y_list)
    Sxx = sum(yh*yh for yh in yhat_list)
    Sxy = sum(yh*y  for yh,y in zip(yhat_list,y_list))
    den = n*Sxx - Sx*Sx
    if abs(den) < 1e-12:
        return 1.0, 0.0
    c = (n*Sxy - Sx*Sy)/den
    b = (Sy - c*Sx)/n
    return c, b

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Nmax", type=int, default=2000000)
    ap.add_argument("--lam", type=float, default=0.05)
    ap.add_argument("--out", type=str, default="RL2b_sweep_cal.csv")
    ap.add_argument("--fit", choices=["scale","affine"], default="scale")
    ap.add_argument("--calx", nargs="+", type=int, default=[1000,10000], help="x points to fit on")
    args = ap.parse_args()

    xs = [10**3, 10**4, 10**5, 10**6]
    grids = [("fejer", 0.03), ("fejer", 0.05), ("gauss", 0.03), ("gauss", 0.05)]

    primes = sieve(args.Nmax)
    pp_list = prime_powers(args.Nmax, primes)

    rows = []
    for kernel, Du in grids:
        # First, compute raw predictions on both calibration xs and eval xs
        yhat_cal, y_cal = [], []
        cache = {}
        for x in xs:
            psi = psi_hat(x, pp_list, Du, kernel, args.lam)
            piH = pi_hat_psum(x, pp_list, Du, kernel, args.lam)
            pix = pi_exact_of_x(x, primes)
            lix = li(x)
            Rx  = R_of_x(x)
            cache[x] = (psi, piH, pix, lix, Rx)
            if x in args.calx:
                yhat_cal.append(piH); y_cal.append(pix)

        # Fit scale or affine on calibration points
        if args.fit == "scale":
            c = fit_scale(yhat_cal, y_cal)
            b = 0.0
        else:
            c, b = fit_affine(yhat_cal, y_cal)

        # Emit calibrated results
        for x in xs:
            psi, piH, pix, lix, Rx = cache[x]
            piHc = c*piH + b
            rows.append({
                "kernel": kernel, "Delta_u": Du, "lam": args.lam,
                "fit": args.fit, "c": c, "b": b, "calset": ",".join(map(str,args.calx)),
                "x": x, "psi_hat": psi, "pi_hat_raw": piH, "pi_hat_cal": piHc,
                "pi_exact": pix, "Li": lix, "R": Rx,
                "pi_err_raw": piH - pix, "pi_err_cal": piHc - pix,
                "Li_err": lix - pix, "R_err": Rx - pix
            })

    # save CSV
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows: w.writerow(r)

    # print compact tables per grid
    for kernel, Du in grids:
        print(f"\n=== {kernel}, Delta_u={Du}, lam={args.lam} (fit={args.fit} on {args.calx}) ===")
        print("x      |  pi_hat_cal |  pi(x)  |  err_cal ||  R(x)  | R-err  |  Li(x)  | Li-err  ||  c      b")
        for r in rows:
            if r["kernel"]==kernel and r["Delta_u"]==Du:
                x, piHc, pix, Rx, lix = r["x"], r["pi_hat_cal"], r["pi_exact"], r["R"], r["Li"]
                print(f"{x:<7}| {piHc:11.3f}| {pix:7.0f}| {piHc-pix:+8.3f} || {Rx:7.3f}| {Rx-pix:+7.3f} | {lix:7.3f}| {lix-pix:+7.3f} || {r['c']:.6f} {r['b']:.3f}")

if __name__ == "__main__":
    main()
