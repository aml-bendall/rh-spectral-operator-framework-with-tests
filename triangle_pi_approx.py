#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Triangle-gap approximation for the prime counting function pi(x).

Idea:
  Between consecutive primes p_k and p_{k+1}, approximate the staircase linearly:
     pi_tilde(x) = (k-1) + (x - p_k) / (p_{k+1} - p_k),    for p_k <= x < p_{k+1}
  where k = number of primes <= p_k.
  This makes the local slope 1 / gap, i.e., exactly the "right triangle" picture.

Outputs a CSV with columns:
  x, pi_exact, pi_triangle, tri_err, Li, Li_err, R, R_err

Optional: --plot overlays pi(x), triangle, Li, R on sample points using matplotlib (no seaborn).

Usage:
  python triangle_pi_approx.py --x 1000 10000 100000 1000000 --out triangle_vs_pi.csv --plot
"""

import argparse, math, csv, bisect

def sieve(n: int):
    bs=[True]*(n+1); bs[:2]=[False,False]
    for p in range(2, int(n**0.5)+1):
        if bs[p]:
            step=p; start=p*p
            bs[start:n+1:step]=[False]*(((n-start)//step)+1)
    return [i for i,b in enumerate(bs) if b]

def mobius(n: int) -> int:
    m=n; p=2; cnt=0
    while p*p<=m:
        if m%p==0:
            m//=p; cnt+=1
            if m%p==0: return 0
        p += 1 if p==2 else 2
    if m>1: cnt+=1
    return -1 if (cnt%2) else 1

def li(x: float) -> float:
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

def pi_exact(x: int, primes) -> int:
    return bisect.bisect_right(primes, x)

def pi_tilde_triangle(x: int, primes):
    k = pi_exact(x, primes)  # primes <= x
    if k == 0:
        return 0.0
    if k >= len(primes):
        # extrapolate beyond last tabulated prime using avg of last few gaps
        p_prev = primes[-1]
        gaps = [primes[i+1]-primes[i] for i in range(max(0,len(primes)-10), len(primes)-1)]
        avg_gap = sum(gaps)/len(gaps)
        return (len(primes) - 1) + max(0, (x - p_prev)/avg_gap)
    p_k = primes[k-1]
    p_next = primes[k]
    gap = p_next - p_k
    frac = (x - p_k) / gap
    return (k-1) + frac

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", nargs="+", type=int, required=True, help="evaluation points (e.g., 1000 10000 100000 1000000)")
    ap.add_argument("--out", type=str, default="triangle_vs_pi.csv")
    ap.add_argument("--plot", action="store_true", help="make a matplotlib plot of points")
    args = ap.parse_args()

    xmax = max(args.x)
    primes = sieve(xmax + 1000)

    rows = []
    for x in args.x:
        piE = pi_exact(x, primes)
        pit = pi_tilde_triangle(x, primes)
        lix = li(x)
        Rx  = R_of_x(x)
        rows.append(dict(
            x=x,
            pi_exact=piE,
            pi_triangle=pit,
            tri_err=pit - piE,
            Li=lix,
            Li_err=lix - piE,
            R=Rx,
            R_err=Rx - piE,
        ))

    # Save CSV
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows: w.writerow(r)

    print(f"Wrote {args.out} with {len(rows)} rows.")
    for r in rows:
        print(f"x={r['x']:<8} pi={r['pi_exact']:<8} tri={r['pi_triangle']:<12.6f}  err={r['tri_err']:+.6f}"
              f"  |  R={r['R']:.3f} (R-err={r['R_err']:+.3f})  Li={r['Li']:.3f} (Li-err={r['Li_err']:+.3f})")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            xs = [r["x"] for r in rows]
            piE = [r["pi_exact"] for r in rows]
            tri = [r["pi_triangle"] for r in rows]
            Rv  = [r["R"] for r in rows]
            Lv  = [r["Li"] for r in rows]

            plt.figure()
            plt.plot(xs, piE, marker="o", label="pi(x) exact")
            plt.plot(xs, tri, marker="o", label="triangle approx")
            plt.plot(xs, Rv, marker="o", label="R(x)")
            plt.plot(xs, Lv, marker="o", label="Li(x)")
            plt.xscale("log")
            plt.xlabel("x (log scale)")
            plt.ylabel("counts")
            plt.title("Prime staircase vs triangle-gap approximation")
            plt.legend()
            plt.tight_layout()
            plt.savefig("triangle_vs_pi.png", dpi=160)
            print("Saved triangle_vs_pi.png")
        except Exception as e:
            print("Plotting failed (matplotlib may not be installed). CSV was written fine.")
            print("Error:", e)

if __name__ == "__main__":
    main()
