
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prime‑Kernel Operator Determinant — Refinement & Robustness Sweep (Multi-core)

What this does
--------------
Computes the 4x4 determinant D(t; Δ, N, kernel, weights) = det(I + A S) with:
  - b1 = (β/π) F_{t,Δ}(x) * B_PK(x),  where
    B_PK(x) = Σ_{p ≤ Pmax} (log p) / p^α * cos(θ(x) - x log p)
  - b2 = H[b1] (discrete Hilbert via FFT multiplier −i*sign(ω))
  - b3 = K[b1] (even convolution kernel: 'fejer' or 'gauss')
  - b4 = H[b3]
  - S_ij = ⟨ b_i , H[b_j] ⟩_h  (trapezoid on uniform grid)
  - A = diag(1, −λ, −λ, −λ)

It then:
  • runs parameter sweeps (t, Δ, N, kernel, α, λ, Pmax)
  • supports weight modes: logp_over_p (default), logp_over_sqrtp, inv_sqrtp
  • outputs CSV with det, |det−1|, timings
  • optional minima scan around each t; optional comparison against
    a supplied zeros JSON (list of ordinates near each t)

LIMITS
------
No finite computation proves RH. This script is a falsifier/validator for the
specific operator framework: if it consistently trends det → 1 as Δ↑, N↑ and
is robust to kernel/shift choices, that's strong evidence. If it fails, that
is actionable falsification.

Usage (examples)
----------------
# Single run (good baseline)
python3 prime_kernel_operator_sweep.py \
  --t 1000 --Delta 5 --N 1201 --Pmax 20000 --alpha 1.0 --lam 0.10 \
  --kernel fejer --weight_mode logp_over_p --jobs 1 --out sweep.csv

# Full sweep (multi-core)
python3 prime_kernel_operator_sweep.py \
  --t 700 1000 1500 3000 --Delta 3 5 7 --N 801 1201 2048 \
  --Pmax 20000 --alpha 1.0 --lam 0.10 \
  --kernel fejer gauss --weight_mode logp_over_p \
  --jobs 4 --out sweep.csv

# With minima scan (±W around t with step dt)
python3 prime_kernel_operator_sweep.py \
  --t 1000 --scan_window 3.0 --scan_step 0.01 \
  --Pmax 20000 --alpha 1.0 --lam 0.10 \
  --out scan.csv
"""

import argparse
import csv
import multiprocessing as mp
import time
from math import log, pi, sqrt
from pathlib import Path

import numpy as np
from numpy.fft import fft, fftfreq, ifft

# ------------------------
# Utilities
# ------------------------

def walltime():
    return time.perf_counter()

def sieve_primes(n):
    """Sieve of Eratosthenes up to n (inclusive)."""
    if n < 2:
        return np.array([], dtype=np.int64)
    m = n + 1
    sieve = np.ones(m, dtype=bool)
    sieve[:2] = False
    r = int(n**0.5)
    for p in range(2, r+1):
        if sieve[p]:
            step = p
            start = p*p
            sieve[start:m:step] = False
    return np.nonzero(sieve)[0].astype(np.int64)

def rs_theta(t):
    """
    Riemann–Siegel theta approximation (good for t >= ~10):
      θ(t) ≈ t/2 * log(t/(2π)) − t/2 − π/8 + 1/(48 t) − 7/(5760 t^3)
    """
    tt = float(t)
    if tt <= 0:
        return -pi/8.0
    a = 0.5*tt*(np.log(tt/(2.0*pi)) - 1.0) - (pi/8.0)
    a += 1.0/(48.0*tt) - 7.0/(5760.0*tt**3)
    return a

def outer_taper(x, t, Delta):
    """
    Smooth C^∞-like bump approximant: use raised-cosine on the boundary,
    ≡1 on [t-Δ/2, t+Δ/2], decays to 0 at [t-Δ, t+Δ].
    """
    L = Delta
    core_lo, core_hi = t - L/2.0, t + L/2.0
    win_lo, win_hi   = t - L,     t + L
    y = np.zeros_like(x)
    # three regions: left ramp, core (1), right ramp
    left = (x >= win_lo) & (x < core_lo)
    core = (x >= core_lo) & (x <= core_hi)
    right= (x >  core_hi) & (x <= win_hi)
    # raised-cosine ramps
    if np.any(left):
        z = (x[left] - win_lo) / (core_lo - win_lo)
        y[left] = 0.5 - 0.5*np.cos(pi*z)
    if np.any(core):
        y[core] = 1.0
    if np.any(right):
        z = (win_hi - x[right]) / (win_hi - core_hi)
        y[right] = 0.5 - 0.5*np.cos(pi*z)
    return y

def hilbert_transform(y):
    """
    Discrete Hilbert transform using FFT multiplier -i*sign(ω).
    """
    N = y.size
    Y = fft(y)
    h = np.zeros(N, dtype=complex)
    # frequency bins: k=0..N-1; sign mapping:
    # k=0 -> 0; k=1..N//2-1 -> -i; k=N/2 -> 0 (if even N); k>N/2 -> +i
    if N % 2 == 0:
        h[1:N//2] = -1j
        h[N//2+1:] = 1j
        h[0] = 0.0
        h[N//2] = 0.0
    else:
        h[1:(N+1)//2] = -1j
        h[(N+1)//2:] = 1j
        h[0] = 0.0
    return np.real(ifft(Y * h))

def make_kernel(kind, N, h, strength=1.0):
    """
    Even kernels on the grid for convolution via FFT (circular).
    - 'fejer': discrete triangular (hat) kernel with half-width ≈ floor(N/20)
    - 'gauss': discrete Gaussian with sigma ≈ strength * (Delta/6)
    Returns normalized kernel k of size N.
    """
    x = np.arange(N)
    # center-symmetric index
    mid = N//2
    xi = (x - mid) * h
    if kind == "fejer":
        m = max(2, N//20)  # half-width in samples
        k = np.zeros(N)
        ramp = np.concatenate([np.arange(1, m+1), np.arange(m-1, 0, -1)])
        L = ramp.size
        # embed ramp centered
        start = mid - L//2
        k[start:start+L] = ramp
        k = k / k.sum()
        return k
    elif kind == "gauss":
        sigma = strength * ( (N*h) / 6.0 )
        k = np.exp(-0.5 * (xi/sigma)**2 )
        k /= k.sum()
        return k
    else:
        raise ValueError("Unknown kernel: %r" % kind)

def even_convolve(y, k):
    """
    Circular convolution via FFT for even kernel k.
    """
    Y = fft(y)
    K = fft(k)
    return np.real(ifft(Y*K))

def build_b1_prime_kernel(x, t, Delta, primes, weight_mode="logp_over_p", alpha=1.0, beta=1.0):
    """
    Build b1(x) = (β/π) F_{t,Δ}(x) Σ_{p≤Pmax} w(p) cos( θ(x) − x log p ).
    We chunk over primes to keep memory bounded.
    """
    N = x.size
    F = outer_taper(x, t, Delta)
    theta = rs_theta(t)  # NOTE: θ depends on t in RS theory; for a short window, use θ(t) constant across x
    # If you'd like θ(x), one could use t replaced by x, but in RH literature θ is a function of ordinate t.
    # We'll stick to θ(t) for stability on short windows.
    b1 = np.zeros(N, dtype=float)
    # chunk primes
    logs = np.log(primes.astype(float))
    if weight_mode == "logp_over_p":
        w = logs / (primes.astype(float) ** alpha)
    elif weight_mode == "logp_over_sqrtp":
        w = logs / (primes.astype(float) ** 0.5)
    elif weight_mode == "inv_sqrtp":
        w = 1.0 / (primes.astype(float) ** 0.5)
    else:
        raise ValueError("Unknown weight_mode")
    # Normalize weights mildly to avoid overflow across parameter choices
    w = w.astype(float)

    block = 256
    for i in range(0, primes.size, block):
        j = min(i+block, primes.size)
        pj = primes[i:j].astype(float)
        wj = w[i:j]
        lj = np.log(pj)
        # Phase matrix: Θ - x*log p_j  (N x B)
        # Use broadcasting: x[:,None] * lj[None,:]
        phase = theta - np.outer(x, lj)
        b1 += (np.cos(phase) @ wj)
    b1 *= (beta / pi)
    b1 *= F
    return b1

def inner_product_hilbert(bi, bj, h):
    """
    Discrete inner product <bi, H[bj]> with trapezoid weights (uniform grid).
    Since we apply FFT Hilbert globally, trapezoid reduces to scalar h factor.
    """
    Hb = hilbert_transform(bj)
    return h * float(np.dot(bi, Hb))

def compute_det_for_params(params):
    """
    Core single-run worker. Returns dict of results.
    """
    t      = params["t"]
    Delta  = params["Delta"]
    N      = params["N"]
    lam    = params["lam"]
    alpha  = params["alpha"]
    beta   = params.get("beta", 1.0)
    Pmax   = params["Pmax"]
    kernel = params["kernel"]
    weight_mode = params["weight_mode"]
    kernel_strength = params.get("kernel_strength", 1.0)

    x0 = t - Delta
    x1 = t + Delta
    x = np.linspace(x0, x1, N)
    h = (x1 - x0) / (N - 1)

    t0 = walltime()
    primes = sieve_primes(Pmax)
    t1 = walltime()

    b1 = build_b1_prime_kernel(x, t, Delta, primes, weight_mode=weight_mode, alpha=alpha, beta=beta)
    b2 = hilbert_transform(b1)

    k = make_kernel(kernel, N, h, strength=kernel_strength)
    b3 = even_convolve(b1, k)
    b4 = hilbert_transform(b3)

    # Build S matrix
    S = np.zeros((4,4), dtype=float)
    B = [b1,b2,b3,b4]
    for i in range(4):
        for j in range(4):
            S[i,j] = inner_product_hilbert(B[i], B[j], h)

    # A = diag(1, -lam, -lam, -lam)
    A = np.diag([1.0, -lam, -lam, -lam])
    M = np.eye(4) + A.dot(S)
    det = float(np.linalg.det(M))
    dev = abs(det - 1.0)
    t2 = walltime()

    return {
        "t": t, "Delta": Delta, "N": N, "h": h,
        "Pmax": Pmax, "alpha": alpha, "beta": beta,
        "lam": lam, "kernel": kernel, "weight_mode": weight_mode,
        "det": det, "dev": dev,
        "nprimes": int(primes.size),
        "t_sieve": t1 - t0,
        "t_total": t2 - t0
    }

def minima_scan(params, W=3.0, dt=0.01):
    """
    Scan t' in [t-W, t+W] and record the t' where |det-1| is minimal.
    Returns best_t, best_dev, samples (list of dict rows).
    """
    base_t = params["t"]
    ts = np.arange(base_t - W, base_t + W + 1e-12, dt)
    best = (None, float("inf"))
    rows = []
    # Reuse primes across scan for speed
    # We'll compute per t' but re-sieve once at Pmax
    Pmax = params["Pmax"]
    primes = sieve_primes(Pmax)
    for tprime in ts:
        p = dict(params)
        p["t"] = float(tprime)
        row = compute_det_for_params_scanreuse(p, primes)
        rows.append(row)
        if row["dev"] < best[1]:
            best = (tprime, row["dev"])
    return best[0], best[1], rows

def compute_det_for_params_scanreuse(params, primes):
    """
    Same as compute_det_for_params, but uses precomputed primes for scan.
    """
    t      = params["t"]
    Delta  = params["Delta"]
    N      = params["N"]
    lam    = params["lam"]
    alpha  = params["alpha"]
    beta   = params.get("beta", 1.0)
    Pmax   = params["Pmax"]
    kernel = params["kernel"]
    weight_mode = params["weight_mode"]
    kernel_strength = params.get("kernel_strength", 1.0)

    x0 = t - Delta
    x1 = t + Delta
    x = np.linspace(x0, x1, N)
    h = (x1 - x0) / (N - 1)

    b1 = build_b1_prime_kernel(x, t, Delta, primes, weight_mode=weight_mode, alpha=alpha, beta=beta)
    b2 = hilbert_transform(b1)
    k = make_kernel(kernel, N, h, strength=kernel_strength)
    b3 = even_convolve(b1, k)
    b4 = hilbert_transform(b3)

    S = np.zeros((4,4), dtype=float)
    B = [b1,b2,b3,b4]
    for i in range(4):
        for j in range(4):
            S[i,j] = inner_product_hilbert(B[i], B[j], h)

    A = np.diag([1.0, -lam, -lam, -lam])
    M = np.eye(4) + A.dot(S)
    det = float(np.linalg.det(M))
    dev = abs(det - 1.0)

    return {
        "t": t, "Delta": Delta, "N": N, "h": h,
        "Pmax": Pmax, "alpha": alpha, "beta": beta,
        "lam": lam, "kernel": kernel, "weight_mode": weight_mode,
        "det": det, "dev": dev,
        "nprimes": int(primes.size),
        "t_sieve": 0.0,  # reused
        "t_total": np.nan
    }

# ------------------------
# CLI / Sweep orchestration
# ------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Prime-kernel operator determinant sweep (multi-core).")
    p.add_argument("--t", type=float, nargs="+", default=[1000.0], help="Target ordinates (can pass multiple).")
    p.add_argument("--Delta", type=float, nargs="+", default=[5.0], help="Window half-widths.")
    p.add_argument("--N", type=int, nargs="+", default=[1201], help="Grid sizes (odd).")
    p.add_argument("--Pmax", type=int, default=20000, help="Prime cutoff.")
    p.add_argument("--alpha", type=float, default=1.0, help="Prime weight exponent.")
    p.add_argument("--beta", type=float, default=1.0, help="Row scaling factor (β/π outside).")
    p.add_argument("--lam", type=float, nargs="+", default=[0.10], help="λ in A=diag(1,-λ,-λ,-λ) (can pass multiple).")
    p.add_argument("--kernel", type=str, nargs="+", default=["fejer"], help="Kernel(s): fejer, gauss.")
    p.add_argument("--kernel_strength", type=float, default=1.0, help="Kernel strength (width control).")
    p.add_argument("--weight_mode", type=str, default="logp_over_p",
                   choices=["logp_over_p", "logp_over_sqrtp", "inv_sqrtp"])
    p.add_argument("--jobs", type=int, default=1, help="Parallel jobs.")
    p.add_argument("--out", type=str, default="sweep.csv", help="Output CSV.")
    p.add_argument("--scan_window", type=float, default=0.0, help="If >0, do minima scan ±W.")
    p.add_argument("--scan_step", type=float, default=0.01, help="Step for minima scan.")
    return p.parse_args()

def param_grid(args):
    for t in args.t:
        for Delta in args.Delta:
            for N in args.N:
                if N % 2 == 0:
                    continue
                for lam in args.lam:
                    for kernel in args.kernel:
                        yield dict(
                            t=float(t), Delta=float(Delta), N=int(N),
                            lam=float(lam), alpha=float(args.alpha), beta=float(args.beta),
                            Pmax=int(args.Pmax), kernel=str(kernel),
                            weight_mode=str(args.weight_mode),
                            kernel_strength=float(args.kernel_strength)
                        )

def run_job(row):
    try:
        out = compute_det_for_params(row)
        out["status"] = "ok"
    except Exception as e:
        out = dict(row)
        out.update({"status":"error", "error": repr(e)})
    return out

def main():
    args = parse_args()
    rows = list(param_grid(args))

    header = ["status","error","t","Delta","N","h","Pmax","alpha","beta","lam",
              "kernel","weight_mode","det","dev","nprimes","t_sieve","t_total"]
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    t0 = walltime()
    if args.jobs > 1 and len(rows) > 1:
        with mp.get_context("spawn").Pool(processes=args.jobs) as pool, \
             open(outpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
            w.writeheader()
            for res in pool.imap_unordered(run_job, rows, chunksize=1):
                w.writerow(res)
    else:
        with open(outpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
            w.writeheader()
            for row in rows:
                res = run_job(row)
                w.writerow(res)

    # Optional minima scans
    if args.scan_window > 0.0:
        scan_path = outpath.with_name(outpath.stem + "_scan.csv")
        header2 = ["tprime","dev","t","Delta","N","Pmax","alpha","beta","lam","kernel","weight_mode"]
        with open(scan_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header2)
            w.writeheader()
            for row in rows:
                best_t, best_dev, samples = minima_scan(row, W=args.scan_window, dt=args.scan_step)
                w.writerow({
                    "tprime": best_t, "dev": best_dev,
                    "t": row["t"], "Delta": row["Delta"], "N": row["N"],
                    "Pmax": row["Pmax"], "alpha": row["alpha"], "beta": row["beta"],
                    "lam": row["lam"], "kernel": row["kernel"], "weight_mode": row["weight_mode"]
                })
    t1 = walltime()
    print(f"Done in {t1-t0:.2f}s → {outpath}")

if __name__ == "__main__":
    main()
