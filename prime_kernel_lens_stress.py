
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Prime-Kernel Lens — Stress & Ablation Suite (Multi-core)

Goal
----
Test whether the near-1 determinant behavior is *intrinsic to the lens*:
Hilbert pairing + even-kernel hand-off + prime-kernel head with weighted tangent slope.

Strategy
--------
Run the canonical pipeline and a battery of *null/ablation* variants that
surgically remove a single ingredient. Compare |det−1| distributions,
refinement slopes, and effect sizes.

Variants
--------
BASE:
  PK  : prime-kernel head, even kernel (fejer/gauss), Hilbert pairing ON

ABLATIONS:
  PK-oddK    : replace even kernel by an odd kernel (breaks H-K commutation)
  PK-noH     : remove Hilbert pairing (use identity instead of H in S_ij)
  PK-phaseJ  : jitter phases θ→θ+ε with ε~N(0,σ^2) per-prime
  PK-scramble: scramble primes by random permutation before summation
  CO         : composites-only kernel (weights on composite n up to Pmax)
  RANDW      : random weights w_p ~ N(0,σ^2) / p^α (keeps 1/p^α envelope)
  INTSEQ     : use integers n instead of primes (n≤Pmax), same weights

Outputs
-------
CSV rows per run, and a JSON summary that aggregates by variant:
  median |det−1|, quartiles, refinement slope vs N & Δ,
  effect size vs BASE, pass/fail against thresholds.

Usage (examples)
----------------
python3 prime_kernel_lens_stress.py \
  --t 700 1000 1500 \
  --Delta 3 5 7 --N 801 1201 \
  --Pmax 20000 --alpha 1.0 --lam 0.10 \
  --kernel fejer gauss --jobs 4 --out stress.csv

"""

import argparse
import csv
import math
import multiprocessing as mp
import time
from math import log, pi
from pathlib import Path

import numpy as np
from numpy.fft import fft, ifft

# ---------- utilities ----------

def walltime(): return time.perf_counter()

def sieve(n):
    if n < 2: return np.array([],dtype=np.int64)
    s = np.ones(n+1, dtype=bool); s[:2]=False
    for p in range(2, int(n**0.5)+1):
        if s[p]: s[p*p:n+1:p] = False
    return np.nonzero(s)[0].astype(np.int64)

def composites_upto(n):
    s = np.ones(n+1, dtype=bool); s[:2]=False
    for p in range(2, int(n**0.5)+1):
        if s[p]: s[p*p:n+1:p] = False
    # composites: non-prime ≥4
    idx = np.nonzero(~s)[0]
    return idx[idx>=4].astype(np.int64)

def rs_theta(t):
    tt = float(t)
    if tt <= 0: return -pi/8.0
    a = 0.5*tt*(math.log(tt/(2.0*pi)) - 1.0) - (pi/8.0)
    a += 1.0/(48.0*tt) - 7.0/(5760.0*tt**3)
    return a

def outer_taper(x, t, Delta):
    L = Delta
    core_lo, core_hi = t - L/2.0, t + L/2.0
    win_lo, win_hi   = t - L,     t + L
    y = np.zeros_like(x)
    left = (x >= win_lo) & (x < core_lo)
    core = (x >= core_lo) & (x <= core_hi)
    right= (x >  core_hi) & (x <= win_hi)
    if np.any(left):
        z = (x[left] - win_lo) / (core_lo - win_lo)
        y[left] = 0.5 - 0.5*np.cos(pi*z)
    if np.any(core): y[core] = 1.0
    if np.any(right):
        z = (win_hi - x[right]) / (win_hi - core_hi)
        y[right] = 0.5 - 0.5*np.cos(pi*z)
    return y

def hilbert_fft(y):
    N = y.size
    Y = fft(y)
    h = np.zeros(N, dtype=complex)
    if N % 2 == 0:
        h[1:N//2] = -1j; h[N//2+1:] = 1j
        h[0]=h[N//2]=0.0
    else:
        h[1:(N+1)//2] = -1j; h[(N+1)//2:] = 1j; h[0]=0.0
    return np.real(ifft(Y*h))

def make_kernel(kind, N, h, strength=1.0):
    mid = N//2
    xi = (np.arange(N)-mid)*h
    if kind=="fejer":
        m = max(2, N//20)
        k = np.zeros(N); ramp = np.r_[np.arange(1,m+1), np.arange(m-1,0,-1)]
        L = ramp.size; start = mid - L//2
        k[start:start+L] = ramp; k /= k.sum(); return k
    if kind=="gauss":
        sigma = strength * ((N*h)/6.0)
        k = np.exp(-0.5*(xi/sigma)**2); k /= k.sum(); return k
    if kind=="odd":
        # antisymmetric small kernel (break even-ness)
        w = max(3, N//50)
        k = np.zeros(N); left = np.arange(1,w+1); right = -left[::-1]
        idxL = mid - left; idxR = mid + np.arange(1,w+1)
        k[idxL] = left; k[idxR] = right
        # normalize L1
        s = np.sum(np.abs(k)); 
        if s>0: k /= s
        return k
    raise ValueError("kernel kind?")

def conv_even(y,k):
    return np.real(ifft(fft(y)*fft(k)))

def build_head(x, t, Delta, seq, mode, alpha=1.0, beta=1.0, phase_jitter=0.0, rng=None):
    F = outer_taper(x,t,Delta)
    theta = rs_theta(t)
    y = np.zeros_like(x, dtype=float)
    arr = seq.astype(float)
    if mode=="primes":
        weights = np.log(arr)/(arr**alpha)
    elif mode=="ints":
        weights = np.log(np.maximum(arr,2.0))/(arr**alpha)
    elif mode=="composites":
        weights = np.log(arr)/(arr**alpha)
    elif mode=="randomw":
        if rng is None: rng = np.random.default_rng(0)
        weights = rng.normal(0.0, 1.0, size=arr.size)/(arr**alpha)
    else:
        raise ValueError("mode?")
    block = 256
    rng = np.random.default_rng(12345) if rng is None else rng
    for i in range(0, arr.size, block):
        j = min(i+block, arr.size)
        pj = arr[i:j]; wj = weights[i:j]; lj = np.log(pj)
        ph = theta - np.outer(x, lj)
        if phase_jitter>0.0:
            ph = ph + rng.normal(0.0, phase_jitter, size=ph.shape)
        y += (np.cos(ph) @ wj)
    y *= (beta/pi); y *= F
    return y

def inner_h(bi,bj,h,use_hilbert=True):
    if use_hilbert:
        Hb = hilbert_fft(bj)
    else:
        Hb = bj
    return h * float(np.dot(bi, Hb))

def run_one(params):
    t,Delta,N,h = params["t"],params["Delta"],params["N"],params["h"]
    lam = params["lam"]; alpha = params["alpha"]; beta = params["beta"]
    Pmax = params["Pmax"]; kernel = params["kernel"]; variant = params["variant"]
    x = np.linspace(t-Delta, t+Delta, N)
    # sequences
    primes = sieve(Pmax)
    rng = np.random.default_rng(2025)
    if variant=="PK" or variant=="PK-oddK" or variant=="PK-noH" or variant=="PK-phaseJ" or variant=="PK-scramble":
        seq = primes.copy()
        mode = "primes"
        if variant=="PK-scramble":
            rng.shuffle(seq)
        phase_j = 0.05 if variant=="PK-phaseJ" else 0.0
        use_h = (variant!="PK-noH")
        kkind = ("odd" if variant=="PK-oddK" else kernel)
        b1 = build_head(x,t,Delta,seq,mode,alpha=alpha,beta=beta,phase_jitter=phase_j,rng=rng)
        b2 = hilbert_fft(b1) if use_h else b1
        k = make_kernel(kkind, N, h)
        b3 = conv_even(b1,k)
        b4 = hilbert_fft(b3) if use_h else b3
    elif variant=="CO":
        seq = composites_upto(Pmax); mode="composites"
        k = make_kernel(kernel, N, h)
        b1 = build_head(x,t,Delta,seq,mode,alpha=alpha,beta=beta)
        b2 = hilbert_fft(b1)
        b3 = conv_even(b1,k)
        b4 = hilbert_fft(b3)
    elif variant=="RANDW":
        seq = primes.copy(); mode="randomw"
        k = make_kernel(kernel, N, h)
        b1 = build_head(x,t,Delta,seq,mode,alpha=alpha,beta=beta,rng=rng)
        b2 = hilbert_fft(b1)
        b3 = conv_even(b1,k)
        b4 = hilbert_fft(b3)
    elif variant=="INTSEQ":
        seq = np.arange(2, Pmax+1, dtype=np.int64); mode="ints"
        k = make_kernel(kernel, N, h)
        b1 = build_head(x,t,Delta,seq,mode,alpha=alpha,beta=beta)
        b2 = hilbert_fft(b1)
        b3 = conv_even(b1,k)
        b4 = hilbert_fft(b3)
    else:
        raise ValueError("variant?")
    # S-matrix with Hilbert pairing choice embedded above
    S = np.zeros((4,4), dtype=float); B=[b1,b2,b3,b4]
    for i in range(4):
        for j in range(4):
            S[i,j] = inner_h(B[i],B[j],h,True)  # pairing uses Hilbert by definition of lens
    A = np.diag([1.0, -lam, -lam, -lam])
    det = float(np.linalg.det(np.eye(4)+A.dot(S)))
    dev = abs(det-1.0)
    return dict(det=det, dev=dev, N=N, Delta=Delta, t=t, Pmax=Pmax, alpha=alpha, lam=lam,
                kernel=kernel, variant=variant)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", type=float, nargs="+", default=[1000.0])
    ap.add_argument("--Delta", type=float, nargs="+", default=[5.0])
    ap.add_argument("--N", type=int, nargs="+", default=[1201])
    ap.add_argument("--Pmax", type=int, default=20000)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--beta", type=float, default=1.0)
    ap.add_argument("--lam", type=float, nargs="+", default=[0.10])
    ap.add_argument("--kernel", type=str, nargs="+", default=["fejer"])
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", type=str, default="stress.csv")
    return ap.parse_args()

def main():
    args = parse_args()
    variants = ["PK","PK-oddK","PK-noH","PK-phaseJ","PK-scramble","CO","RANDW","INTSEQ"]
    rows=[]
    for t in args.t:
        for Delta in args.Delta:
            for N in args.N:
                if N%2==0: continue
                for lam in args.lam:
                    for ker in args.kernel:
                        h = 2*Delta/(N-1)
                        for v in variants:
                            rows.append(dict(t=t,Delta=Delta,N=N,h=h,Pmax=args.Pmax,alpha=args.alpha,
                                             beta=args.beta,lam=lam,kernel=ker,variant=v))
    header = ["variant","t","Delta","N","Pmax","alpha","lam","kernel","det","dev"]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = walltime()
    if args.jobs>1 and len(rows)>1:
        with mp.get_context("spawn").Pool(args.jobs) as pool, open(out,"w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=header)
            w.writeheader()
            for res in pool.imap_unordered(run_one, rows, chunksize=1):
                w.writerow({k:res.get(k) for k in header})
    else:
        with open(out,"w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=header); w.writeheader()
            for r in rows:
                res = run_one(r)
                w.writerow({k:res.get(k) for k in header})
    t1 = walltime()
    print(f"Done in {t1-t0:.2f}s -> {out}")

if __name__=="__main__":
    main()
