#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runs RL-2b with a small grid and prints a compact comparison table vs Li(x), R(x).
"""
import csv, subprocess, sys, os, tempfile

def run(cmd):
    print("$", " ".join(cmd)); sys.stdout.flush()
    return subprocess.check_output(cmd, text=True)

def parse_csv(path):
    rows=[]
    with open(path, newline="") as f:
        r=csv.DictReader(f)
        for d in r:
            rows.append({k: float(d[k]) if k not in ("x",) else int(float(d[k])) for k in d})
    return rows

def main():
    xs = ["1000","10000","100000","1000000"]
    grids = [
        ("fejer", "0.05", "0.05"),
        ("fejer", "0.03", "0.05"),
        ("gauss", "0.05", "0.05"),
    ]
    for kernel, Delta_u, lam in grids:
        out = f"RL2b_{kernel}_Du{Delta_u}_lam{lam}.csv"
        cmd = ["python3","prime_kernel_pi_lens_logdomain.py","--x",*xs,
               "--Pmax","2000000","--kernel",kernel,"--Delta_u",Delta_u,"--lam",lam,"--out",out]
        run(cmd)
        rows = parse_csv(out)
        print(f"\n=== {kernel}, Δu={Delta_u}, λ={lam} ===")
        print("x   |  pi_hat  |  pi(x)  |  err   ||  R(x)  |  R-err  |  Li(x)  | Li-err")
        for d in rows:
            x   = int(d["x"])
            piH = d["pi_hat"]; piE = d["pi_exact"]
            Rx  = d["R"]; lix = d["Li"]
            print(f"{x:<6}| {piH:8.3f}| {piE:7.0f}| {piH-piE:+7.3f} || {Rx:8.3f}| {Rx-piE:+7.3f} | {lix:8.3f}| {lix-piE:+7.3f}")
        print()

if __name__=="__main__":
    main()