#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick small-x benchmark to sanity check normalization and partial summation."""
import subprocess, sys, csv

def run(cmd):
    print("$", " ".join(cmd)); sys.stdout.flush()
    return subprocess.check_output(cmd, text=True)

def read_csv(path):
    rows=[]
    with open(path, newline="") as f:
        r=csv.DictReader(f)
        for d in r:
            rows.append({k: float(d[k]) if k!="x" else int(float(d[k])) for k in d})
    return rows

def main():
    out = "RL2b_corrected_psum_demo.csv"
    cmd = ["python3", "prime_kernel_pi_lens_logdomain_psum.py",
           "--x","1000","10000","100000","1000000",
           "--Nmax","2000000",
           "--kernel","fejer","--Delta_u","0.05","--lam","0.05",
           "--out", out]
    try:
        run(cmd)
        rows = read_csv(out)
        print("\n x      |  pi_hat   |  pi(x)  |  err   ||   R(x)   | R-err   |   Li(x)  | Li-err ")
        for d in rows:
            x   = d["x"]; piH = d["pi_hat"]; piE = d["pi_exact"]
            Rx  = d["R"]; lix = d["Li"]
            print(f"{x:<7}| {piH:8.3f}| {piE:7.0f}| {piH-piE:+7.3f} || {Rx:8.3f}| {Rx-piE:+7.3f} | {lix:8.3f}| {lix-piE:+7.3f}")
    except Exception as e:
        print("Note: This demo may exceed the time budget here. Run locally if it times out.")
        print("Error:", e)

if __name__ == "__main__":
    main()
