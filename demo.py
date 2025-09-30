#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RH Spectral Operator Framework - Quick Demo

This script demonstrates the core determinant computation and validates
the theoretical prediction det(I₄ + AS) ≈ 1 for the prime-kernel construction.

Run this first to verify the framework works on your system.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    print("🔬 RH Spectral Operator Framework - Demo")
    print("=" * 50)
    print()
    
    print("📋 Testing Requirements...")
    try:
        import numpy as np
        import pandas as pd
        print("✅ Core dependencies available")
    except ImportError as e:
        print(f"❌ Missing dependency: {e}")
        print("Please install: pip install -r requirements.txt")
        return
    
    print()
    print("🧮 Running Basic Computation...")
    print("Parameters: t=1000, Δ=5, N=1201, Pmax=20000")
    
    # Run a simple test case
    try:
        result = subprocess.run([
            sys.executable, "prime_kernel_operator_sweep.py",
            "--t", "1000",
            "--Delta", "5", 
            "--N", "1201",
            "--Pmax", "20000",
            "--alpha", "1.0",
            "--lam", "0.10",
            "--kernel", "fejer",
            "--jobs", "1",
            "--out", "demo_result.csv"
        ], capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print("✅ Computation completed successfully")
            
            # Read and display results
            if Path("demo_result.csv").exists():
                df = pd.read_csv("demo_result.csv")
                row = df.iloc[0]
                
                print()
                print("📊 Results:")
                print(f"   Determinant: {row['det']:.6f}")
                print(f"   Deviation from 1: {row['dev']:.6f}")
                print(f"   Computation time: {row['t_total']:.3f}s")
                print(f"   Primes used: {row['nprimes']}")
                
                print()
                if row['dev'] < 0.01:
                    print("🎯 SUCCESS: Determinant very close to 1 (as expected)")
                    print("   The theoretical prediction det(I₄ + AS) ≈ 1 is confirmed")
                else:
                    print("⚠️  WARNING: Deviation larger than expected")
                    print("   This might indicate a computational issue")
                
                # Clean up
                Path("demo_result.csv").unlink(missing_ok=True)
                
            else:
                print("❌ Output file not created")
        else:
            print(f"❌ Computation failed: {result.stderr}")
            
    except subprocess.TimeoutExpired:
        print("❌ Computation timed out (>60s)")
    except Exception as e:
        print(f"❌ Error running computation: {e}")
        return
    
    print()
    print("🧪 Running Stress Test Sample...")
    
    try:
        result = subprocess.run([
            sys.executable, "prime_kernel_lens_stress.py",
            "--t", "1000",
            "--Delta", "5",
            "--N", "1201", 
            "--Pmax", "20000",
            "--alpha", "1.0",
            "--lam", "0.10",
            "--kernel", "fejer",
            "--jobs", "1",
            "--out", "demo_stress.csv"
        ], capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            print("✅ Stress test completed")
            
            if Path("demo_stress.csv").exists():
                df = pd.read_csv("demo_stress.csv")
                
                print()
                print("📈 Stress Test Results:")
                print("   Variant Performance:")
                
                for variant in ['PK', 'PK-noH', 'RANDW']:
                    if variant in df['variant'].values:
                        dev = df[df['variant'] == variant]['dev'].iloc[0]
                        print(f"     {variant:8s}: dev = {dev:.6f}")
                
                pk_dev = df[df['variant'] == 'PK']['dev'].iloc[0]
                randw_dev = df[df['variant'] == 'RANDW']['dev'].iloc[0]
                
                if pk_dev < randw_dev:
                    ratio = randw_dev / pk_dev
                    print(f"   🎯 PK method {ratio:.1f}x better than random (as expected)")
                else:
                    print("   ⚠️  Unexpected: Random method performing better")
                
                # Clean up
                Path("demo_stress.csv").unlink(missing_ok=True)
        else:
            print(f"❌ Stress test failed: {result.stderr}")
            
    except Exception as e:
        print(f"❌ Error running stress test: {e}")
    
    print()
    print("🎉 Demo Complete!")
    print()
    print("📚 Next Steps:")
    print("   1. Read the full paper: rh_operator_framework_submission_v3.tex")
    print("   2. Explore parameter sweeps with prime_kernel_operator_sweep.py")
    print("   3. Run comprehensive validation with prime_kernel_lens_stress.py")
    print("   4. Examine the validation data in *.csv files")
    print()
    print("❓ Questions or Issues:")
    print("   - Check README.md for detailed documentation")
    print("   - Review the theoretical framework in the LaTeX paper")
    print("   - Open GitHub issues for bugs or mathematical questions")

if __name__ == "__main__":
    main()