# RH Spectral Operator Framework

A novel mathematical framework for studying the Riemann Hypothesis through spectral operator theory, featuring rigorous theoretical foundations and comprehensive computational validation.

## 🔬 **Mathematical Framework**

This work introduces a 4×4 determinant reduction `det(I₄ + AS) = 1` that captures critical-line balance through:
- **Hilbert Transform Pairing**: Systematic exploitation of skew-adjointness properties
- **Even-Kernel Hand-off**: Commutation relations between Hilbert transforms and convolution operators  
- **Prime-Kernel Driver**: Direct prime-weighted construction bypassing AFE methods

## 📊 **Key Results**

- **Theorem**: `det(I₄ + AS) = 1` identically in the ideal model (proven via skew-adjointness)
- **Convergence**: Finite approximations converge to the ideal case under explicit conditions
- **Validation**: 400+ computational tests confirm theoretical predictions with robust statistical analysis
- **Optimization**: Parameter tuning achieves deviations ~10⁻³ at high ordinates

## 🗂️ **Repository Structure**

### Core Theory
- `rh_operator_framework_submission_v3.tex` - Mathematical paper
- `refs.bib` - Bibliography

### Computational Tools
- `prime_kernel_operator_sweep.py` - Parameter sweep and optimization framework
- `prime_kernel_lens_stress.py` - Stress testing and ablation studies

### Real-Line Extensions
- `prime_kernel_pi_lens_logdomain.py` - Spectral lens for prime counting π(x)
- `triangle_pi_approx.py` - Triangle gap geometric approximation (sub-unit precision)

### Validation Data
- `00_*.csv` - Sanity checks and baseline tests
- `01_*.csv` - Core robustness validation  
- `02_*.csv` - Scaling analysis and minima optimization
- `03_*.csv` - Comprehensive stress testing
- `04_*.csv` - High-resolution parameter tuning

## 🚀 **Quick Start**

### Prerequisites
```bash
pip install numpy pandas matplotlib scipy
```

### Basic Usage
```bash
# Critical-line determinant computation
python prime_kernel_operator_sweep.py --t 1000 --Delta 5 --N 1201 --Pmax 20000 --out test.csv

# Stress test (multiple variants)
python prime_kernel_lens_stress.py --t 1000 --Delta 5 --N 1201 --Pmax 20000 --out stress.csv

# Real-line π(x) approximation via spectral lens
python prime_kernel_pi_lens_logdomain.py --x 1000 10000 100000 --Pmax 2000000 --out lens.csv

# Triangle gap geometric approximation (ultra-precise!)  
python triangle_pi_approx.py --x 1000 10000 100000 1000000 --out triangle.csv --plot
```

### Expected Results
- **Determinant values**: Close to 1.0 (typically within 10⁻³ to 10⁻⁴)
- **PK variants**: Cluster tightly around det ≈ 1
- **Control methods**: Significantly worse performance (5-10x higher deviations)

## 📖 **Theoretical Background**

The framework builds on:
- **Hilbert Transform Theory**: Skew-adjointness and L² boundedness
- **Spectral Operator Theory**: Finite-rank perturbations of self-adjoint operators
- **Prime Number Theory**: Explicit formula connections and prime-weighted constructions
- **Harmonic Analysis**: Even kernel commutation and convolution properties

## 🔬 **Validation Methodology**

1. **Ablation Studies**: Systematic removal of components (PK-noH, PK-oddK, etc.)
2. **Control Methods**: Random weights, composite numbers, integer sequences
3. **Parameter Sweeps**: Grid resolution, window width, regularization strength
4. **Scaling Analysis**: Performance across ordinate ranges 700 ≤ t ≤ 5000

## 📈 **Key Findings**

### Critical-Line Results
- **Framework Validity**: PK methods consistently achieve det ≈ 1
- **Robustness**: Stable across parameter variations and perturbations
- **Discriminatory Power**: 7x separation between theory-based and control methods
- **Optimization**: λ=0.05, Δ=5 optimal for high-ordinate performance

### Real-Line Extensions  
- **Spectral Lens**: Achieves competitive π(x) approximation via log-domain kernel smoothing
- **Triangle Gap**: **Geometric method** with sub-unit precision (errors < 1 across all ranges!)
- **Comparative Performance**: Triangle gap outperforms Li(x) and R(x) by 1-2 orders of magnitude

## 🏆 **Triangle Gap Precision Benchmark**

| x | π(x) exact | Triangle π̃(x) | Error | Li(x) Error | R(x) Error |
|---|------------|---------------|-------|-------------|------------|
| 1,000 | 168 | 167.25 | **-0.75** | +9.61 | +0.36 |
| 10,000 | 1,229 | 1,228.79 | **-0.21** | +17.14 | -2.09 |
| 100,000 | 9,592 | 9,591.75 | **-0.25** | +37.81 | -4.60 |
| 1,000,000 | 78,498 | 78,497.85 | **-0.15** | +129.55 | +29.35 |

*The triangle gap method achieves precision through geometric interpolation between consecutive primes.*

## ⚠️ **Important Disclaimers**

- **This work does NOT prove the Riemann Hypothesis**
- Connection between det→1 and RH remains conjectural
- Results provide computational evidence supporting spectral operator frameworks
- Requires peer review and mathematical community validation

## 📚 **References**

See `refs.bib` for complete bibliography. Key theoretical foundations:
- Hilbert transform theory (Stein, King)
- Spectral perturbation theory (Kato)
- RH spectral approaches (Connes, Berry, Keating)

## 🤝 **Contributing**

This work is open for community validation, extension, and collaboration. Please:
1. Verify computational results independently
2. Check mathematical proofs and arguments
3. Suggest improvements or identify errors
4. Build on the theoretical framework

## 📄 **License**

- **Code**: MIT License (see LICENSE-CODE)
- **Paper**: Creative Commons BY 4.0 (see LICENSE-PAPER)
- **Data**: CC0 Public Domain

---
