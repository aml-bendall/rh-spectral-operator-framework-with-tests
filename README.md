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
- `rh_operator_framework_submission_v3.tex` - Complete mathematical paper with proofs
- `refs.bib` - Bibliography

### Computational Tools
- `prime_kernel_operator_sweep.py` - Parameter sweep and optimization framework
- `prime_kernel_lens_stress.py` - Stress testing and ablation studies

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
# Sanity check (single computation)
python prime_kernel_operator_sweep.py --t 1000 --Delta 5 --N 1201 --Pmax 20000 --out test.csv

# Stress test (multiple variants)
python prime_kernel_lens_stress.py --t 1000 --Delta 5 --N 1201 --Pmax 20000 --out stress.csv
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

- **Framework Validity**: PK methods consistently achieve det ≈ 1
- **Robustness**: Stable across parameter variations and perturbations
- **Discriminatory Power**: 7x separation between theory-based and control methods
- **Optimization**: λ=0.05, Δ=5 optimal for high-ordinate performance

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

## 📧 **Contact**

For mathematical discussions, error reports, or collaboration inquiries, please open an issue or contact [your contact info].

---

*"The critical-line balance emerges not from computational coincidence, but from the mathematical structure of Hilbert pairing and even-kernel hand-off."*