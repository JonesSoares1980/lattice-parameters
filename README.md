# Deterministic Combinatorial Calculation of Lattice Parameters Using Python

Official repository containing the Python implementation of a deterministic combinatorial framework designed to calculate unit-cell parameters directly from indexed powder X-ray diffraction data (d_hkl), bypassing complex global profile fitting routines such as Rietveld or Le Bail refinement.

## 🚀 Key Features & Architectural Stabilization

* **High Performance:** Resolves millions of algebraic matrix combinations in seconds utilizing high-performance NumPy arrays.
* **Algebraic Soundness:** Migrated from non-linear symbolic systems (SymPy) to a linearized reciprocal-metric tensor solver framework.
* **Massive Datasets Support:** Benchmark-tested with expanded arrays of up to 124 reflections and millions of simultaneous permutations.

## 📈 Reproducibility & Technical Highlights

* **Peak Selection Context:** Manual conditioning of hkl reflections (such as the SrCO3 dataset filtered down to 66 reflections to isolate true geometric metrics).
* **The 3x3 Filtering Update:** Adoption of full-row determinant evaluation (Bareiss elimination) to ensure structural rank checking, replacing the previous restrictive geometric pre-filter.
* **Statistical Dispersion (sigma):** Robust statistical validation with global averages converging smoothly to the ideal database constants and an absolute chemical accuracy of Z_exact = 11.9985.

## 🔬 Mathematical Foundation

The relationship between the interplanar spacing (d_hkl) and the reciprocal lattice constants is governed by the quadratic form:
1/d_hkl^2 = h^2 Q_{11} + k^2 Q_{22} + l^2 Q_{33} + 2hk Q_{12} + 2hl Q_{13} + 2kl Q_{23}

## 🤝 Acknowledgments & Citation

Special thanks to **Lloyd McDonald (Mac)** for his outstanding peer review, mathematical validation, and contributions to the reciprocal-metric tensor solver framework.

```bibtex
@article{soares2026deterministic,
  author  = {Soares, Jones and Soares, Jobson},
  title   = {Deterministic combinatorial calculation of lattice parameters using Python},
  journal = {Academia Materials Science},
  volume  = {3},
  number  = {3},
  year    = {2026},
  pages   = {1--11},
  doi     = {10.20935/AcadMatSci8544}
}
```
