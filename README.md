# Deterministic Combinatorial Calculation of Lattice Parameters Using Python

Official repository for calculating unit-cell parameters from powder X-ray diffraction data using Python.

## 🚀 Key Features

* **Performance & Verification:** Utilizes vectorized NumPy arrays (`lattice_metric_trial.py`) and a reciprocal-metric tensor solver.
* **Universal CSV Input (`csv_input.py`):** Automatically detects delimiters, decimal commas, named columns, and validates Miller indices and $2\theta$ bounds.
* **Console Output:** Streams subset combinations with optional listing options.

## 📦 Repository Structure

* `updated_code.py`: Interactive GUI script with visual file selection and report export.
* `lattice_metric_trial.py`: Linear reciprocal-metric tensor solver kernel.
* `csv_input.py`: Named-column CSV parser and input validator.
* `verify_soares.py`: Standalone verification script.

```bash
pip install -r requirements.txt
```
*(Requires Python 3.12+, NumPy, SciPy, SymPy, and Tkinter).*

## 🧪 Reproducibility

* Run `python lattice_metric_trial.py --self-test` to validate reference benchmarks ($\text{SrSiO}_3$ and $\text{ZnO}$).

## 🤝 Acknowledgments

Special thanks to **Lloyd H. McDonald (Mac)** for peer review and contributions to the solver kernel and CSV parser.
