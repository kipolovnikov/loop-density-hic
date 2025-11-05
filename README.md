# 🧬 Loop Density Inference from Hi-C Data

### Code and data accompanying the paper  
**“Universal contact statistics of looped polymers resolve cohesin density and stoichiometry in vivo”**  
*K. Polovnikov & D. Starkov (2025)*

---

## 📘 Overview

This repository provides the **data and analysis code** used in our study of short-range Hi-C contact probability curves \( P(s) \).  

It allows the user to **reproduce all figures** in the paper and to **infer the loop density (period \( T \))** and **effective fragment length (\( v_0^{\mathrm{eff}} \))** from any Hi-C dataset.

---

## 📂 Repository Structure
```
├── data/
│ ├── full_logder_x_.pickle # Mids (genomic distances)
│ ├── full_logder_y_.pickle # Log-derivatives (smoothed slopes of log P(s))
│ └── ... # All datasets used in the paper
│
├── notebooks/
│ ├── fig2A_logders.ipynb # Reproduces Fig. 2A (experimental log-derivatives)
│ ├── fig4A_RAD21_degron_mESC.ipynb.ipynb # Reproduces Fig. 4A (partial RAD21 degradation)
│ ├── fig4B_protocol_variation_ESC.ipynb # Reproduces Fig. 4B (protocol change)
│ ├── fig5A_fountain.ipynb # Reproduces Fig. 5A (fountain diagram)
│ ├── fig5A_fountain_aux.ipynb # Additional script with inferred parameters for Fig. 5A 
│ └── infer_params.ipynb # Example: fit model to Hi-C scaling or log-derivative
│
├── notebooks/src/
│ ├── infer_density.py # Helper functions for fitting T and v₀
│ |── utils.py # Helper functions to fit the data to the theoretical curves
| |── data_load.py # Helper functions to load the data
│
│
├── LICENSE
└── README.md
```

---

## ⚙️ Data Description

All data used in the paper are stored in the **`data/`** folder as Python pickle files.  
Each dataset contains:

- `full_logder_x_<dataset>.pickle` → mids (genomic distances, kb)  
- `full_logder_y_<dataset>.pickle` → log-derivative (smoothed derivative of log P(s))

These files correspond to the datasets analyzed in the manuscript (Abramo, Bonev, Rao, Zhang, Wutz, Schwarzer, etc.).

No preprocessing from raw `.cool` files is required — these intermediate products are **already included** for reproducibility.

---

## 🧩 Reproducing Paper Figures

All figure notebooks are provided in the `notebooks/` folder:

| Notebook | Description |
|-----------|-------------|
| **`fig2A_logders.ipynb`** | Plots all experimental log-derivatives used in Fig. 2A. |
| **`fig4A_RAD21_degron_mESC.ipynb`** | Reproduces Fig. 4A — RAD21 degradation (loop density reduction). |
| **`fig4B_protocol_variation_ESC.ipynb`** | Reproduces Fig. 4B (protocol variation). |
| **`fig5A_fountain.ipynb`** | Reproduces Fig. 5A (cross-species comparison). |
| **`infer_params.ipynb`** | Demonstrates parameter inference from new data. |

Each notebook loads pickled mids and log-derivatives directly from `data/` and plots publication-ready figures.

---

## 🔍 Inferring Loop Density from Your Own Data

You can infer \( T \) and \( v_0^{\mathrm{eff}} \) for any new dataset using the provided command-line script or notebook.

### 🧠 From a Jupyter notebook
```python
!python src/infer_density.py \
  --x data/full_logder_x_rao_GM12878_inSitu_DpnII.hg38.mapq_30.1000.mcool.pickle \
  --y data/full_logder_y_rao_GM12878_inSitu_DpnII.hg38.mapq_30.1000.mcool.pickle \
  --mode slope \
  --output-plot results/gm12878_fit.png


---

## ⚙️ Installation

### 1. Clone the repository
```bash
git clone https://github.com/kipolovnikov/loop-density-hic.git
cd loop-density-hic
```

### 2. Create the conda environment
```bash
conda env create -f environment.yml
conda activate loop-density-hic
```

### 3. (Optional) Install manually via pip
```bash
pip install cooler numpy scipy matplotlib pandas jupyter
```

---

## 🚀 Quick Start
Compute $P(s)$ and its log-derivative
```bash
python src/compute_ps_curve.py --cool data/example_coolers/GM12878.cool --binsize 1000
```

## Fit loop-density parameters
```bash
python src/fit_loop_density.py --input results/GM12878_ps.csv
```

## Or run the notebook interactively
```bash
jupyter notebook notebooks/01_compute_ps_and_derivative.ipynb
```

## 📊 Example Output

Each analysis produces:
- Contact probability $P(s)$
- Logarithmic derivative $\frac{d\log P}{d\log s}$
- Fitted parameters: loop period $T$ and efffective fragment length $v_0^{\mathrm{eff}}$
- Overlay plots comparing experimental and theoretical curves.

## 📈 Reproducing Figures

Each major figure in the paper can be regenerated via notebooks in notebooks/:
| Figure    | Notebook                          | Description                                             |
| --------- | --------------------------------- | ------------------------------------------------------- |
| Fig. 2A   | `03_reproduce_fig2.ipynb`         | Theory vs. data: characteristic “dip” in log-derivative |
| Fig. 4A   | `04_reproduce_fig4.ipynb`         | RAD21 degron perturbation: loop density reduction       |
| Fig. 5A–B | `05_stoichiometry_analysis.ipynb` | Cross-dataset inference and cohesin stoichiometry       |


## 🧠 Model Summary

At short genomic separations, the slope of the contact probability curve \(P(s)\) exhibits a conserved local minimum (“dip”) whose position scales with the geometric mean of loop period $T$ and $v_0^{\mathrm{eff}}$ and depth depends only on $v_0^{\mathrm{eff}}/T$:

$$
s_{\min} \sim \sqrt{T v_0^{\mathrm{eff}}}, \qquad
y_{\min} = f\left(\frac{v_0^{\mathrm{eff}}}{T}\right)
$$

This **two-parameter reduction** provides a universal fingerprint of cohesin loop density and experimental resolution across Hi-C protocols.

## 📚 Reference

If you use this repository, please cite:

Polovnikov, K., & Starkov, D. (2025).
Universal contact statistics of looped polymers resolve cohesin density and stoichiometry in vivo.
bioRxiv: 10.1101/2025.09.04.674214v1

## 🤝 Acknowledgments

Developed by Kirill Polovnikov and Dmitry Starkov.
We thank Job Dekker and Leonid Mirny for valuable feedback.
Supported by the Russian Science Foundation (Grant No. 25-13-00277) and the Alexander von Humboldt Foundation.
