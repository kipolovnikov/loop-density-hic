# 🧬 Loop Density Inference from Hi-C Data

### Code and data accompanying the paper  
**“A universal polymer signature in Hi-C resolves cohesin loop density and supports monomeric extrusion”**  
*K. Polovnikov & D. Starkov (2025)*

---

## 📘 Overview

This repository provides the **data and analysis code** used in our study of short-range Hi-C contact probability curves \( P(s) \).  

**Reproduces the main data–theory figures** from the paper and provides a computational framework to infer the **cohesin loop density** (inverse loop period $T^{-1}$) and the **effective fragment length** $v_0^{\text{eff}}$ from a precomputed Hi-C contact scaling curve $P(s)$.

---

## 📂 Repository Structure
```
├── data/
│ ├── full_logder_x_.pickle # Mids (genomic distances)
│ ├── full_logder_y_.pickle # Log-derivatives (smoothed slopes of log P(s))
│ ├── ... # Main datasets used in the paper
│ └── data_table_info.xlsx # The descriptive table of the datasets & the inferred parameters for each data
│
├── data/data_fig4
│ ├── full_logder_x_.pickle # Mids (genomic distances)
│ ├── full_logder_y_.pickle # Log-derivatives (smoothed slopes of log P(s))
│ └── ... # The datasets used for Fig. 4A and Fig. 4B
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
│ ├── utils.py # Helper functions to fit the data to the theoretical curves
| └── data_load.py # Helper functions to load the data
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

You can infer the **loop period** \(T\), **loop density** \(T^{-1}\), and **effective fragment length** \(v_0^{\mathrm{eff}}\) for any new dataset using either:

- the **interactive notebook** `notebooks/infer_params.ipynb`, or  
- the **command-line tool** `src/infer_density.py`.

Below is an example of using the command-line tool on two pickle files containing **mids** and **slopes** (log-derivatives) precomputed from a Hi-C `.cool` file:

```bash
python src/infer_density.py \
  --x data/full_logder_x_<dataset>.pickle \
  --y data/full_logder_y_<dataset>.pickle \
  --mode slope \
  --output-plot results/<dataset>_fit.png
```
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
