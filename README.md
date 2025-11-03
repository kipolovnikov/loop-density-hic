# Loop Density Inference from Hi-C Data

### Code and data accompanying the paper:
**“Universal contact statistics of looped polymers resolve cohesin density and stoichiometry in vivo”**  
by *K. Polovnikov & D. Starkov (2025)*

---

## 🧩 Overview
This repository provides the full computational framework and example data used in our analysis of short-scale Hi-C contact statistics.  

It includes tools to:
- Compute **contact probability curves** $P(s)$ and their **log-derivatives** directly from `.cool` Hi-C files.  
- Fit **loop-density and protocol parameters** — the loop period $T$ and the effective fragment length $v_0^{\mathrm{eff}}$ — to experimental data.  
- Reproduce all **data-intensive figures** from the paper using interactive Jupyter notebooks.  

Our goal is to make the inference of **cohesin loop density** from Hi-C data transparent, reproducible, and accessible so that anyone can apply it to their own datasets.  
If you have questions or encounter issues, please don’t hesitate to contact us.

---

## 📁 Repository Structure
```
loop-density-hic/
│
├── data/ # Example Hi-C data and metadata
│ ├── example_coolers/
│ │ ├── GM12878.cool
│ │ ├── mESC_WT.cool
│ │ └── ...
│ └── metadata/
│ └── sample_info.csv
│
├── src/ # Core analysis scripts
│ ├── compute_ps_curve.py # Compute P(s) and log-derivative from cooler files
│ ├── fit_loop_density.py # Fit T and v0_eff parameters from P(s)
│ └── utils.py # Helper functions
│
├── notebooks/ # Jupyter notebooks reproducing paper figures
│ ├── 01_compute_ps_and_derivative.ipynb
│ ├── 02_fit_parameters.ipynb
│ ├── 03_reproduce_fig2.ipynb
│ ├── 04_reproduce_fig4.ipynb
│ ├── 05_stoichiometry_analysis.ipynb
│ └── ...
│
├── results/ # Output data and fitted parameters
│ ├── inferred_parameters.csv
│ └── figures/
│ ├── Fig2A_fit.png
│ ├── Fig4B_protocol_comparison.png
│ └── ...
│
├── environment.yml # Conda environment specification
├── LICENSE # License information (MIT for code, CC-BY 4.0 for data)
└── README.md
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
