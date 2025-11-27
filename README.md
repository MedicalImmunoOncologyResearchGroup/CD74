# 🧬 CD74

[![DOI](https://zenodo.org/badge/DOI/your-doi-here.svg)](https://doi.org/your-doi-here)

> **Paper Citation:** Maranou et al. (2026). , *Journal Name*. DOI: [link]

## 📋 Overview

This repository contains scripts for analyzing single-cell RNA-seq data presented in the manuscript. 

## 🗂️ Repository Structure

```
.
├── data/           # Data files
├── R/              # Analysis scripts in R
└── Jupyter/        # Analysis scripts in Python as Jupyter notebbook files
   

## 🔧 Requirements

Python codes use mostly standard libraries such as

- pandas
- numpy
- seaborn
- matplotlib
- scikit-learn
- scanpy
- celltypist
- scanorama
- pyDESEq2

Jupyter notebooks were tested with Python version 3.10.12.

## 📊 Analysis Workflow

1. **R Codes** (`R/`)
   - A script for DEA with interactions, using DESeq2 

2. **Jupyter Notebooks** (`Jupyter/`)
   - Preprocessing, QC and integration of scRNA-seq data
   - Annotation of cell types 
   - Evaluation of gene expression in cell types
   - Pseudobulk DEA
   - Gene ontology analysis

## 🔍 Data Availability

Data are available at XXX.


## 📜 License

This project is licensed under MIT licence, see LICENSE file for details. 

## ✍️ Citation

If you use this code or data, please cite:

@article{Maranou2026,
  title={CD74 regulates antitumor immunity in melanoma by reprogramming dendritic cell immunogenicity and migration},
  author={E. Maranou, G. Park, P. Weinzettl, J. Alanko, O. Pulkkinen, S. E. Coupland,, M. Salmi, C.R. Figueiredo},
  journal={Journal Name},
  year={2025},
  doi={your-doi}
}
```


## 📫 Contact

* **Principal Investigator:** Carlos Rogerio de Figueiredo (rogerio.defigueiredo@utu.fi)
* **Lead Author:** Eleftheria Maranou (eleftheria.maranou@utu.fi)

---
*Repository maintained by Otto Pulkkinen*
# CD74
