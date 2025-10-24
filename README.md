# 🧬 Transcriptomics-Based Drug Target Discovery in Diabetic Cardiomyopathy

## 📖 Project Overview
This project explores **differential gene expression** in Type 2 Diabetes–induced **diabetic cardiomyopathy (DCM)** using **iPSC-derived cardiomyocytes** under insulin resistance and hypoxia.  
The ultimate goal is to identify **potential drug targets** and understand molecular signatures driving disease progression.  

This work is aligned with my preparation for the **Erasmus Mundus S-DISCO Master's program** in Sustainable Drug Discovery.  

---

## 🖥️ Environment Setup
- **System**: Windows + VS Code  
- **Python**: via Conda environment (`t2d-env`)  
- **R**: installed separately with custom library path (`E:/Rlibs`)  
- **Jupyter**: integrated with IRkernel for running R + Python in the same workspace  

---

## 📂 Project Structure

```

project-root/
│
├── data/
│   ├── raw/                 # Original downloaded data (e.g., GEO series matrix, CSV)
│   ├── processed/           # Cleaned/processed outputs (aligned expression, design matrix)
│
├── results/
│   ├── plots/               # Figures/plots generated during analysis
│
├── notebooks/               # Jupyter notebooks for exploration & analysis
│
├── scripts/                 # Reusable Python/R scripts
│
└── README.md                # Project documentation

```

## 📊 Phase 3 Progress

### **Step 1: Data Access**
- Downloaded raw dataset **GSE288708** from GEO (series matrix + supplementary CSV).  
- Stored under `../data/raw/`.

### **Step 2: Metadata Parsing**
- Parsed **sample metadata** (`sample_metadata.csv`).  
- Extracted meaningful columns:  
  - `Sample` (GSM IDs)  
  - `Treatment` (Control / Insulin Resistant)  
  - `Oxygen` (Normoxia / Hypoxia)  
  - `Replicate`  
  - `Group` (combined condition)  
- Created a `ShortCode` (e.g., `CH1`, `IRN2`) for mapping samples.

### **Step 3: Expression Matrix Preparation**
- Loaded raw expression matrix.  
- Verified column names (short codes).  
- Checked consistency with metadata.  

### **Step 4: Code → GSM Mapping**
- Built dictionary to map **short codes → GSM IDs**.  
- Ensured no duplicates and no mismatches.

### **Step 5: Alignment**
- Renamed expression matrix columns to **GSM IDs**.  
- Reordered columns to exactly match metadata sample order.  
- Confirmed alignment with assertions.

### **Step 6: Export Outputs**
- Exported:  
  - **Aligned Expression Matrix** → `../data/processed/expression_matrix_aligned.csv`  
  - **Design Matrix** → `../data/processed/design_matrix.csv`  
- Sanity check done: dimensions matched, ready for downstream DE analysis.

---

## ✅ Current Status
- ✔️ Environment fully set up  
- ✔️ Dataset accessed and parsed  
- ✔️ Expression + metadata successfully aligned  
- ⏭️ Next step: **Data cleaning / normalization** before DESeq2  

---

## 🚀 How to Use
1. Clone repository / open in VS Code.  
2. Ensure `t2d-env` conda environment is activated.  
3. Run preprocessing script(s) or Jupyter notebooks in order.  
4. Outputs will be generated under `../data/processed/`.  

---

## 📝 Notes
- This README reflects progress **up to Step 6 of Phase 3**.  
- Further updates (data cleaning, DE analysis, visualization) will be documented in future commits.  
