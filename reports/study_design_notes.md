# 🧠 Study Design Notes — Phase 3, Step 7
**Date:** 2025-08-20  
**Project:** Transcriptomics-Based Drug Target Discovery in Diabetic Cardiomyopathy (T2D)

---

## 1. Biological Context
Diabetic cardiomyopathy (DCM) is a major complication of type 2 diabetes (T2D), characterized by cardiac dysfunction in the absence of coronary artery disease or hypertension. Two key drivers of DCM pathology are:  

- **Insulin resistance (IR):** impairs glucose uptake, alters cardiac metabolism, and causes toxic lipid accumulation in cardiomyocytes.  
- **Hypoxia (low O₂):** occurs in diabetic hearts due to microvascular damage, reducing oxygen supply to the myocardium.  

In this project, we study **iPSC-derived cardiomyocytes** under controlled perturbations of IR and hypoxia. This provides a disease-relevant **2 × 2 factorial design** that mimics critical stress conditions in the diabetic heart.

---

## 2. Main Biological Question
How do **insulin resistance** and **hypoxia**, separately and together, reprogram the cardiomyocyte transcriptome?  
- Which **molecular pathways** are perturbed?  
- Which **genes** may serve as potential **therapeutic targets** for preventing or treating DCM?  

---

## 3. Study Design Logic
We have **4 experimental groups**:  

1. **Control + Normoxia (CN):** Baseline cardiomyocytes under normal conditions.  
2. **Control + Hypoxia (CH):** Effect of oxygen deprivation without IR.  
3. **IR + Normoxia (IRN):** Effect of insulin resistance alone.  
4. **IR + Hypoxia (IRH):** Combined stress condition — the closest model to diabetic cardiomyopathy.  

This factorial design allows us to measure both **main effects** (individual contributions of IR and hypoxia) and **interaction effects** (whether IR modifies the hypoxia response).  

---

## 4. Planned Contrasts (for DE analysis)

- **IR effect (normoxia):** `IRN vs CN`  
  - Identifies genes dysregulated by **insulin resistance alone**.  
  - Biologically: reveals how impaired insulin signaling rewires metabolism and stress pathways in cardiomyocytes under otherwise healthy oxygen levels.  

- **Hypoxia effect (control):** `CH vs CN`  
  - Captures the **direct effect of low oxygen** in cardiomyocytes.  
  - Biologically: highlights hypoxia-induced stress responses, angiogenesis signals, and survival/death pathways independent of IR.  

- **Combined effect:** `IRH vs CN`  
  - Tests the **full disease-relevant condition** against baseline.  
  - Biologically: identifies the genes/pathways most dysregulated in “diabetic cardiomyopathy-like” conditions.  

- **Interaction effect:** `(IRH − IRN) vs (CH − CN)`  
  - A true **statistical interaction contrast**: tests whether the effect of hypoxia is **different in IR cells compared to controls**.  
  - Biologically: reveals **synergistic or sensitization effects**, e.g., does insulin resistance make cells more vulnerable to hypoxia?  

---

## 5. Relevance to Drug Discovery
Each contrast offers a different therapeutic insight:  

- **IR-specific genes:** potential targets for metabolic therapies improving insulin sensitivity.  
- **Hypoxia-specific genes:** targets for cardioprotective or pro-angiogenic interventions.  
- **Combined stress response:** high-priority genes that may represent true drivers of DCM pathogenesis.  
- **Interaction genes:** unique candidates where dual targeting of IR and hypoxia-related pathways may yield the strongest therapeutic effect.  

By integrating these results with pathway enrichment, we aim to **prioritize candidate targets and pathways** that can be further validated for **drug discovery in DCM**.

---

✅ **Next step:** Perform **data cleaning and QC** to ensure the expression matrix is ready for robust DESeq2 differential expression analysis.
