
# DHS Sequence-Based Classification Using Interpretable Machine Learning

This project explores whether accessible DNA sequences (DNase I hypersensitive sites, or DHSs) can predict:
1. **Cell Type Identity** (e.g., GM12878, HepG2, K562, hESC)
2. **Cancer vs. Non-Cancer Phenotype**

We use k-mer frequency representations, GC content, and entropy as sequence features to train Random Forest classifiers. Performance is evaluated across a range of k-mer sizes (`k = 2,3,4,5,6,7,8`).

---

## 📁 Project Structure

```
.
├── cancer_classifier.py           # Binary classifier (Cancer vs. Non-Cancer)
├── celltype_classifier.py         # Multiclass classifier (Cell type)
├── feature_extraction.py          # Functions to extract sequence features
├── *.fa                           # Input DHS sequence files (FASTA)
├── *.txt                          # Output ROC/AUC data (saved during run)
├── *.png                          # Plots: Confusion matrix, ROC curves
└── README.md                      # This file
```

---

## 📥 Input FASTA Files

Each FASTA file contains DHS sequences from one cell line:
- `GM12878_DHS_sequences.fa`
- `hESCT0_DHS_sequences.fa`
- `HepG2_DHS_sequences.fa`
- `K562_DHS_sequences.fa`

Each sequence is 150–200 bp long and corresponds to an open chromatin region.

---

## ⚙️ Dependencies
---

## 🧪 Virtual Environment Setup (Recommended)

To isolate dependencies and avoid conflicts, use a Python virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```


To install required Python packages:
```bash
pip install -r requirements.txt
```

Or manually:
```bash
pip install numpy pandas matplotlib seaborn scikit-learn biopython
```

---

## 🚀 How to Run

### 1️⃣ Run Cell Type Classifier (4-class classification)
```bash
python celltype_classifier.py
```

- Trains Random Forest on k-mer features
- Evaluates precision, recall, F1-score
- Saves confusion matrices & ROC curves per `k`
- Exports `rocdata_celltype_k{k}.txt` and `auc_celltype.txt`

---

### 2️⃣ Run Cancer vs Non-Cancer Classifier (binary classification)
```bash
python cancer_classifier.py
```

- Maps GM12878, hESC → non-cancer (0), K562, HepG2 → cancer (1)
- Trains and evaluates across k-mer sizes
- Outputs confusion matrices, ROC plots
- Saves `rocdata_cancer_k{k}.txt` and `auc_cancer.txt`

---

## 📊 Outputs

| File                             | Description                                |
|----------------------------------|--------------------------------------------|
| `confusion_matrix_k{k}.png`      | Confusion matrix for cell type prediction  |
| `cancer_confusion_matrix_k{k}.png` | Confusion matrix for binary cancer test   |
| `celltype_all_roc.png`           | Overlay of ROC curves for cell types       |
| `cancer_all_roc.png`             | Overlay of ROC curves for cancer task      |
| `rocdata_*.txt`                  | Tab-delimited FPR/TPR data (for replotting)|
| `auc_*.txt`                      | Summary of AUC for each `k`                |

---

## 📘 Citation / Background

This project was developed as part of a functional genomics study to assess the predictive capacity of DNA sequence features from DHSs. Our goal was to build interpretable models with high biological relevance using only sequence data.

---

## 🧬 Contact

For any questions, contact **[Your Name]** at **[your-email@domain.com]**  
(Replace with your actual contact info)

---

## 🧠 License

MIT License — feel free to reuse and modify.
