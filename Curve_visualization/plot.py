# Update the user's code to include AUC values in the plot legends for both cell-type and cancer classification

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from sklearn.metrics import auc

plt.figure(figsize=(8, 6))
roc_files = sorted(glob.glob("rocdata_celltype_k*.txt"))

for filepath in roc_files:
    try:
        data = np.loadtxt(filepath, skiprows=1)
        if data.size == 0:
            continue
        fpr, tpr = data[:, 0], data[:, 1]
        auc_score = auc(fpr, tpr)
        fname = os.path.basename(filepath)
        k = fname.split("_k")[-1].replace(".txt", "")
        plt.plot(fpr, tpr, label=f"k={k} (AUC={auc_score:.2f})")
    except Exception as e:
        print(f"Error reading {filepath}: {e}")

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curves for Cell-Type Classification (k-mer sizes)")
plt.legend(title="k-mer", loc="lower right")
plt.grid(False)
plt.tight_layout()
plt.savefig("celltype_roc_comparison.png", dpi=300)

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from sklearn.metrics import auc

plt.figure(figsize=(8, 6))
roc_files = sorted(glob.glob("rocdata_cancer_k*.txt"))

for filepath in roc_files:
    try:
        data = np.loadtxt(filepath, skiprows=1)
        if data.size == 0:
            continue
        fpr, tpr = data[:, 0], data[:, 1]
        auc_score = auc(fpr, tpr)
        fname = os.path.basename(filepath)
        k = fname.split("_k")[-1].replace(".txt", "")
        plt.plot(fpr, tpr, label=f"k={k} (AUC={auc_score:.2f})")
    except Exception as e:
        print(f"Error reading {filepath}: {e}")

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curves for Cancer vs. Non-Cancer Classification")
plt.legend(title="k-mer", loc="lower right")
plt.grid(False)
plt.tight_layout()
plt.savefig("cancer_roc_comparison.png", dpi=300)

