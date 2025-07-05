import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    classification_report, confusion_matrix,
    roc_auc_score, RocCurveDisplay, roc_curve
)
from sklearn.preprocessing import label_binarize
from sklearn.utils import resample
from feature_extraction import extract_features

# Input FASTA files
fasta_files = {
    "GM12878": "GM12878_DHS_sequences.fa",
    "HepG2": "HepG2_DHS_sequences.fa",
    "hESC": "hESCT0_DHS_sequences.fa",
    "K562": "K562_DHS_sequences.fa"
}

kmers_to_try = [5]
classes = ["GM12878", "HepG2", "hESC", "K562"]
all_roc = []

# Reset summary AUC file
with open("auc_celltype.txt", "w") as f:
    f.write("k-mer\tAUC\n")

for k in kmers_to_try:
    print(f"\n=== Processing k={k} ===")
    df = extract_features(fasta_files, k)

    grouped = [resample(g, replace=False, n_samples=6000, random_state=42)
               for _, g in df.groupby("cell_type")]
    df_balanced = pd.concat(grouped).sample(frac=1, random_state=42)

    X = df_balanced.drop(columns=["sequence_id", "cell_type"])
    y = df_balanced["cell_type"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=0.2, random_state=42)

    clf = RandomForestClassifier(n_estimators=100, class_weight="balanced", random_state=42)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=clf.classes_)
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap="Blues",
                xticklabels=clf.classes_, yticklabels=clf.classes_)
    plt.title(f"Confusion Matrix (k={k})")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.tight_layout()
    plt.savefig(f"confusion_matrix_k{k}.png")
    plt.close()

    # ROC-AUC
    y_test_bin = label_binarize(y_test, classes=classes)
    y_score = clf.predict_proba(X_test)
    auc_score = roc_auc_score(y_test_bin, y_score, average="macro", multi_class="ovr")
    print(f"Multiclass AUC (k={k}): {auc_score:.4f}")
    all_roc.append((k, y_test_bin, y_score))

    # Save ROC curve data
    fpr, tpr, _ = roc_curve(y_test_bin.ravel(), y_score.ravel())
    roc_array = np.vstack((fpr, tpr)).T
    np.savetxt(f"rocdata_celltype_k{k}.txt", roc_array, header="FPR\tTPR", fmt="%.6f", delimiter="\t")
    with open("auc_celltype.txt", "a") as f:
        f.write(f"{k}\t{auc_score:.4f}\n")

# Combined ROC plot
plt.figure(figsize=(8, 6))
for k, y_test_bin, y_score in all_roc:
    auc_score = roc_auc_score(y_test_bin, y_score, average="macro", multi_class="ovr")
    RocCurveDisplay.from_predictions(
        y_test_bin.ravel(), y_score.ravel(),
        name=f"k={k} (AUC={auc_score:.2f})", ax=plt.gca()
    )
plt.title("Multiclass ROC Curves across k-mer sizes")
plt.tight_layout()
plt.savefig("celltype_all_roc.png")
