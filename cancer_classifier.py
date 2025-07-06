from Module import *
from feature_extraction import extract_features

# Input FASTA files
fasta_files = {
    "GM12878": "GM12878_DHS_sequences.fa",
    "HepG2": "HepG2_DHS_sequences.fa",
    "hESC": "hESCT0_DHS_sequences.fa",
    "K562": "K562_DHS_sequences.fa"
}
cancer_map = {"GM12878": 0, "hESC": 0, "HepG2": 1, "K562": 1}
kmers_to_try = [5]
all_roc = []

# Reset summary AUC file
with open("auc_cancer.txt", "w") as f:
    f.write("k-mer\tAUC\n")

for k in kmers_to_try:
    print(f"\n=== Processing k={k} ===")
    df = extract_features(fasta_files, k)
    df["cancer_label"] = df["cell_type"].map(cancer_map)

    df_cancer = df[df.cancer_label == 1]
    df_normal = df[df.cancer_label == 0]
    n = min(len(df_cancer), len(df_normal))

    df_balanced = pd.concat([
        resample(df_cancer, n_samples=n, random_state=42),
        resample(df_normal, n_samples=n, random_state=42)
    ]).sample(frac=1, random_state=42)

    X = df_balanced.drop(columns=["sequence_id", "cell_type", "cancer_label"])
    y = df_balanced["cancer_label"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=0.2, random_state=42)

    clf = RandomForestClassifier(n_estimators=100, class_weight="balanced", random_state=42)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=[0, 1])
    plt.figure(figsize=(5, 4))
    sns.heatmap(cm, annot=True, fmt='d', cmap="Reds",
                xticklabels=["Non-Cancer", "Cancer"],
                yticklabels=["Non-Cancer", "Cancer"])
    plt.title(f"Cancer Confusion Matrix (k={k})")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.tight_layout()
    plt.savefig(f"cancer_confusion_matrix_k{k}.png")
    plt.close()

    # ROC & AUC
    y_proba = clf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_proba)
    print(f"Binary AUC (k={k}): {auc:.4f}")
    all_roc.append((k, y_test, y_proba))

    # Save ROC curve data
    fpr, tpr, _ = roc_curve(y_test, y_proba)
    roc_array = np.vstack((fpr, tpr)).T
    np.savetxt(f"rocdata_cancer_k{k}.txt", roc_array, header="FPR\tTPR", fmt="%.6f", delimiter="\t")
    with open("auc_cancer.txt", "a") as f:
        f.write(f"{k}\t{auc:.4f}\n")

# Combined ROC plot
plt.figure(figsize=(8, 6))
for k, y_test, y_proba in all_roc:
    auc_score = roc_auc_score(y_test, y_proba)
    RocCurveDisplay.from_predictions(
        y_test, y_proba, name=f"k={k} (AUC={auc_score:.2f})", ax=plt.gca()
    )
plt.title("ROC Curves: Cancer vs Non-Cancer across k-mer sizes")
plt.tight_layout()
plt.savefig("cancer_all_roc.png")
