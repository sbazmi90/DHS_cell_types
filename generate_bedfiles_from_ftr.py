import argparse
import gzip
import os
import shutil
from pathlib import Path
import marimo as mo

import gdown
import numpy as np
import pandas as pd
import requests
from Bio import SeqIO



# Load master dataset
df = pd.read_feather("master_dataset.ftr")

# Check if target column is present
cell_type = "GM12878_ENCLB441ZZZ"
assert cell_type in df.columns, f"{cell_type} not found in dataset"

# Filter rows where GM12878 is active
active_dhs = df[df[cell_type] > 0].copy()

# Save filtered DHSs for GM12878
active_dhs.to_csv("GM12878_active_DHS.txt", sep="\t", index=False)
print(f"{len(active_dhs)} DHSs found for {cell_type} and saved to GM12878_active_DHS.txt")


# Load filtered DHSs for GM12878
df_active = pd.read_csv("GM12878_active_DHS.txt", sep="\t")

# Load DHS index (vocabulary)
vocab = pd.read_csv("DHS_Index_and_Vocabulary_hg38_WM20190703.txt", sep="\t")

print("Active DHS columns:\n", df_active.columns.tolist())
print("Vocabulary columns:\n", vocab.columns.tolist())

# Ensure all column names are clean
vocab.columns = vocab.columns.str.strip()
df_active.columns = df_active.columns.str.strip()

# Rename for merge
df_active = df_active.rename(columns={"dhs_id": "DHS_ID"})
vocab = vocab.rename(columns={"identifier": "DHS_ID", "seqname": "chr"})

# Check which coordinate columns are actually present
print("Available coordinate columns in vocab:", vocab.columns.tolist())

# Merge
merged = df_active.merge(vocab, on="DHS_ID")

# Print sample row to inspect column structure
print(merged.head(1))

# Try to select coordinate columns
try:
    bed_df = merged[["chr", "start", "end", "DHS_ID"]]
    bed_df.to_csv("GM12878_DHS.bed", sep="\t", index=False, header=False)
    print(f"{len(bed_df)} DHS coordinates written to GM12878_DHS.bed")
except KeyError:
    print("ERROR: Coordinate columns not found in merged dataframe.")


print("Sample DHS_IDs from df_active:")
print(df_active["DHS_ID"].astype(str).head())

print("\nSample DHS_IDs from vocab:")
print(vocab["DHS_ID"].astype(str).head())


df_active["DHS_ID"] = df_active["DHS_ID"].astype(str)
vocab["DHS_ID"] = vocab["DHS_ID"].astype(str)

merged = df_active.merge(vocab, on="DHS_ID")


# Step 1: Standardize and clean column names
df_active.columns = df_active.columns.str.strip()
vocab.columns = vocab.columns.str.strip()

# Step 2: Merge on chr + start + end instead of DHS_ID
merged = pd.merge(
    df_active,
    vocab,
    how="inner",
    left_on=["chr", "start", "end"],
    right_on=["chr", "start", "end"]
)

# Step 3: Create BED-format dataframe
bed_df = merged[["chr", "start", "end"]].copy()
bed_df["name"] = ["DHS_" + str(i) for i in range(len(bed_df))]

# Step 4: Save to BED file
bed_df.to_csv("GM12878_DHS.bed", sep="\t", index=False, header=False)

print(f"{len(bed_df)} DHS regions written to GM12878_DHS.bed")
