# feature_extraction.py

from Module import *

def gc_content(seq):
    gc = sum(1 for b in seq if b in "GCgc")
    return gc / len(seq) if len(seq) else 0

def shannon_entropy(seq):
    base_counts = Counter(seq.upper())
    total = sum(base_counts.values())
    probs = [count / total for count in base_counts.values()]
    return -sum(p * math.log2(p) for p in probs if p > 0)

def kmer_frequencies(seq, k):
    seq = seq.upper()
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    total = sum(counts.values())
    return {f'kmer_{k}_{kmer}': counts[kmer]/total for kmer in counts}

def extract_features(fasta_files, k):
    from pathlib import Path
    import pandas as pd

    data = []
    for cell_type, path in fasta_files.items():
        for record in SeqIO.parse(path, "fasta"):
            seq = str(record.seq)
            features = {
                "sequence_id": record.id,
                "cell_type": cell_type,
                "length": len(seq),
                "gc_content": gc_content(seq),
                "entropy": shannon_entropy(seq)
            }
            features.update(kmer_frequencies(seq, k=k))
            data.append(features)
    return pd.DataFrame(data).fillna(0)
