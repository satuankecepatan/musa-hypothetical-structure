"""Sliding-window comparative genomics: nucleotide diversity (pi), Shannon entropy (H), and GC skew."""

from dataclasses import asdict, dataclass
from pathlib import Path
import json
import sys
from Bio import AlignIO, SeqIO
import numpy as np
import pandas as pd
from scipy import stats

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import DATA_DIR, OUTPUT_DIR, PROCESSED_DIR


@dataclass
class WindowMetric:
    midpoint: int
    region: str
    pi: float
    entropy: float
    gc_skew: float
    gc_content: float
    gap_freq: float


def calculate_shannon_entropy(column: list[str]) -> float:
    """Calculates Shannon entropy H(x) = -sum(p_i * log2(p_i)) for a column."""
    counts = {}
    total = len(column)
    for char in column:
        counts[char] = counts.get(char, 0) + 1
    
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * np.log2(p)
    return float(entropy)


def calculate_pairwise_pi(column: list[str]) -> float:
    """Calculates Nei & Li nucleotide diversity (pi) for a single alignment site."""
    valid_bases = [b for b in column if b in "ACGT"]
    n = len(valid_bases)
    if n < 2:
        return 0.0
    
    diff_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if valid_bases[i] != valid_bases[j]:
                diff_count += 1
    total_pairs = (n * (n - 1)) / 2
    return float(diff_count / total_pairs)


def compute_sliding_window_metrics(
    alignment_file: Path,
    manifest_path: Path,
    window_size: int = 600,
    step_size: int = 100,
) -> pd.DataFrame:
    """Computes sliding-window genomic diversity metrics across the aligned plastomes."""
    alignment = AlignIO.read(str(alignment_file), "fasta")
    aln_len = alignment.get_alignment_length()
    n_taxa = len(alignment)

    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    # Precompute per-site metrics
    site_entropy = np.zeros(aln_len)
    site_pi = np.zeros(aln_len)
    site_g = np.zeros(aln_len)
    site_c = np.zeros(aln_len)
    site_gaps = np.zeros(aln_len)

    for i in range(aln_len):
        col = [alignment[t, i].upper() for t in range(n_taxa)]
        site_entropy[i] = calculate_shannon_entropy(col)
        site_pi[i] = calculate_pairwise_pi(col)
        site_g[i] = col.count("G")
        site_c[i] = col.count("C")
        site_gaps[i] = col.count("-") / n_taxa

    # Sliding window sweep
    windows = []
    lsc_end = manifest["lsc"]["length"]
    irb_end = lsc_end + manifest["irb"]["length"]
    ssc_end = irb_end + manifest["ssc"]["length"]

    for start in range(0, aln_len - window_size + 1, step_size):
        end = start + window_size
        mid = (start + end) // 2

        # Assign region
        if mid < lsc_end:
            region = "LSC"
        elif mid < irb_end:
            region = "IR_B"
        elif mid < ssc_end:
            region = "SSC"
        else:
            region = "IR_A"

        win_pi = float(np.mean(site_pi[start:end]))
        win_entropy = float(np.mean(site_entropy[start:end]))
        total_g = np.sum(site_g[start:end])
        total_c = np.sum(site_c[start:end])
        gc_skew = float((total_g - total_c) / (total_g + total_c)) if (total_g + total_c) > 0 else 0.0
        gc_content = float((total_g + total_c) / (window_size * n_taxa))
        gap_freq = float(np.mean(site_gaps[start:end]))

        windows.append({
            "midpoint": mid,
            "region": region,
            "pi": round(win_pi, 5),
            "entropy": round(win_entropy, 5),
            "gc_skew": round(gc_skew, 5),
            "gc_content": round(gc_content, 5),
            "gap_freq": round(gap_freq, 5),
        })

    df = pd.DataFrame(windows)
    
    # Statistical tests across regions
    print("\n--- Non-Parametric Tests across Quadripartite Compartments ---")
    lsc_pi = df[df["region"] == "LSC"]["pi"]
    ssc_pi = df[df["region"] == "SSC"]["pi"]
    ir_pi = df[df["region"].isin(["IR_B", "IR_A"])]["pi"]

    # Mann-Whitney U test between Single Copy vs IR
    stat_u, p_val = stats.mannwhitneyu(lsc_pi, ir_pi, alternative="greater")
    print(f"[Mann-Whitney U] LSC vs IR nucleotide diversity pi: U = {stat_u:.1f}, p = {p_val:.2e}")

    stat_kruskal, p_kruskal = stats.kruskal(lsc_pi, ssc_pi, ir_pi)
    print(f"[Kruskal-Wallis] Across LSC, SSC, and IR: H = {stat_kruskal:.1f}, p = {p_kruskal:.2e}")

    # Identify top 5% divergence hotspots
    pi_threshold = df["pi"].quantile(0.95)
    hotspots = df[df["pi"] >= pi_threshold]
    print(f"\n[Divergence Hotspots] Top 5% pi threshold = {pi_threshold:.5f} ({len(hotspots)} windows)")

    # Save metrics
    metrics_path = DATA_DIR / "musa_sliding_window_metrics.csv"
    df.to_csv(metrics_path, index=False)
    print(f"Saved sliding window metrics to: {metrics_path}")

    return df


def create_full_concatenated_alignment() -> Path:
    """Concatenates aligned partitions (LSC, IR_B, SSC, revcomp(IR_B)) into a full alignment."""
    full_aln_path = PROCESSED_DIR / "aligned_FULL.fasta"
    if full_aln_path.is_file() and full_aln_path.stat().st_size > 1000:
        return full_aln_path

    lsc_aln = AlignIO.read(str(PROCESSED_DIR / "aligned_LSC.fasta"), "fasta")
    irb_aln = AlignIO.read(str(PROCESSED_DIR / "aligned_IR_B.fasta"), "fasta")
    ssc_aln = AlignIO.read(str(PROCESSED_DIR / "aligned_SSC.fasta"), "fasta")

    # Map records by taxon ID
    lsc_map = {r.id: str(r.seq) for r in lsc_aln}
    irb_map = {r.id: str(r.seq) for r in irb_aln}
    ssc_map = {r.id: str(r.seq) for r in ssc_aln}

    concatenated_records = []
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rc_map = str.maketrans("ACGTURYKMSWBDHVN-", "TGCAAYRMKSWVHDBN-")

    for r in lsc_aln:
        tid = r.id
        l_seq = lsc_map[tid]
        b_seq = irb_map[tid]
        s_seq = ssc_map[tid]
        a_seq = b_seq.translate(rc_map)[::-1]
        full_seq = l_seq + b_seq + s_seq + a_seq
        concatenated_records.append(SeqRecord(Seq(full_seq), id=tid, description="Concatenated aligned plastome"))

    SeqIO.write(concatenated_records, full_aln_path, "fasta")
    print(f"[Analytics] Created full quadripartite alignment: {full_aln_path} ({len(concatenated_records[0].seq)} columns)")
    return full_aln_path


if __name__ == "__main__":
    aln_file = create_full_concatenated_alignment()
    manifest_file = DATA_DIR / "musa_mrca_manifest.json"
    if aln_file.is_file() and manifest_file.is_file():
        compute_sliding_window_metrics(aln_file, manifest_file)
    else:
        print("Prerequisite files missing.")
