"""Exact mathematical Inverted Repeat (IR) boundary detection and quadripartite segmentation."""

from dataclasses import asdict, dataclass
from pathlib import Path
import json
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import DATA_DIR, RAW_FASTA_DIR, RAW_GB_DIR


@dataclass
class IRJunctions:
    accession: str
    genome_length: int
    lsc_start: int
    lsc_end: int
    lsc_len: int
    irb_start: int
    irb_end: int
    irb_len: int
    ssc_start: int
    ssc_end: int
    ssc_len: int
    ira_start: int
    ira_end: int
    ira_len: int
    j_lb_gene: str = "unknown"
    j_sb_gene: str = "unknown"
    j_sa_gene: str = "unknown"
    j_la_gene: str = "unknown"


def find_exact_ir_boundaries(seq: str, min_ir_len: int = 15_000, seed_size: int = 50) -> tuple[tuple[int, int], tuple[int, int]]:
    """Pinpoints exact [IR_B, IR_A] coordinates using a seed-and-extend reverse-complement search.
    
    Returns:
        ((ir1_start, ir1_end), (ir2_start, ir2_end)) in 0-based indexing.
    """
    seq = seq.upper()
    n = len(seq)
    rc_map = str.maketrans("ACGTURYKMSWBDHVN", "TGCAAYRMKSWVHDBN")

    # Sample seeds and locate reverse-complement matches
    hit = None
    for i in range(0, n - seed_size, 200):
        seed_rc = seq[i : i + seed_size].translate(rc_map)[::-1]
        pos = seq.find(seed_rc)
        if pos != -1 and abs(pos - i) >= min_ir_len:
            hit = (i, pos)
            break

    if not hit:
        raise ValueError(f"Could not locate Inverted Repeat with length >= {min_ir_len} bp in genome of length {n} bp.")

    i, j = hit
    k = seed_size

    # In forward sequence: s[i + offset] == revcomp(s[j + k - 1 - offset])
    # Extend right from seed: i moves +1, j moves -1
    ext_r = k
    cur_j = j
    while (i + ext_r < n) and (cur_j - 1 >= 0):
        char_i = seq[i + ext_r]
        char_j_rc = seq[cur_j - 1].translate(rc_map)
        if char_i == char_j_rc:
            ext_r += 1
            cur_j -= 1
        else:
            break

    # Extend left from seed: i moves -1, partner_end moves +1
    partner_end = j + k
    ext_l = 0
    while (i - ext_l - 1 >= 0) and (partner_end < n):
        char_i = seq[i - ext_l - 1]
        char_partner_rc = seq[partner_end].translate(rc_map)
        if char_i == char_partner_rc:
            ext_l += 1
            partner_end += 1
        else:
            break

    start1 = i - ext_l
    end1 = i + ext_r
    start2 = cur_j
    end2 = partner_end

    # Ensure canonical linear order (start1 < start2)
    if start1 > start2:
        start1, end1, start2, end2 = start2, end2, start1, end1

    return (start1, end1), (start2, end2)


def segment_plastome(record: SeqRecord) -> IRJunctions:
    """Segments a full plastome into LSC, IRb, SSC, IRa with 0-based coordinates."""
    seq_str = str(record.seq).upper()
    n = len(seq_str)

    (ir1_start, ir1_end), (ir2_start, ir2_end) = find_exact_ir_boundaries(seq_str)
    ir_len = ir1_end - ir1_start

    # Between ir1 and ir2 there is a single-copy spacer (spacer 1)
    # Between ir2 and the end/beginning of circular genome is another single-copy spacer (spacer 2)
    spacer1_len = ir2_start - ir1_end
    spacer2_len = (n - ir2_end) + ir1_start

    # Determine which spacer is LSC and which is SSC
    # LSC is canonically ~80-95 kb in Musa; SSC is ~10-15 kb
    if spacer1_len > spacer2_len:
        # spacer 1 is LSC
        lsc_start = ir1_end
        lsc_end = ir2_start
        lsc_len = spacer1_len

        # ir2 is IR_B, preceding SSC
        irb_start = ir2_start
        irb_end = ir2_end
        irb_len = ir_len

        # spacer 2 is SSC (wrap around or at the ends)
        # Note: if spacer 2 is split across the origin, record this
        ssc_start = ir2_end
        ssc_end = ir1_start  # circular wrap
        ssc_len = spacer2_len

        ira_start = ir1_start
        ira_end = ir1_end
        ira_len = ir_len
    else:
        # spacer 2 is LSC (wrapped around origin)
        lsc_start = ir2_end
        lsc_end = ir1_start  # circular wrap
        lsc_len = spacer2_len

        # ir1 is IR_B, preceding spacer 1 (which is SSC)
        irb_start = ir1_start
        irb_end = ir1_end
        irb_len = ir_len

        ssc_start = ir1_end
        ssc_end = ir2_start
        ssc_len = spacer1_len

        ira_start = ir2_start
        ira_end = ir2_end
        ira_len = ir_len

    # Identify flanking junction genes if GenBank features are present
    j_genes = {"j_lb": "unknown", "j_sb": "unknown", "j_sa": "unknown", "j_la": "unknown"}
    for feat in record.features:
        if feat.type == "CDS" and "gene" in feat.qualifiers:
            g_name = feat.qualifiers["gene"][0]
            f_start = int(feat.location.start)
            f_end = int(feat.location.end)

            if abs(f_end - irb_start) < 250 or abs(f_start - irb_start) < 250:
                j_genes["j_lb"] = g_name
            if abs(f_start - irb_end) < 250 or abs(f_end - irb_end) < 250:
                j_genes["j_sb"] = g_name
            if abs(f_end - ira_start) < 250 or abs(f_start - ira_start) < 250:
                j_genes["j_sa"] = g_name
            if abs(f_start - ira_end) < 250 or abs(f_end - ira_end) < 250:
                j_genes["j_la"] = g_name

    return IRJunctions(
        accession=record.id,
        genome_length=n,
        lsc_start=lsc_start,
        lsc_end=lsc_end,
        lsc_len=lsc_len,
        irb_start=irb_start,
        irb_end=irb_end,
        irb_len=irb_len,
        ssc_start=ssc_start,
        ssc_end=ssc_end,
        ssc_len=ssc_len,
        ira_start=ira_start,
        ira_end=ira_end,
        ira_len=ira_len,
        j_lb_gene=j_genes["j_lb"],
        j_sb_gene=j_genes["j_sb"],
        j_sa_gene=j_genes["j_sa"],
        j_la_gene=j_genes["j_la"],
    )


def process_all_curated_junctions() -> list[IRJunctions]:
    """Processes all Tier 1 GenBank files and calculates exact IR boundaries."""
    summary_path = DATA_DIR / "tier1_metadata_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Tier 1 summary not found at {summary_path}. Run curation.py first.")

    with open(summary_path, "r", encoding="utf-8") as f:
        meta_list = json.load(f)

    junctions_list = []
    print("\n--- Computing Exact IR Boundaries & Quadripartite Partitions ---")
    for meta in meta_list:
        clean_acc = meta["accession"].split(".")[0]
        gb_file = RAW_GB_DIR / f"{clean_acc}.gb"
        if not gb_file.is_file():
            continue
        record = SeqIO.read(gb_file, "genbank")
        junc = segment_plastome(record)
        junctions_list.append(junc)
        print(
            f"[{junc.accession:<12}] LSC: {junc.lsc_len:>6} bp | IR_B: {junc.irb_len:>5} bp | "
            f"SSC: {junc.ssc_len:>5} bp | IR_A: {junc.ira_len:>5} bp | Total: {junc.genome_length:>6} bp"
        )

    # Save junctions
    out_path = DATA_DIR / "tier1_ir_junctions.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump([asdict(j) for j in junctions_list], f, indent=2)
    print(f"\nSaved IR junctions to: {out_path}")
    return junctions_list


if __name__ == "__main__":
    process_all_curated_junctions()
