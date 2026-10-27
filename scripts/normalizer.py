"""Circular sequence linearization, SSC flip-flop resolution, and quadripartite partition extraction."""

from pathlib import Path
import json
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import (
    DATA_DIR,
    PARTITIONS_DIR,
    PROCESSED_DIR,
    RAW_FASTA_DIR,
    RAW_GB_DIR,
)
from scripts.ir_detector import IRJunctions, segment_plastome


def extract_quadripartite_sequences(record: SeqRecord, junctions: IRJunctions) -> tuple[str, str, str, str]:
    """Extracts raw sequence substrings for (LSC, IR_B, SSC, IR_A)."""
    seq = str(record.seq).upper()
    n = len(seq)

    # Extract IR_B and IR_A
    irb_seq = seq[junctions.irb_start : junctions.irb_end]
    ira_seq = seq[junctions.ira_start : junctions.ira_end]

    # Extract LSC (handle circular wrap if split across origin)
    if junctions.lsc_start < junctions.lsc_end:
        lsc_seq = seq[junctions.lsc_start : junctions.lsc_end]
    else:
        lsc_seq = seq[junctions.lsc_start :] + seq[: junctions.lsc_end]

    # Extract SSC (handle circular wrap if split across origin)
    if junctions.ssc_start < junctions.ssc_end:
        ssc_seq = seq[junctions.ssc_start : junctions.ssc_end]
    else:
        ssc_seq = seq[junctions.ssc_start :] + seq[: junctions.ssc_end]

    return lsc_seq, irb_seq, ssc_seq, ira_seq


def check_and_fix_ssc_orientation(reference_ssc: str, target_ssc: str) -> tuple[str, bool]:
    """Compares target SSC orientation against reference using anchor k-mers.
    
    Returns:
        (rectified_ssc_seq, was_flipped)
    """
    # Sample anchor seeds from 5' end of reference SSC (typically ndhF region)
    seed = reference_ssc[:120]
    rc_seed = str(Seq(seed).reverse_complement())

    # Count matching occurrences or substrings
    pos_forward = target_ssc.find(seed[:40])
    pos_rc = target_ssc.find(rc_seed[:40])

    if pos_rc != -1 and pos_forward == -1:
        # Target SSC is in flipped isomeric phase
        rectified = str(Seq(target_ssc).reverse_complement())
        return rectified, True
    return target_ssc, False


def normalize_dataset() -> dict[str, dict]:
    """Loads all curated genomes, standardizes linear phase, resolves SSC flip-flop,
    and writes out partitioned multi-FASTA files for LSC, SSC, and single-copy IR_B.
    """
    summary_path = DATA_DIR / "tier1_metadata_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Tier 1 summary not found at {summary_path}")

    with open(summary_path, "r", encoding="utf-8") as f:
        meta_list = json.load(f)

    # Reference genome is Musa_acuminata (HF677508)
    ref_gb = RAW_GB_DIR / "HF677508.gb"
    if not ref_gb.is_file():
        raise FileNotFoundError("Reference GenBank file HF677508.gb missing. Run curation.py first.")

    ref_record = SeqIO.read(ref_gb, "genbank")
    ref_junc = segment_plastome(ref_record)
    ref_lsc, ref_irb, ref_ssc, ref_ira = extract_quadripartite_sequences(ref_record, ref_junc)

    records_by_partition: dict[str, list[SeqRecord]] = {
        "LSC": [],
        "SSC": [],
        "IR_B": [],
        "FULL_NORMALIZED": [],
    }

    normalized_info = {}
    print("\n--- Structural Normalization & SSC Flip-Flop Resolution ---")

    for meta in meta_list:
        clean_acc = meta["accession"].split(".")[0]
        gb_file = RAW_GB_DIR / f"{clean_acc}.gb"
        if not gb_file.is_file():
            continue

        record = SeqIO.read(gb_file, "genbank")
        taxon_id = meta["taxon"].replace(" ", "_").replace(".", "")
        junc = segment_plastome(record)
        lsc_seq, irb_seq, ssc_seq, ira_seq = extract_quadripartite_sequences(record, junc)

        # Check and resolve SSC flip-flop relative to reference
        rectified_ssc, was_flipped = check_and_fix_ssc_orientation(ref_ssc, ssc_seq)
        flip_status = "FLIPPED (Corrected)" if was_flipped else "CONCORDANT"

        # Assemble normalized collinear pseudo-molecule: LSC + IR_B + SSC + revcomp(IR_B)
        # Note: revcomp(IR_B) ensures exact quadripartite continuity
        reconstituted_full = lsc_seq + irb_seq + rectified_ssc + str(Seq(irb_seq).reverse_complement())

        normalized_info[clean_acc] = {
            "taxon": taxon_id,
            "section": meta["section"],
            "lsc_len": len(lsc_seq),
            "irb_len": len(irb_seq),
            "ssc_len": len(rectified_ssc),
            "full_len": len(reconstituted_full),
            "ssc_flip": was_flipped,
        }

        print(f"[{clean_acc:<12}] {taxon_id:<24} | SSC Phase: {flip_status:<20} | Length: {len(reconstituted_full)} bp")

        # Create SeqRecords for partitioned multi-FASTA
        records_by_partition["LSC"].append(SeqRecord(Seq(lsc_seq), id=taxon_id, description="LSC partition"))
        records_by_partition["SSC"].append(SeqRecord(Seq(rectified_ssc), id=taxon_id, description="SSC partition"))
        records_by_partition["IR_B"].append(SeqRecord(Seq(irb_seq), id=taxon_id, description="IR_B single copy partition"))
        records_by_partition["FULL_NORMALIZED"].append(SeqRecord(Seq(reconstituted_full), id=taxon_id, description="Collinear normalized plastome"))

    # Write out partitioned multi-FASTAs
    for part_name, rec_list in records_by_partition.items():
        out_fasta = PARTITIONS_DIR / f"tier1_{part_name}.fasta"
        SeqIO.write(rec_list, out_fasta, "fasta")
        print(f"Saved {len(rec_list)} records to: {out_fasta}")

    # Save normalization summary
    summary_out = DATA_DIR / "tier1_normalization_summary.json"
    with open(summary_out, "w", encoding="utf-8") as f:
        json.dump(normalized_info, f, indent=2)

    return normalized_info


if __name__ == "__main__":
    normalize_dataset()
