"""Feature annotation lift-over, trans-splicing resolution, and Table 11 translation audit."""

from pathlib import Path
import json
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import (
    DATA_DIR,
    OUTPUT_DIR,
    PROCESSED_DIR,
    RAW_GB_DIR,
)


def align_and_locate_feature(feature_seq: str, target_genome: str) -> tuple[int, int, int] | None:
    """Finds optimal occurrence of a feature sequence in target genome.
    
    Returns:
        (start, end, strand) or None
    """
    feature_seq = feature_seq.upper()
    target_genome = target_genome.upper()
    
    # Try exact forward match
    pos = target_genome.find(feature_seq)
    if pos != -1:
        return pos, pos + len(feature_seq), +1

    # Try exact reverse complement match
    rc_feat = str(Seq(feature_seq).reverse_complement())
    pos_rc = target_genome.find(rc_feat)
    if pos_rc != -1:
        return pos_rc, pos_rc + len(rc_feat), -1

    # Try seed match (first 100 bp or full length if shorter)
    k = min(80, len(feature_seq))
    seed = feature_seq[:k]
    pos = target_genome.find(seed)
    if pos != -1:
        return pos, pos + len(feature_seq), +1

    rc_seed = str(Seq(feature_seq[-k:]).reverse_complement())
    pos_rc = target_genome.find(rc_seed)
    if pos_rc != -1:
        return pos_rc, pos_rc + len(feature_seq), -1

    return None


def annotate_ancestral_plastome() -> tuple[Path, dict]:
    """Transfers annotations from reference Musa acuminata to ancestral plastome,
    resolving trans-spliced rps12, boundary pseudogenes, and auditing Table 11 ORFs.
    """
    ref_gb = RAW_GB_DIR / "HF677508.gb"
    anc_fasta = DATA_DIR / "musa_mrca_plastome.fasta"

    if not ref_gb.is_file() or not anc_fasta.is_file():
        raise FileNotFoundError("Prerequisite files HF677508.gb or musa_mrca_plastome.fasta missing.")

    ref_record = SeqIO.read(ref_gb, "genbank")
    anc_record = SeqIO.read(anc_fasta, "fasta")
    anc_seq_str = str(anc_record.seq).upper()
    total_len = len(anc_seq_str)

    new_features: list[SeqFeature] = []
    audit_results: list[dict] = []

    # Manifest boundaries
    manifest_path = DATA_DIR / "musa_mrca_manifest.json"
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    irb_start = manifest["irb"]["start"]
    irb_end = manifest["irb"]["end"]
    ira_start = manifest["ira"]["start"]
    ira_end = manifest["ira"]["end"]

    print("\n--- Annotating Inferred Ancestral Plastome & ORF Audit ---")

    for feat in ref_record.features:
        if feat.type not in ["CDS", "tRNA", "rRNA"]:
            continue

        gene_name = feat.qualifiers.get("gene", ["unnamed"])[0]

        # Extract sequence of feature from reference
        try:
            feat_seq = str(feat.extract(ref_record.seq)).upper()
        except Exception:
            continue

        # Handle trans-spliced rps12 explicitly
        if gene_name == "rps12" and isinstance(feat.location, CompoundLocation):
            # Exon 1 in LSC, Exons 2-3 in IRs
            locs = []
            for part in feat.location.parts:
                part_seq = str(ref_record.seq[part.start : part.end]).upper()
                hit = align_and_locate_feature(part_seq, anc_seq_str)
                if hit:
                    locs.append(FeatureLocation(hit[0], hit[1], strand=hit[2]))
            if len(locs) >= 2:
                new_loc = CompoundLocation(locs)
                new_feat = SeqFeature(new_loc, type=feat.type, qualifiers=dict(feat.qualifiers))
                new_features.append(new_feat)
                continue

        # Single-location features
        hit = align_and_locate_feature(feat_seq, anc_seq_str)
        if not hit:
            continue

        start, end, strand = hit
        new_loc = FeatureLocation(start, end, strand=strand)
        new_feat = SeqFeature(new_loc, type=feat.type, qualifiers=dict(feat.qualifiers))

        # Check for boundary-induced pseudogene psi-ycf1
        is_pseudogene = False
        if gene_name == "ycf1" and (len(feat_seq) < 3000):
            new_feat.qualifiers["pseudogene"] = ["unclassified"]
            new_feat.qualifiers["note"] = ["Truncated boundary-induced pseudogene (psi-ycf1) at IR junction"]
            is_pseudogene = True

        # ORF Integrity Audit for CDS
        orf_valid = True
        orf_notes = []
        if feat.type == "CDS" and not is_pseudogene:
            cds_seq = Seq(str(new_feat.extract(anc_record.seq)).upper())
            
            # Check length multiple of 3
            if len(cds_seq) % 3 != 0:
                orf_valid = False
                orf_notes.append(f"Non-triplet length ({len(cds_seq)} bp)")

            # Check start codon (allow ACG for ndhD RNA-editing)
            start_codon = str(cds_seq[:3])
            canonical_starts = {"ATG", "GTG", "TTG"}
            if gene_name.lower() == "ndhd" and start_codon == "ACG":
                orf_notes.append("Canonical C-to-U RNA-edited start codon (ACG -> AUG)")
            elif start_codon not in canonical_starts:
                orf_notes.append(f"Atypical start codon ({start_codon})")

            # Translation using Table 11
            try:
                translation = cds_seq.translate(table=11)
                # Check for premature internal stops
                internal_stops = translation[:-1].count("*")
                if internal_stops > 0:
                    orf_valid = False
                    orf_notes.append(f"{internal_stops} premature stop codons")
            except Exception as e:
                orf_valid = False
                orf_notes.append(f"Translation error: {e}")

            audit_results.append({
                "gene": gene_name,
                "length": len(cds_seq),
                "start": start,
                "end": end,
                "strand": strand,
                "valid_orf": orf_valid,
                "notes": "; ".join(orf_notes) if orf_notes else "Intact functional ORF",
            })

        new_features.append(new_feat)

    # Attach features to record
    anc_record.features = new_features
    anc_record.annotations = {
        "molecule_type": "DNA",
        "topology": "circular",
        "data_file_division": "PLN",
        "organism": "Musa hypothetical ancestor",
        "taxonomy": ["Eukaryota", "Viridiplantae", "Streptophyta", "Liliopsida", "Zingiberales", "Musaceae", "Musa"],
    }

    # Save GenBank record
    data_gbk = DATA_DIR / "musa_mrca_plastome.gbk"
    public_gbk = OUTPUT_DIR / "musa_mrca_plastome.gbk"
    SeqIO.write(anc_record, data_gbk, "genbank")
    SeqIO.write(anc_record, public_gbk, "genbank")

    # Save ORF audit report
    audit_path = DATA_DIR / "musa_mrca_orf_audit.json"
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_results, f, indent=2)

    passed_count = sum(1 for a in audit_results if a["valid_orf"])
    print(f"[Annotation] Successfully mapped {len(new_features)} features.")
    print(f"[ORF Audit]  Intact coding ORFs: {passed_count}/{len(audit_results)}")
    print(f"Saved annotated GenBank to: {data_gbk}")

    # Generate PGA and GeSeq submission package
    export_validation_packages(anc_record)

    return data_gbk, audit_results


def export_validation_packages(record: SeqRecord) -> None:
    """Exports structured directories for external validation via PGA and GeSeq."""
    pga_dir = DATA_DIR / "pga_validation"
    geseq_dir = DATA_DIR / "geseq_validation"
    pga_dir.mkdir(parents=True, exist_ok=True)
    geseq_dir.mkdir(parents=True, exist_ok=True)

    # Save target FASTA
    SeqIO.write(record, pga_dir / "target_ancestral.fasta", "fasta")
    SeqIO.write(record, geseq_dir / "target_ancestral.fasta", "fasta")

    # Copy reference GenBank into PGA reference folder
    ref_gb = RAW_GB_DIR / "HF677508.gb"
    if ref_gb.is_file():
        with open(ref_gb, "r", encoding="utf-8") as src, open(pga_dir / "reference.gb", "w", encoding="utf-8") as dst:
            dst.write(src.read())

    print(f"[Export] Prepared PGA validation bundle at: {pga_dir}")
    print(f"[Export] Prepared GeSeq submission bundle at: {geseq_dir}")


if __name__ == "__main__":
    annotate_ancestral_plastome()
