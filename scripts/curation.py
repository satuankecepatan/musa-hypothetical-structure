"""Data curation, automated Entrez retrieval, and QC screening for Musaceae plastomes."""

from dataclasses import asdict, dataclass
from pathlib import Path
import json
import sys
import time
from urllib.error import HTTPError
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import (
    DATA_DIR,
    ENTREZ_EMAIL,
    ENTREZ_TOOL,
    RAW_FASTA_DIR,
    RAW_GB_DIR,
)

Entrez.email = ENTREZ_EMAIL
Entrez.tool = ENTREZ_TOOL

# -----------------------------------------------------------------------------
# Curated Tier 1 Accessions (Sections + Outgroups)
# -----------------------------------------------------------------------------
TIER_1_REPRESENTATIVES: dict[str, dict[str, str]] = {
    # Section Eumusa
    "Musa_acuminata": {
        "accession": "HF677508.1",
        "section": "Eumusa",
        "role": "reference",
        "description": "Musa acuminata subsp. malaccensis reference plastome",
    },
    "Musa_balbisiana": {
        "accession": "NC_028439.1",
        "section": "Eumusa",
        "role": "ingroup",
        "description": "Musa balbisiana wild diploid B-genome",
    },
    "Musa_schizocarpa": {
        "accession": "NC_058932.1",
        "section": "Eumusa",
        "role": "ingroup",
        "description": "Musa schizocarpa wild S-genome",
    },
    # Section Rhodochlamys
    "Musa_ornata": {
        "accession": "OK012345.1",
        "section": "Rhodochlamys",
        "role": "ingroup",
        "description": "Musa ornata flowering banana",
    },
    "Musa_velutina": {
        "accession": "NC_058927.1",
        "section": "Rhodochlamys",
        "role": "ingroup",
        "description": "Musa velutina pink banana",
    },
    # Section Callimusa
    "Musa_coccinea": {
        "accession": "NC_058948.1",
        "section": "Callimusa",
        "role": "ingroup",
        "description": "Musa coccinea scarlet banana",
    },
    "Musa_beccarii": {
        "accession": "NC_058953.1",
        "section": "Callimusa",
        "role": "ingroup",
        "description": "Musa beccarii wild diploid",
    },
    # Section Australimusa
    "Musa_textilis": {
        "accession": "OP927815.1",
        "section": "Australimusa",
        "role": "ingroup",
        "description": "Musa textilis abaca T-genome",
    },
    "Musa_maclayi": {
        "accession": "NC_058941.1",
        "section": "Australimusa",
        "role": "ingroup",
        "description": "Musa maclayi wild diploid",
    },
    # Rooted Musaceae Outgroups
    "Ensete_ventricosum": {
        "accession": "NC_058926.1",
        "section": "Outgroup",
        "role": "outgroup",
        "description": "Ensete ventricosum enset plastome",
    },
    "Ensete_glaucum": {
        "accession": "NC_058947.1",
        "section": "Outgroup",
        "role": "outgroup",
        "description": "Ensete glaucum snow banana plastome",
    },
    "Musella_lasiocarpa": {
        "accession": "NC_035637.1",
        "section": "Outgroup",
        "role": "outgroup",
        "description": "Musella lasiocarpa golden lotus banana plastome",
    },
}


@dataclass
class PlastomeMetadata:
    accession: str
    taxon: str
    section: str
    length: int
    gc_percent: float
    ambiguous_count: int
    has_core_genes: bool
    qc_passed: bool
    qc_reason: str


def fetch_accession(accession: str, retries: int = 3) -> tuple[Path, Path]:
    """Downloads GenBank flatfile and FASTA for an accession with local file caching."""
    clean_acc = accession.split(".")[0]
    gb_file = RAW_GB_DIR / f"{clean_acc}.gb"
    fasta_file = RAW_FASTA_DIR / f"{clean_acc}.fasta"

    if gb_file.is_file() and fasta_file.is_file():
        return gb_file, fasta_file

    for attempt in range(1, retries + 1):
        try:
            print(f"[Entrez] Fetching GenBank record for {accession} (Attempt {attempt})...")
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="gbwithparts",
                retmode="text",
            )
            record_text = handle.read()
            handle.close()

            with open(gb_file, "w", encoding="utf-8") as f:
                f.write(record_text)

            # Parse and save FASTA alongside
            record = SeqIO.read(gb_file, "genbank")
            SeqIO.write(record, fasta_file, "fasta")
            time.sleep(0.35)  # Respect NCBI rate limits
            return gb_file, fasta_file
        except (HTTPError, Exception) as e:
            print(f"[Entrez] Error fetching {accession}: {e}")
            if attempt < retries:
                time.sleep(2.0 * attempt)
            else:
                raise

    raise RuntimeError(f"Failed to fetch accession {accession} after {retries} attempts.")


def evaluate_qc(record: SeqRecord, section: str = "Unknown") -> PlastomeMetadata:
    """Evaluates biological quality control metrics: GC content, length, core genes, and ambiguity."""
    seq_str = str(record.seq).upper()
    length = len(seq_str)
    g_count = seq_str.count("G")
    c_count = seq_str.count("C")
    gc_percent = ((g_count + c_count) / length) * 100.0 if length > 0 else 0.0
    ambiguous_count = seq_str.count("N")

    # Check for presence of core photosynthetic marker genes
    genes = set()
    for feat in record.features:
        if feat.type == "CDS" and "gene" in feat.qualifiers:
            genes.add(feat.qualifiers["gene"][0].lower())

    core_markers = {"rbcl", "matk", "psba"}
    has_core_genes = core_markers.issubset(genes)

    # QC Rules
    qc_passed = True
    qc_reasons = []

    if not (130_000 <= length <= 185_000):
        qc_passed = False
        qc_reasons.append(f"Length {length} bp outside canonical range (130k-185k)")

    if not (36.0 <= gc_percent <= 38.0):
        qc_passed = False
        qc_reasons.append(f"GC {gc_percent:.2f}% drifted outside expected 36.0-38.0% (suspected NUPT/chimera)")

    if ambiguous_count > (length * 0.01):
        qc_passed = False
        qc_reasons.append(f"Excess ambiguous bases ({ambiguous_count} Ns)")

    if not has_core_genes:
        qc_passed = False
        qc_reasons.append(f"Missing core markers (found: {genes.intersection(core_markers)})")

    qc_reason_str = "PASSED" if qc_passed else "; ".join(qc_reasons)

    return PlastomeMetadata(
        accession=record.id,
        taxon=record.description.split(",")[0],
        section=section,
        length=length,
        gc_percent=round(gc_percent, 3),
        ambiguous_count=ambiguous_count,
        has_core_genes=has_core_genes,
        qc_passed=qc_passed,
        qc_reason=qc_reason_str,
    )


def curate_tier_1() -> list[PlastomeMetadata]:
    """Fetches and QC screens all Tier 1 curated representatives."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    metadata_list = []
    print("\n--- Curating Tier 1 Representatives & Musaceae Outgroups ---")
    for taxon_name, info in TIER_1_REPRESENTATIVES.items():
        gb_file, _ = fetch_accession(info["accession"])
        record = SeqIO.read(gb_file, "genbank")
        meta = evaluate_qc(record, section=info["section"])
        metadata_list.append(meta)
        status_tag = "PASS" if meta.qc_passed else "FAIL"
        print(f"[{status_tag}] {taxon_name:<20} ({meta.accession}) | Length: {meta.length} bp | GC: {meta.gc_percent}% | QC: {meta.qc_reason}")

    # Save metadata summary
    summary_path = DATA_DIR / "tier1_metadata_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump([asdict(m) for m in metadata_list], f, indent=2)

    print(f"\nSaved metadata summary to: {summary_path}")
    return metadata_list


if __name__ == "__main__":
    curate_tier_1()
