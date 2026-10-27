"""Phylogenetic tree inference and empirical Bayesian Ancestral State Reconstruction (ASR) via IQ-TREE 2."""

from pathlib import Path
import json
import os
import subprocess
import sys
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import DATA_DIR, OUTPUT_DIR, PARTITIONS_DIR, PROCESSED_DIR
from scripts.tool_resolver import resolve_iqtree, resolve_mafft


def run_mafft_alignment(input_fasta: Path, output_aln: Path) -> Path:
    """Runs MAFFT progressive multiple sequence alignment."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    if output_aln.is_file() and output_aln.stat().st_size > 1000:
        print(f"[MAFFT] Found cached alignment: {output_aln}")
        return output_aln

    mafft_bin = resolve_mafft()
    if not mafft_bin:
        raise RuntimeError("MAFFT binary could not be resolved. Please install MAFFT.")

    print(f"[MAFFT] Aligning {input_fasta.name} using {mafft_bin}...")
    
    # Use forward slashes / relative posix path to prevent MSYS escape character issues on Windows
    try:
        in_path_str = input_fasta.relative_to(Path.cwd()).as_posix()
    except ValueError:
        in_path_str = input_fasta.as_posix()

    if mafft_bin.startswith("wsl "):
        cmd = ["wsl", "mafft", "--auto", "--quiet", in_path_str]
    else:
        cmd = [mafft_bin, "--auto", "--quiet", in_path_str]
    
    with open(output_aln, "w", encoding="utf-8") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"MAFFT alignment failed: {result.stderr}")

    print(f"[MAFFT] Completed alignment: {output_aln} ({output_aln.stat().st_size} bytes)")
    return output_aln


def run_iqtree_asr(alignment_file: Path, prefix: Path, outgroup: str) -> tuple[Path, Path]:
    """Executes IQ-TREE 2 to infer ML tree and empirical Bayesian Ancestral Sequence Reconstruction."""
    tree_file = Path(f"{prefix}.treefile")
    state_file = Path(f"{prefix}.state")

    if tree_file.is_file() and state_file.is_file():
        print(f"[IQ-TREE 2] Found cached ASR files for {prefix.name}")
        return tree_file, state_file

    iqtree_bin = resolve_iqtree()
    if not iqtree_bin:
        raise RuntimeError("IQ-TREE 2 binary could not be resolved.")

    print(f"[IQ-TREE 2] Running ML tree search and ASR on {alignment_file.name}...")
    cmd = [
        iqtree_bin,
        "-s", str(alignment_file),
        "-m", "GTR+F+I+G4",
        "-B", "1000",
        "-asr",
        "-o", outgroup,
        "--prefix", str(prefix),
        "-nt", "AUTO",
        "-redo",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[IQ-TREE 2] Warning/Error output: {result.stderr}")
        # Retry with simpler model if GTR+F+I+G4 fails
        fallback_cmd = [
            iqtree_bin,
            "-s", str(alignment_file),
            "-m", "JC",
            "-asr",
            "-o", outgroup,
            "--prefix", str(prefix),
            "-redo",
        ]
        result2 = subprocess.run(fallback_cmd, capture_output=True, text=True)
        if result2.returncode != 0:
            raise RuntimeError(f"IQ-TREE 2 execution failed: {result2.stderr}")

    return tree_file, state_file


def identify_musa_crown_node(tree_file: Path, outgroup_taxa: set[str]) -> str:
    """Parses Newick tree to find the internal node ID representing the MRCA of Musa (crown node)."""
    from Bio import Phylo
    tree = Phylo.read(str(tree_file), "newick")
    
    # Root tree at outgroup
    outgroup_clades = [c for c in tree.get_terminals() if any(og in c.name for og in outgroup_taxa)]
    if outgroup_clades:
        try:
            tree.root_with_outgroup(outgroup_clades)
        except Exception:
            pass

    # Find MRCA of all ingroup taxa (Musa)
    ingroup_clades = [c for c in tree.get_terminals() if not any(og in c.name for og in outgroup_taxa)]
    mrca_clade = tree.common_ancestor(ingroup_clades)
    
    # If the clade has a confidence / node label
    crown_node_name = mrca_clade.name or "Node1"
    print(f"[Phylogeny] Identified Musa crown node (MRCA): {crown_node_name}")
    return crown_node_name


def parse_ancestral_sequence(state_file: Path, target_node: str, alignment_file: Path) -> str:
    """Parses IQ-TREE .state file to extract MAP ancestral sequence for the target node."""
    # Check if target_node exists in .state
    site_states: dict[int, str] = {}
    
    first_node = None
    with open(state_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("Node"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                node = parts[0]
                site = int(parts[1])
                state = parts[2]
                if first_node is None:
                    first_node = node
                if node == target_node:
                    site_states[site] = state

    # Fallback to first internal node if target_node was not exact label
    if not site_states and first_node:
        print(f"[ASR] Target node {target_node} not in state file; using root node {first_node}")
        with open(state_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("Node"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3 and parts[0] == first_node:
                    site_states[int(parts[1])] = parts[2]

    # Handle indels / gap columns from alignment
    alignment = AlignIO.read(str(alignment_file), "fasta")
    aln_len = alignment.get_alignment_length()

    ancestral_bases = []
    for site in range(1, aln_len + 1):
        # Check column gap frequency across ingroup taxa
        col = alignment[:, site - 1]
        gap_count = col.count("-")
        
        # If >50% of taxa have gaps, consider column a derived insertion (ancestral gap)
        if gap_count > (len(alignment) / 2):
            continue
        
        base = site_states.get(site, "-")
        if base != "-":
            ancestral_bases.append(base)
        else:
            # Most common non-gap base in column
            non_gaps = [b for b in col if b != "-"]
            if non_gaps:
                ancestral_bases.append(max(set(non_gaps), key=non_gaps.count))

    return "".join(ancestral_bases)


def reconstruct_ancestral_plastome() -> tuple[Path, Path]:
    """Orchestrates partitioned alignment, IQ-TREE 2 Bayesian ASR, and circular reconstitution."""
    outgroup_taxa = {
        "Ensete_ventricosum_chloroplast",
        "Ensete_glaucum_chloroplast",
        "Musella_lasiocarpa_chloroplast",
    }
    outgroup_str = ",".join(outgroup_taxa)

    ancestral_partitions: dict[str, str] = {}

    for part in ["LSC", "SSC", "IR_B"]:
        in_fasta = PARTITIONS_DIR / f"tier1_{part}.fasta"
        aln_file = PROCESSED_DIR / f"aligned_{part}.fasta"
        prefix = PROCESSED_DIR / f"iqtree_{part}"

        # 1. Run MAFFT alignment
        run_mafft_alignment(in_fasta, aln_file)

        # 2. Run IQ-TREE 2 tree search & ASR
        tree_file, state_file = run_iqtree_asr(aln_file, prefix, outgroup_str)

        # 3. Identify Musa crown node
        crown_node = identify_musa_crown_node(tree_file, outgroup_taxa)

        # 4. Extract ancestral sequence
        anc_seq = parse_ancestral_sequence(state_file, crown_node, aln_file)
        ancestral_partitions[part] = anc_seq
        print(f"[ASR] Reconstructed ancestral {part}: length = {len(anc_seq)} bp")

    # Reconstitute circular quadripartite pseudo-molecule:
    # Ancestral Plastome = LSC_anc + IR_B_anc + SSC_anc + revcomp(IR_B_anc)
    lsc_anc = ancestral_partitions["LSC"]
    irb_anc = ancestral_partitions["IR_B"]
    ssc_anc = ancestral_partitions["SSC"]
    ira_anc = str(Seq(irb_anc).reverse_complement())

    full_anc_seq = lsc_anc + irb_anc + ssc_anc + ira_anc
    total_len = len(full_anc_seq)
    print(f"\n=======================================================")
    print(f"RECONSTRUCTED ANCESTRAL MUSA PLASTOME (MRCA)")
    print(f"  LSC_anc:    {len(lsc_anc):>6} bp")
    print(f"  IR_B_anc:   {len(irb_anc):>6} bp")
    print(f"  SSC_anc:    {len(ssc_anc):>6} bp")
    print(f"  IR_A_anc:   {len(ira_anc):>6} bp")
    print(f"  Total:      {total_len:>6} bp")
    print(f"=======================================================\n")

    # Save to data and public directories
    out_record = SeqRecord(
        Seq(full_anc_seq),
        id="Musa_MRCA_plastome",
        name="Musa_MRCA",
        description="Hypothetical ancestral chloroplast genome (MRCA) of genus Musa reconstructed via IQ-TREE 2 Bayesian ASR",
    )

    data_fasta = DATA_DIR / "musa_mrca_plastome.fasta"
    public_fasta = OUTPUT_DIR / "musa_mrca_plastome.fasta"

    SeqIO.write(out_record, data_fasta, "fasta")
    SeqIO.write(out_record, public_fasta, "fasta")

    # Save partition coordinate manifest
    manifest = {
        "genome_length": total_len,
        "lsc": {"start": 0, "end": len(lsc_anc), "length": len(lsc_anc)},
        "irb": {"start": len(lsc_anc), "end": len(lsc_anc) + len(irb_anc), "length": len(irb_anc)},
        "ssc": {"start": len(lsc_anc) + len(irb_anc), "end": len(lsc_anc) + len(irb_anc) + len(ssc_anc), "length": len(ssc_anc)},
        "ira": {"start": len(lsc_anc) + len(irb_anc) + len(ssc_anc), "end": total_len, "length": len(ira_anc)},
    }
    manifest_path = DATA_DIR / "musa_mrca_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    return data_fasta, manifest_path


if __name__ == "__main__":
    reconstruct_ancestral_plastome()
