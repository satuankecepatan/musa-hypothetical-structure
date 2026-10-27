"""Global configuration, path resolution, and aesthetic defaults for Musa plastome analysis."""

from pathlib import Path
import sys
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Environment & Path Resolution
# -----------------------------------------------------------------------------
IS_COLAB: bool = "google.colab" in sys.modules

if IS_COLAB:
    BASE_DIR = Path("/content/musa-hypothetical-structure")
else:
    BASE_DIR = Path(__file__).resolve().parent.parent

DATA_DIR = BASE_DIR / "data"
RAW_GB_DIR = DATA_DIR / "raw_genbank"
RAW_FASTA_DIR = DATA_DIR / "raw_fasta"
PROCESSED_DIR = DATA_DIR / "processed"
PARTITIONS_DIR = DATA_DIR / "partitions"
OUTPUT_DIR = BASE_DIR / "public"
BIN_DIR = BASE_DIR / "bin"
NOTEBOOKS_DIR = BASE_DIR / "notebooks"

for directory in [
    DATA_DIR,
    RAW_GB_DIR,
    RAW_FASTA_DIR,
    PROCESSED_DIR,
    PARTITIONS_DIR,
    OUTPUT_DIR,
    BIN_DIR,
    NOTEBOOKS_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Entrez Email Configuration
# -----------------------------------------------------------------------------
ENTREZ_EMAIL: str = "antigravity.musa.research@example.com"
ENTREZ_TOOL: str = "MusaPlastomeReconstruction"

# -----------------------------------------------------------------------------
# Publication Aesthetic Configuration (per AESTHETICS.md)
# -----------------------------------------------------------------------------
PALETTE: dict[str, str] = {
    # Quadripartite compartments
    "LSC": "#1B365D",        # Deep Slate Navy
    "SSC": "#E05A47",        # Vibrant Coral / Vermilion
    "IR": "#D99B26",         # Muted Amber / Ochre
    "IR_A": "#D99B26",
    "IR_B": "#D99B26",
    "BACKBONE": "#222222",   # Charcoal Black
    "MUTED": "#B0B5B9",      # Cool Stone Grey
    # Taxonomic section accents
    "Eumusa": "#008080",        # Electric Teal
    "Rhodochlamys": "#B8336A",  # Crimson Magenta
    "Callimusa": "#2E5EAA",     # Ultramarine Blue
    "Australimusa": "#D16014",  # Burnt Sienna
    "Outgroup": "#5A6B7C",      # Slate Outgroup Grey
}


def apply_publication_theme() -> None:
    """Configures matplotlib rcParams globally to enforce high data-ink ratio principles."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 11,
        "lines.linewidth": 1.2,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "figure.facecolor": "#FFFFFF",
        "axes.facecolor": "#FFFFFF",
        "savefig.facecolor": "#FFFFFF",
        "savefig.edgecolor": "#FFFFFF",
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })
