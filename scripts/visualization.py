"""Scientific visualization adhering to Edward Tufte minimalist design and AESTHETICS.md."""

from pathlib import Path
import json
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from pycirclize import Circos

# Ensure workspace root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.config import (
    DATA_DIR,
    OUTPUT_DIR,
    PALETTE,
    apply_publication_theme,
)


def plot_circular_plastome(manifest_path: Path, output_prefix: Path) -> Path:
    """Generates a publication-grade circular plastome map with pycirclize."""
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    total_len = manifest["genome_length"]
    lsc_len = manifest["lsc"]["length"]
    irb_len = manifest["irb"]["length"]
    ssc_len = manifest["ssc"]["length"]
    ira_len = manifest["ira"]["length"]

    sectors = {
        "LSC": lsc_len,
        "IR_B": irb_len,
        "SSC": ssc_len,
        "IR_A": ira_len,
    }

    circos = Circos(sectors, space=1.5)

    # Track 1: Quadripartite Compartments
    region_colors = {
        "LSC": PALETTE["LSC"],
        "IR_B": PALETTE["IR"],
        "SSC": PALETTE["SSC"],
        "IR_A": PALETTE["IR"],
    }

    for sector in circos.sectors:
        # Outer backbone
        track = sector.add_track((92, 100))
        track.axis(fc=region_colors[sector.name], ec="none")
        track.text(sector.name, color="#FFFFFF", size=9, fontweight="bold")

    fig = circos.plotfig()
    
    # Save at 600 DPI and SVG
    png_out = Path(f"{output_prefix}_circular.png")
    svg_out = Path(f"{output_prefix}_circular.svg")

    fig.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="#FFFFFF")
    fig.savefig(svg_out, format="svg", bbox_inches="tight", facecolor="#FFFFFF")
    plt.close(fig)

    print(f"[Visualization] Saved circular plastome map to: {png_out}")
    return png_out


def plot_divergence_landscape(metrics_df: pd.DataFrame, manifest_path: Path, output_prefix: Path) -> Path:
    """Plots sliding-window nucleotide diversity (pi), Shannon entropy (H), and GC skew."""
    apply_publication_theme()

    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    lsc_end = manifest["lsc"]["length"]
    irb_end = lsc_end + manifest["irb"]["length"]
    ssc_end = irb_end + manifest["ssc"]["length"]
    total_len = manifest["genome_length"]

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9.5, 6.5), sharex=True, dpi=300)

    x = metrics_df["midpoint"]
    pi = metrics_df["pi"]
    entropy = metrics_df["entropy"]
    gc_skew = metrics_df["gc_skew"]

    # Boundary lines and compartment shading
    for ax in [ax1, ax2, ax3]:
        ax.axvline(lsc_end, linestyle="--", linewidth=0.8, color="#888888", zorder=0)
        ax.axvline(irb_end, linestyle="--", linewidth=0.8, color="#888888", zorder=0)
        ax.axvline(ssc_end, linestyle="--", linewidth=0.8, color="#888888", zorder=0)

    # Panel a: Nucleotide Diversity (pi)
    ax1.plot(x, pi, color=PALETTE["BACKBONE"], linewidth=1.1, label=r"Nucleotide Diversity ($\pi$)")
    mean_pi = pi.mean()
    ax1.axhline(mean_pi, linestyle=":", color="#777777", linewidth=0.8, label=rf"Mean $\pi$ ({mean_pi:.4f})")
    ax1.set_ylabel(r"Diversity ($\pi$)")
    ax1.text(-0.06, 1.05, "a", transform=ax1.transAxes, fontsize=11, fontweight="bold", va="top")
    ax1.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    # Panel b: Shannon Entropy (H)
    ax2.plot(x, entropy, color=PALETTE["LSC"], linewidth=1.0, label="Shannon Entropy ($H$)")
    ax2.set_ylabel(r"Entropy ($H$)")
    ax2.text(-0.06, 1.05, "b", transform=ax2.transAxes, fontsize=11, fontweight="bold", va="top")
    ax2.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    # Panel c: GC Skew
    ax3.plot(x, gc_skew, color=PALETTE["SSC"], linewidth=0.9, label="GC Skew")
    ax3.axhline(0, linestyle="-", color="#AAAAAA", linewidth=0.6)
    ax3.set_ylabel("GC Skew")
    ax3.set_xlabel("Ancestral Coordinate (bp)")
    ax3.text(-0.06, 1.05, "c", transform=ax3.transAxes, fontsize=11, fontweight="bold", va="top")
    ax3.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    # Label Quadripartite regions at the top
    y_top = 1.15
    ax1.text(lsc_end / 2, ax1.get_ylim()[1] * 1.05, "LSC", ha="center", fontsize=9, fontweight="bold", color=PALETTE["LSC"])
    ax1.text((lsc_end + irb_end) / 2, ax1.get_ylim()[1] * 1.05, r"$\mathrm{IR_B}$", ha="center", fontsize=9, fontweight="bold", color=PALETTE["IR"])
    ax1.text((irb_end + ssc_end) / 2, ax1.get_ylim()[1] * 1.05, "SSC", ha="center", fontsize=9, fontweight="bold", color=PALETTE["SSC"])
    ax1.text((ssc_end + total_len) / 2, ax1.get_ylim()[1] * 1.05, r"$\mathrm{IR_A}$", ha="center", fontsize=9, fontweight="bold", color=PALETTE["IR"])

    plt.subplots_adjust(right=0.78, hspace=0.25)

    png_out = Path(f"{output_prefix}_landscape.png")
    svg_out = Path(f"{output_prefix}_landscape.svg")

    fig.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="#FFFFFF")
    fig.savefig(svg_out, format="svg", bbox_inches="tight", facecolor="#FFFFFF")
    plt.close(fig)

    print(f"[Visualization] Saved divergence landscape to: {png_out}")
    return png_out


def plot_junction_shifts(junctions_path: Path, output_prefix: Path) -> Path:
    """Generates a publication-grade IR junction shift diagram across Musa sections and outgroups."""
    apply_publication_theme()

    with open(junctions_path, "r", encoding="utf-8") as f:
        juncs = json.load(f)

    fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)

    y_labels = []
    y_positions = []

    for idx, item in enumerate(reversed(juncs)):
        y = idx
        y_positions.append(y)
        acc = item["accession"].split(".")[0]
        y_labels.append(acc)

        lsc = item["lsc_len"]
        irb = item["irb_len"]
        ssc = item["ssc_len"]
        ira = item["ira_len"]

        # Stacked horizontal bar
        ax.barh(y, lsc, left=0, color=PALETTE["LSC"], height=0.55, edgecolor="none")
        ax.barh(y, irb, left=lsc, color=PALETTE["IR"], height=0.55, edgecolor="none")
        ax.barh(y, ssc, left=lsc + irb, color=PALETTE["SSC"], height=0.55, edgecolor="none")
        ax.barh(y, ira, left=lsc + irb + ssc, color=PALETTE["IR"], height=0.55, edgecolor="none")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel("Plastome Length (bp)")

    # Legend
    legend_patches = [
        Patch(facecolor=PALETTE["LSC"], label="LSC"),
        Patch(facecolor=PALETTE["IR"], label=r"$\mathrm{IR_B} / \mathrm{IR_A}$"),
        Patch(facecolor=PALETTE["SSC"], label="SSC"),
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    png_out = Path(f"{output_prefix}_junction_shifts.png")
    svg_out = Path(f"{output_prefix}_junction_shifts.svg")

    fig.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="#FFFFFF")
    fig.savefig(svg_out, format="svg", bbox_inches="tight", facecolor="#FFFFFF")
    plt.close(fig)

    print(f"[Visualization] Saved junction shifts plot to: {png_out}")
    return png_out


if __name__ == "__main__":
    pass
