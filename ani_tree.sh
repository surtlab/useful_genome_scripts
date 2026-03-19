#!/usr/bin/env bash
# =============================================================================
# ani_tree.sh
# Compare bacterial genome assemblies by ANI using pyANI, then build and
# visualise a distance tree from the results.
#
# Accepts .gbk, .fna, .fa, or .fasta input files — .gbk files are
# automatically converted to .fna before running pyANI.
#
# Outputs:
#   - <outdir>/ani/                  pyANI working directory
#   - <outdir>/ani_matrix.tsv        full pairwise ANI matrix
#   - <outdir>/ani_tree.nwk          Newick tree (UPGMA from ANI distance)
#   - <outdir>/ani_tree.pdf          PDF visualisation of the tree
#
# Install dependencies:
#   pip install pyani scipy numpy matplotlib
#
# Usage:
#   chmod +x ani_tree.sh
#   ./ani_tree.sh -i ./genomes -o ./ani_output
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
INDIR=""
OUTDIR="./ani_output"
THREADS=4
METHOD="ANIb"   # ANIb (BLAST-based) or ANIm (MUMmer-based)

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'

usage() {
    echo "Usage: $0 -i <genome_dir> [-o <outdir>] [-t <threads>] [-m <method>]"
    echo ""
    echo "  -i  Directory containing genome files (.gbk, .fna, .fa, .fasta) (required)"
    echo "  -o  Output directory                                              (default: ./ani_output)"
    echo "  -t  Number of threads                                             (default: 4)"
    echo "  -m  ANI method: ANIb (BLAST) or ANIm (MUMmer)                    (default: ANIb)"
    echo ""
    echo "Example:"
    echo "  $0 -i ./genomes -o ./ani_output -t 8"
    exit 1
}

while getopts ":i:o:t:m:h" opt; do
    case $opt in
        i) INDIR="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: Option -$OPTARG requires an argument."; usage ;;
        \?) echo "ERROR: Unknown option -$OPTARG"; usage ;;
    esac
done

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ -z "$INDIR" ]]; then
    echo -e "${RED}ERROR:${NC} No input directory specified. Use -i ./genomes"
    usage
fi
if [[ ! -d "$INDIR" ]]; then
    echo -e "${RED}ERROR:${NC} Input directory not found: $INDIR"
    exit 1
fi

case "$METHOD" in
    ANIb|ANIm) ;;
    *) echo -e "${RED}ERROR:${NC} Invalid method '$METHOD'. Use ANIb or ANIm."; exit 1 ;;
esac

# ── Dependency check ──────────────────────────────────────────────────────────
missing=()
for cmd in python3 average_nucleotide_identity.py; do
    command -v "$cmd" &>/dev/null || missing+=("$cmd")
done

if [[ ${#missing[@]} -gt 0 ]]; then
    echo -e "${RED}ERROR:${NC} Missing required dependencies: ${missing[*]}"
    echo "       Install with: pip install pyani"
    echo "       Also required: pip install scipy numpy matplotlib biopython"
    exit 1
fi

if [[ "$METHOD" == "ANIb" ]] && ! command -v blastn &>/dev/null; then
    echo -e "${RED}ERROR:${NC} ANIb requires blastn. Install BLAST+ from:"
    echo "       https://www.ncbi.nlm.nih.gov/books/NBK569861/"
    exit 1
fi

if [[ "$METHOD" == "ANIm" ]] && ! command -v nucmer &>/dev/null; then
    echo -e "${RED}ERROR:${NC} ANIm requires MUMmer (nucmer). Install from:"
    echo "       https://github.com/mummer4/mummer"
    exit 1
fi

# Check required Python packages
python3 -c "import scipy, numpy, matplotlib" 2>/dev/null || {
    echo -e "${RED}ERROR:${NC} Missing required Python packages."
    echo "       Install with: pip install scipy numpy matplotlib"
    exit 1
}

mkdir -p "$OUTDIR"
FNADIR="${OUTDIR}/fna"
ANIDIR="${OUTDIR}/ani"
mkdir -p "$FNADIR"
mkdir -p "$ANIDIR"

# ── Convert .gbk → .fna if needed ────────────────────────────────────────────
convert_gbk() {
    local gbk="$1"
    local base
    base=$(basename "$gbk" | sed 's/\.[^.]*$//')
    local out="${FNADIR}/${base}.fna"

    python3 << PYEOF
import re

gbk_file = "$gbk"
out_file = "$out"

records = []
current = []
with open(gbk_file) as f:
    for line in f:
        current.append(line)
        if line.strip() == "//":
            records.append("".join(current))
            current = []

with open(out_file, "w") as out:
    for rec in records:
        acc  = ""
        defn = ""
        for line in rec.splitlines():
            if line.startswith("ACCESSION"):
                acc = line.split()[1]
            if line.startswith("DEFINITION"):
                defn = line[12:].strip()
        if not acc:
            continue
        in_origin = False
        nuc_seq = []
        for line in rec.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if line.strip() == "//":
                break
            if in_origin:
                nuc_seq.append(re.sub(r'[^a-zA-Z]', '', line))
        nuc_str = "".join(nuc_seq).upper()
        if nuc_str:
            out.write(f">{acc} {defn}\n")
            for i in range(0, len(nuc_str), 60):
                out.write(nuc_str[i:i+60] + "\n")
PYEOF
    echo "$out"
}

# ── Collect genome files ──────────────────────────────────────────────────────
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  ANI Tree Builder"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Input dir  : $INDIR"
echo "  Output dir : $OUTDIR"
echo "  Method     : $METHOD"
echo "  Threads    : $THREADS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo -e "${CYAN}[PREP]${NC} Collecting and converting genome files..."

gbk_count=0
fna_count=0

for f in "${INDIR}"/*.gbk "${INDIR}"/*.gbff; do
    [[ -f "$f" ]] || continue
    base=$(basename "$f" | sed 's/\.[^.]*$//')
    echo -e "  ${YELLOW}[GBK]${NC}  Converting ${base}..."
    convert_gbk "$f"
    ((gbk_count++)) || true
done

for f in "${INDIR}"/*.fna "${INDIR}"/*.fa "${INDIR}"/*.fasta; do
    [[ -f "$f" ]] || continue
    base=$(basename "$f" | sed 's/\.[^.]*$//')
    cp "$f" "${FNADIR}/${base}.fna"
    ((fna_count++)) || true
done

total=$(find "$FNADIR" -name "*.fna" | wc -l | tr -d ' ')

if [[ "$total" -lt 2 ]]; then
    echo -e "${RED}ERROR:${NC} Need at least 2 genome files to compare. Found: ${total}"
    exit 1
fi

echo -e "  ${GREEN}[OK]${NC}   ${total} genomes ready (${gbk_count} converted from .gbk, ${fna_count} already .fna)"
echo ""

# ── Run pyANI ─────────────────────────────────────────────────────────────────
echo -e "${CYAN}[ANI]${NC}  Running pyANI ${METHOD} on ${total} genomes..."

average_nucleotide_identity.py \
    -i "$FNADIR" \
    -o "$ANIDIR" \
    -m "$METHOD" \
    --workers "$THREADS" \
    --force \
    2>/dev/null

# Locate the ANI percentage matrix output from pyANI
ANI_MATRIX_CSV=""
for f in "${ANIDIR}"/*percentage_identity*.tab \
         "${ANIDIR}"/*_percentage_identity*.csv \
         "${ANIDIR}"/*ANIb_percentage_identity* \
         "${ANIDIR}"/*ANIm_percentage_identity*; do
    [[ -f "$f" ]] && ANI_MATRIX_CSV="$f" && break
done

if [[ -z "$ANI_MATRIX_CSV" ]]; then
    # Fall back: find any .tab file in the output
    ANI_MATRIX_CSV=$(find "$ANIDIR" -name "*.tab" | grep -i "identity" | head -1 || true)
fi

if [[ -z "$ANI_MATRIX_CSV" || ! -f "$ANI_MATRIX_CSV" ]]; then
    echo -e "${RED}ERROR:${NC} Could not find pyANI percentage identity output in $ANIDIR"
    echo "       Files present:"
    ls "$ANIDIR" | sed 's/^/         /'
    exit 1
fi

echo -e "  ${GREEN}[OK]${NC}   pyANI complete → ${ANI_MATRIX_CSV}"
echo ""

# ── Build tree and PDF from pyANI matrix ─────────────────────────────────────
echo -e "${CYAN}[TREE]${NC} Building distance matrix and tree..."

python3 << PYEOF
import os, sys, re
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ani_csv     = "$ANI_MATRIX_CSV"
matrix_file = "$OUTDIR/ani_matrix.tsv"
nwk_file    = "$OUTDIR/ani_tree.nwk"
pdf_file    = "$OUTDIR/ani_tree.pdf"

# ── Parse pyANI matrix (tab-delimited, first row and col are genome names) ────
genomes = []
rows    = []

with open(ani_csv) as f:
    for i, line in enumerate(f):
        cols = line.rstrip("\n").split("\t")
        if i == 0:
            # Header row — genome names (skip first empty cell)
            genomes = [os.path.basename(c).replace(".fna","") for c in cols[1:]]
        else:
            label = os.path.basename(cols[0]).replace(".fna","")
            vals  = []
            for v in cols[1:]:
                try:
                    vals.append(float(v) * 100 if float(v) <= 1.0 else float(v))
                except ValueError:
                    vals.append(0.0)
            rows.append(vals)

n            = len(genomes)
ani_matrix   = np.array(rows)

# Ensure self-comparisons are 100
np.fill_diagonal(ani_matrix, 100.0)

# Symmetrise
ani_matrix = (ani_matrix + ani_matrix.T) / 2
np.fill_diagonal(ani_matrix, 100.0)

# ── Write ANI matrix ──────────────────────────────────────────────────────────
with open(matrix_file, "w") as out:
    out.write("\t" + "\t".join(genomes) + "\n")
    for i, g in enumerate(genomes):
        row = "\t".join(f"{ani_matrix[i][j]:.4f}" for j in range(n))
        out.write(f"{g}\t{row}\n")

print(f"  ANI matrix written → {matrix_file}")

# ── UPGMA tree from ANI distance ──────────────────────────────────────────────
dist_matrix = 100.0 - ani_matrix
np.fill_diagonal(dist_matrix, 0.0)
dist_matrix = np.clip(dist_matrix, 0, None)
condensed       = squareform(dist_matrix)
linkage_matrix  = linkage(condensed, method="average")

# ── Newick ────────────────────────────────────────────────────────────────────
def linkage_to_newick(Z, labels):
    n = len(labels)
    nodes = {i: (labels[i], 0.0) for i in range(n)}
    for i, (left, right, dist, _) in enumerate(Z):
        left, right = int(left), int(right)
        left_label,  left_acc  = nodes[left]
        right_label, right_acc = nodes[right]
        lb = max(dist / 2 - left_acc,  0.0)
        rb = max(dist / 2 - right_acc, 0.0)
        newick = f"({left_label}:{lb:.6f},{right_label}:{rb:.6f})"
        nodes[n + i] = (newick, dist / 2)
    return nodes[n + len(Z) - 1][0] + ";"

newick = linkage_to_newick(linkage_matrix, genomes)
with open(nwk_file, "w") as out:
    out.write(newick + "\n")
print(f"  Newick tree written → {nwk_file}")

# ── PDF: dendrogram + heatmap ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, max(6, n * 0.6 + 2)),
                         gridspec_kw={"width_ratios": [2, 3]})
fig.suptitle("Pairwise ANI Tree", fontsize=14, fontweight="bold")

ax_tree = axes[0]
dend = dendrogram(
    linkage_matrix,
    labels=genomes,
    orientation="left",
    ax=ax_tree,
    color_threshold=0,
    above_threshold_color="steelblue",
    leaf_font_size=9,
)
ax_tree.set_xlabel("ANI distance (100 - ANI %)", fontsize=9)
ax_tree.set_title("UPGMA tree", fontsize=10)
ax_tree.spines["top"].set_visible(False)
ax_tree.spines["right"].set_visible(False)

order          = dend["leaves"]
ordered_genomes = [genomes[i] for i in order]
ordered_matrix  = ani_matrix[np.ix_(order, order)]

ax_heat = axes[1]
vmin = max(80, float(np.min(ani_matrix[ani_matrix > 0])) - 1)
im = ax_heat.imshow(ordered_matrix, aspect="auto", cmap="RdYlGn",
                    vmin=vmin, vmax=100)
ax_heat.set_xticks(range(n))
ax_heat.set_xticklabels(ordered_genomes, rotation=45, ha="right", fontsize=8)
ax_heat.set_yticks(range(n))
ax_heat.set_yticklabels(ordered_genomes, fontsize=8)
ax_heat.set_title("ANI heatmap (%)", fontsize=10)

for i in range(n):
    for j in range(n):
        val   = ordered_matrix[i][j]
        color = "white" if val < 92 else "black"
        ax_heat.text(j, i, f"{val:.1f}", ha="center", va="center",
                     fontsize=6, color=color)

plt.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.04, label="ANI %")
plt.tight_layout()
plt.savefig(pdf_file, bbox_inches="tight", dpi=150)
print(f"  Tree PDF written    → {pdf_file}")
PYEOF

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Genomes compared : $total"
echo "  ANI matrix       : ${OUTDIR}/ani_matrix.tsv"
echo "  Newick tree      : ${OUTDIR}/ani_tree.nwk"
echo "  Tree PDF         : ${OUTDIR}/ani_tree.pdf"
echo "  pyANI output     : ${ANIDIR}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
