# ani_tree.sh

A bash script that compares bacterial genome assemblies by Average Nucleotide Identity (ANI) using pyANI, then builds and visualizes a distance tree from the results. Accepts `.gbk`, `.fna`, `.fa`, or `.fasta` input files — GenBank `.gbk` files are automatically converted to `.fna` before running pyANI.

---

## Outputs

| File | Description |
|---|---|
| `<outdir>/ani/` | pyANI working directory |
| `<outdir>/ani_matrix.tsv` | Full pairwise ANI matrix |
| `<outdir>/ani_tree.nwk` | Newick tree (UPGMA from ANI distance) |
| `<outdir>/ani_tree.pdf` | PDF visualization of the tree and heatmap |

---

## Dependencies

Install Python dependencies via pip:

```bash
pip install pyani scipy numpy matplotlib biopython
```

**Additional dependencies depending on ANI method:**

| Method | Requirement |
|---|---|
| `ANIb` (default) | BLAST+ (`blastn`) — [installation guide](https://github.com/surtlab/Baltrus_Lab_Tutorials/blob/main/Blast-Tutorial/blast-install-guide.md) |
| `ANIm` | MUMmer (`nucmer`) — install from [https://github.com/mummer4/mummer](https://github.com/mummer4/mummer) |

---

## Usage

**Make the script executable (one-time setup):**

```bash
chmod +x ani_tree.sh
```

**Basic usage:**

```bash
./ani_tree.sh -i ./genomes -o ./ani_output
```

**All options:**

```bash
./ani_tree.sh -i <genome_dir> [-o <outdir>] [-t <threads>] [-m <method>]
```

| Flag | Description | Default |
|---|---|---|
| `-i` | Directory containing genome files (required) | — |
| `-o` | Output directory | `./ani_output` |
| `-t` | Number of threads | `4` |
| `-m` | ANI method: `ANIb` (BLAST-based) or `ANIm` (MUMmer-based) | `ANIb` |

---

## Examples

**Run with default settings (ANIb, 4 threads):**

```bash
./ani_tree.sh -i ./genomes -o ./ani_output
```

**Run with MUMmer-based ANI and 8 threads:**

```bash
./ani_tree.sh -i ./genomes -o ./ani_output -t 8 -m ANIm
```

**Run with a mix of .gbk and .fna files:**

```bash
./ani_tree.sh -i ./mixed_genomes -o ./ani_output
```

GenBank `.gbk` files will be automatically converted to `.fna` before the analysis runs.

---

## Input requirements

- The input directory must contain at least **2 genome files** for comparison
- Accepted formats: `.gbk`, `.gbff`, `.fna`, `.fa`, `.fasta`
- All genome files should be placed in the same input directory

---

## Troubleshooting

**"Missing required dependencies"**
Install pyANI and its dependencies: `pip install pyani scipy numpy matplotlib biopython`

**"ANIb requires blastn"**
BLAST+ is not installed or not in your PATH. See the [BLAST installation guide](https://github.com/surtlab/Baltrus_Lab_Tutorials/blob/main/Blast-Tutorial/blast-install-guide.md).

**"Need at least 2 genome files to compare"**
The script found fewer than 2 valid genome files in the input directory. Check that your files are in a supported format and located directly in the input directory (not in subdirectories).

**"Could not find pyANI percentage identity output"**
pyANI may have failed silently. Check that all dependencies are installed correctly and that your genome files are valid FASTA/GenBank format.
