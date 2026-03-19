# pull_gene_nuc_seq_from_gbk.py

Extract nucleotide sequences from one or more GenBank files, either by gene annotation or by genomic coordinates derived from a BLAST search.

---

## What's new in v1.3

### Bug fix: duplicate sequences in gene mode
GenBank files often contain both a `gene` feature and a `CDS` feature at identical coordinates for the same locus. Previous versions returned both, producing duplicate sequences in the output. Gene mode now deduplicates by `(record_id, start, end, strand)` — identical coordinates are only extracted once regardless of feature type. Use `--feature-type CDS` if you specifically want only CDS entries.

### Bug fix: duplicate sequences in BLAST mode
BLAST output files can contain multiple hits that map to the same genomic region — for example when the same contig is hit by more than one query sequence, or when overlapping hits pass the filters. BLAST mode now deduplicates by `(sseqid, start, end, strand)` after applying padding, so each unique region is only extracted once. A `[SKIP]` message is printed when a duplicate is detected.

### Bug fix: annotation-only GenBank files (no sequence data)
Some NCBI GenBank files contain annotation but no nucleotide sequence (the `ORIGIN` block is absent). Previous versions crashed with a `Bio.Seq.UndefinedSequenceError`. The script now detects these records upfront, prints a warning, and skips them gracefully rather than exiting.

### Bug fix: BLAST sseqid version suffix mismatch
BLAST databases built from FASTA files sometimes strip the version suffix from sequence IDs (e.g. `NZ_JABEBQ010000003` instead of `NZ_JABEBQ010000003.1`). This caused BLAST mode to silently find no hits when the GenBank records used versioned IDs. The script now builds a secondary lookup with version suffixes stripped, so both `NZ_JABEBQ010000003` and `NZ_JABEBQ010000003.1` resolve to the correct record.

### Feature: folder input
`-i` now accepts either a single GenBank file or a directory. When a directory is given, all `.gbk`/`.gb`/`.genbank` files inside are processed and per-genome output is enabled automatically.

---

## Requirements

- Python 3.7+
- [Biopython](https://biopython.org/)

```bash
pip install biopython
```

---

## Usage

```
python pull_gene_nuc_seq_from_gbk.py -i <input> (-g <gene> | --blast <file>) [options]
```

Exactly one of `-g` / `--blast` must be supplied. They cannot be used together.

---

## Modes

### Gene mode (`-g`)

Searches GenBank features for a matching gene name and extracts the annotated nucleotide sequence. The query is matched case-insensitively against three qualifiers: `/gene`, `/locus_tag`, and `/product`.

```bash
# Extract all features matching "rpoB" from a single genome
python pull_gene_nuc_seq_from_gbk.py -i genome.gbk -g rpoB

# Restrict to CDS features only
python pull_gene_nuc_seq_from_gbk.py -i genome.gbk -g rpoB --feature-type CDS

# Run across a folder of genomes
python pull_gene_nuc_seq_from_gbk.py -i genomes/ -g rpoB
```

### BLAST region mode (`--blast`)

Takes a BLAST `-outfmt 6` tabular file and extracts the subject (genome) region for each hit directly from the GenBank sequence. Reverse-strand hits (`sstart > send`) are automatically reverse-complemented. Duplicate regions are skipped with a `[SKIP]` message.

```bash
# Extract hit regions from a single genome
python pull_gene_nuc_seq_from_gbk.py -i genome.gbk --blast hits.blast6

# Add 200 bp of flanking sequence on each side
python pull_gene_nuc_seq_from_gbk.py -i genome.gbk --blast hits.blast6 --padding 200

# Filter hits and run across a folder of genomes
python pull_gene_nuc_seq_from_gbk.py -i genomes/ --blast hits.blast6 \
    --min-pident 95 --min-length 200 --padding 100 --outdir rpoB_regions/
```

> **Note on sseqid matching:** The `sseqid` values in the BLAST output (column 2) must match contig record IDs in the GenBank file(s). The script handles both versioned (`NZ_JABEBQ010000003.1`) and unversioned (`NZ_JABEBQ010000003`) IDs automatically.

---

## Input

| Flag | Description |
|------|-------------|
| `-i` / `--input` | A single `.gbk` / `.gb` file, **or** a directory containing multiple GenBank files. Recognised extensions: `.gbk`, `.gb`, `.genbank` |

---

## Output

By default, all hits are written to a single pooled FASTA file. When `-i` points to a directory, **per-genome mode is enabled automatically** and one `.fna` file is written per input GenBank.

| Flag | Description |
|------|-------------|
| `-o` / `--outfile` | Output FASTA filename (pooled mode only). Default: `<gene>.fna` or `blast_regions.fna` |
| `--per-genome` | Force per-genome output even when `-i` is a single file |
| `--outdir` | Directory for per-genome output files. Default: `extracted_seqs/` |

### FASTA header format

**Gene mode:**
```
>RECORD_ID|gene_val|locus_tag  product=<product> location=<start>..<end>(<strand>) type=<feature_type> source=<filename>
```

**BLAST mode:**
```
>sseqid|start..end(strand)|query=qseqid  pident=<pident> aln_len=<len> padding=<bp> source=<filename>
```

---

## All options

| Flag | Default | Description |
|------|---------|-------------|
| `-i` / `--input` | *(required)* | Input GenBank file or directory |
| `-g` / `--gene` | — | Gene name to search for (gene mode) |
| `--feature-type` | any | Restrict gene search to a specific feature type (e.g. `CDS`, `rRNA`) |
| `--blast` | — | BLAST outfmt 6 file (BLAST mode) |
| `--padding` | `0` | Bases to add on each side of a BLAST hit |
| `--min-pident` | `0` | Minimum percent identity filter for BLAST hits |
| `--min-length` | `0` | Minimum alignment length filter for BLAST hits |
| `-o` / `--outfile` | auto | Output FASTA path (pooled mode) |
| `--per-genome` | off | Write one FASTA per input genome |
| `--outdir` | `extracted_seqs/` | Output directory for per-genome FASTAs |

---

## Examples

```bash
# Single genome, gene search, CDS only
python pull_gene_nuc_seq_from_gbk.py \
    -i GCF_000001405.gbk \
    -g recA \
    --feature-type CDS \
    -o recA_hits.fna

# Folder of genomes, gene search, one output file per genome
python pull_gene_nuc_seq_from_gbk.py \
    -i genomes/ \
    -g 16S \
    --feature-type rRNA \
    --outdir 16S_seqs/

# Single genome, BLAST hits with flanking sequence and quality filters
python pull_gene_nuc_seq_from_gbk.py \
    -i genome.gbk \
    --blast blastn_results.txt \
    --min-pident 90 \
    --min-length 300 \
    --padding 500 \
    -o flanked_regions.fna

# Folder of genomes, BLAST mode — hits distributed across genomes automatically
python pull_gene_nuc_seq_from_gbk.py \
    -i genomes/ \
    --blast blastn_results.txt \
    --min-pident 95 \
    --outdir blast_regions/
```

---

## Troubleshooting

**"No features found matching gene name"**
Check spelling and capitalisation. Try without `--feature-type` first to see if the gene exists under a different feature type. Gene names in GenBank files are not always consistent — try searching a partial name (e.g. `-g rpo` instead of `-g rpoD`).

**"No sequences extracted" in BLAST mode**
Run this to compare IDs between your BLAST output and your GenBank file:
```bash
# IDs in your BLAST hits (column 2)
cut -f2 hits.blast6 | sort -u

# IDs in your GenBank file
grep "^LOCUS" genome.gbk | awk '{print $2}'
```
The script handles version suffix mismatches (`.1`) automatically, but other differences in naming will still cause misses.

**Warnings about records with no sequence data**
Some NCBI GenBank files contain annotation without nucleotide sequence. These records are skipped automatically. Re-download using the full GenBank flat file format (`.gbff`) to get files with sequence included. With the NCBI `datasets` CLI:
```bash
datasets download genome accession GCF_017921755.1 --include genome,gbff
```
