# download_genbank.sh

A Bash script for downloading bacterial genome sequences from NCBI in GenBank (`.gbk`) format, given a list of accession numbers. Optionally builds BLAST databases from each downloaded genome.

---

## Requirements

### Always required
- `bash` 4.0+
- `curl` or `wget`
- `python3`
- Standard Unix tools: `grep`, `sed`, `tr`, `mktemp`, `head`

### Required only with `-b` (BLAST database building)
- `makeblastdb` from the [NCBI BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) package

The script checks all relevant dependencies at startup and exits with a clear error if anything is missing.

---

## Usage

```bash
chmod +x download_genbank.sh
./download_genbank.sh -i accessions.txt -o ./sequences -e your@email.com
```

To also build BLAST databases:

```bash
./download_genbank.sh -i accessions.txt -o ./genomes -e your@email.com -b
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Path to accessions text file **(required)** | — |
| `-o` | Output directory | `./sequences` |
| `-e` | Your email address for NCBI Entrez | `your.email@example.com` |
| `-b` | Build BLAST databases after downloading | off |
| `-h` | Show help | — |

> **Note:** NCBI's Entrez API requires a real email address to identify requests. Always pass your email with `-e`.

---

## Input file format

A plain text file with one accession per line. Lines starting with `#` and blank lines are ignored.

```
# My genomes
NC_000913.3
GCF_017921755.1
JABEBS01
```

---

## Supported accession types

| Type | Example | Method |
|------|---------|--------|
| Standard RefSeq / GenBank | `NC_000913`, `NZ_CP012345`, `CP001234` | Direct `efetch` |
| WGS master records | `JABEBS01`, `ABCD02` | `esearch` → `efetch` |
| Assembly accessions | `GCF_017921755.1`, `GCA_017921755.1` | `esearch` (by assembly tag) → `efetch` |

**WGS records:** The script queries all contig UIDs associated with the WGS prefix and downloads them as a single concatenated `.gbk` file.

**Assembly accessions (GCF/GCA):** The script uses the `[Assembly Accession]` field tag in `esearch` to retrieve only sequences directly associated with that assembly, avoiding spurious linked records from other organisms.

---

## Output

Each accession is saved as `<ACCESSION>.gbk` in the output directory. For multi-contig assemblies, all sequences are concatenated into a single file.

If any downloads fail after all retries, their accessions are saved to `failed_accessions.txt` in the output directory.

To retry only the failed accessions:

```bash
./download_genbank.sh -i ./sequences/failed_accessions.txt -o ./sequences -e your@email.com
```

---

## BLAST database building (`-b`)

When `-b` is passed, the script inspects each successfully downloaded `.gbk` and builds whichever databases the annotation supports:

1. The `.gbk` is parsed to extract:
   - **Nucleotide sequences** from `ORIGIN` blocks → `<accession>.fna`
   - **Protein sequences** from `/translation=` tags in CDS features → `<accession>.faa`
2. `makeblastdb` is run on each non-empty file:
   - `.fna` → nucleotide database (`_nucl.*`) — usable with `blastn`, `tblastn`, `tblastx`
   - `.faa` → protein database (`_prot.*`) — usable with `blastp`, `blastx`

Both databases are built if the genome is fully annotated. If a genome has no CDS annotations (e.g. a draft WGS assembly), only the nucleotide database is built and the protein step is silently skipped. Partial annotation is handled gracefully — if only some contigs have `/translation=` tags, proteins are extracted from those contigs and the rest are ignored.

Duplicate protein IDs (common in WGS assemblies where the same `WP_` accession appears across multiple contigs) are deduplicated automatically — only the first occurrence is written to the `.faa`.

All BLAST databases and intermediate FASTA files are written to `<outdir>/blast_db/`.

Example BLAST usage after building:

```bash
# Nucleotide search
blastn -query my_query.fna -db ./genomes/blast_db/GCF_017921755.1_nucl -out results.txt

# Protein search
blastp -query my_query.faa -db ./genomes/blast_db/GCF_017921755.1_prot -out results.txt
```

---

## Retry behaviour

The script retries failed downloads up to **5 times** using exponential backoff with jitter (roughly 5s, 20s, 45s, 80s between attempts). This handles transient NCBI timeouts without hammering the server.

Already-downloaded files are skipped automatically, so it is safe to re-run the script on the same output directory.

---

## NCBI API rate limits

NCBI allows up to 3 requests/second without an API key, and 10 requests/second with one. The script uses a 0.4s inter-request delay by default to stay within the unauthenticated limit. If you have an [NCBI API key](https://www.ncbi.nlm.nih.gov/account/), you can reduce the delay and speed up large batch downloads by editing the `DELAY` variable at the top of the script.
