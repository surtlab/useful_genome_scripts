# download_genbank.sh

A Bash script for downloading bacterial genome sequences from NCBI in GenBank (`.gbk`) format, given a list of accession numbers. Optionally builds BLAST databases from each downloaded genome.

---

## What's new in v1.3

### Bug fix: spaces in output directory paths
`makeblastdb` and `blastn` have a known bug where they tokenize their own argument strings on spaces, ignoring the OS-level quoting that bash provides. This means paths containing spaces — such as iCloud Drive paths (`~/Library/Mobile Documents/...`) — would cause both database building and querying to silently fail with errors like:

```
BLAST options error: File /Users/you/Library/Mobile does not exist
```

**Fix:** The script now stages all `makeblastdb` operations through a guaranteed space-free temp directory (`/tmp/blastdb_XXXXXX`), builds the databases there, then moves the finished index files to the real output directory. This means `-b` now works correctly regardless of whether your output path contains spaces.

**For `blastn`/`blastp` queries:** If your blast_db directory is inside a path with spaces, you must still work around the same BLAST bug yourself. The recommended approach is to `cd` into the blast_db directory first and use a relative db name:

```bash
cd "/Users/you/Library/Mobile Documents/.../blast_db"
blastn -db GCF_017921755.1 -query my_query.fna -evalue 1e-5 -outfmt 6
```

Or use a symlink to avoid the issue permanently:

```bash
ln -s "/Users/you/Library/Mobile Documents/.../blast_db" ~/blast_db
blastn -db ~/blast_db/GCF_017921755.1 -query my_query.fna -evalue 1e-5 -outfmt 6
```

### Change: BLAST database naming (no `_nucl` / `_prot` suffix)

Previous versions named BLAST databases with a `_nucl` or `_prot` suffix:
```
GCF_017921755.1_nucl.nhr   ← old
GCF_017921755.1_prot.phr   ← old
```

As of v1.3, databases are named using just the accession, matching the convention of building manually with `makeblastdb`:
```
GCF_017921755.1.nhr   ← new
GCF_017921755.1.phr   ← new
```

This means `-db` now takes the plain accession name with no suffix:
```bash
# v1.3 and later
blastn -db ./blast_db/GCF_017921755.1 -query my_query.fna -outfmt 6

# Previous versions required
blastn -db ./blast_db/GCF_017921755.1_nucl -query my_query.fna -outfmt 6
```

> **Note:** If you built databases with a previous version of this script, rebuild them by re-running with `-b` against your existing `.gbk` files. Old `_nucl`/`_prot` databases will not be automatically removed.

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

## Running on Windows

This script requires a bash environment and Unix tools that are not available natively on Windows. The recommended solution is to use **Windows Subsystem for Linux (WSL)**, which provides a full Linux environment inside Windows.

### Step 1: Install WSL (one-time setup)

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This installs WSL with Ubuntu by default. Restart your computer when prompted.

### Step 2: Open WSL

Search for **Ubuntu** in the Start menu and open it. The first launch will ask you to create a Linux username and password.

### Step 3: Install dependencies inside WSL

```bash
sudo apt update
sudo apt install curl python3
```

BLAST+ (only needed if using `-b`):
```bash
sudo apt install ncbi-blast+
```

### Step 4: Access your files

Your Windows files are accessible inside WSL under `/mnt/c/`. For example, if your accessions file is on your Desktop:

```bash
cd /mnt/c/Users/YourWindowsUsername/Desktop
```

### Step 5: Run the script

```bash
chmod +x download_genbank.sh
./download_genbank.sh -i accessions.txt -o ./sequences -e your@email.com
```

> **Note:** If you have conda installed on Windows, you can also run the script from **Git Bash** (included with Git for Windows). However, WSL is more reliable for scripts that use `mktemp` and `/tmp/` paths.

---

## Running on macOS

macOS comes with bash and most required Unix tools pre-installed, so setup is straightforward.

### Step 1: Check your bash version

macOS ships with bash 3.2 by default (due to licensing), but this script requires bash 4.0+. Check your version:

```bash
bash --version
```

If the version is below 4.0, install an updated bash via Homebrew:

```bash
brew install bash
```

If you don't have Homebrew installed, you can install it with:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### Step 2: Check for required tools

`curl` and `python3` are pre-installed on modern macOS. Verify with:

```bash
curl --version
python3 --version
```

If `python3` is missing, install it via Homebrew:

```bash
brew install python3
```

### Step 3: Install BLAST+ (only needed if using `-b`)

Via Homebrew:

```bash
brew install blast
```

Or via conda:

```bash
conda install -c bioconda blast
```

### Step 4: Run the script

Open Terminal, navigate to the folder containing the script, and run:

```bash
chmod +x download_genbank.sh
./download_genbank.sh -i accessions.txt -o ./sequences -e your@email.com
```

> **Note for iCloud Drive users:** If your output directory is inside iCloud Drive (e.g. `~/Library/Mobile Documents/com~apple~CloudDocs/...`), the `-b` flag will still work correctly as of v1.3. However, running `blastn`/`blastp` against databases in that path requires the `cd` workaround described in the "What's new" section above, or a symlink to a space-free location.

> **Note:** If you are using conda, make sure your environment is activated before running the script so that all dependencies are available: `conda activate your-environment`

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
   - `.fna` → nucleotide database (`<accession>.n*`) — usable with `blastn`, `tblastn`, `tblastx`
   - `.faa` → protein database (`<accession>.p*`) — usable with `blastp`, `blastx`

Both databases are built if the genome is fully annotated. If a genome has no CDS annotations (e.g. a draft WGS assembly), only the nucleotide database is built and the protein step is silently skipped.

Duplicate protein IDs (common in WGS assemblies where the same `WP_` accession appears across multiple contigs) are deduplicated automatically — only the first occurrence is written to the `.faa`.

All BLAST databases and intermediate FASTA files are written to `<outdir>/blast_db/`.

Example BLAST usage after building:

```bash
# Nucleotide search
blastn -query my_query.fna -db ./genomes/blast_db/GCF_017921755.1 -out results.txt

# Protein search
blastp -query my_query.faa -db ./genomes/blast_db/GCF_017921755.1 -out results.txt
```

---

## Retry behaviour

The script retries failed downloads up to **5 times** using exponential backoff with jitter (roughly 5s, 20s, 45s, 80s between attempts). This handles transient NCBI timeouts without hammering the server.

Already-downloaded files are skipped automatically, so it is safe to re-run the script on the same output directory.

---

## NCBI API rate limits

NCBI allows up to 3 requests/second without an API key, and 10 requests/second with one. The script uses a 0.4s inter-request delay by default to stay within the unauthenticated limit. If you have an [NCBI API key](https://www.ncbi.nlm.nih.gov/account/), you can reduce the delay and speed up large batch downloads by editing the `DELAY` variable at the top of the script.
