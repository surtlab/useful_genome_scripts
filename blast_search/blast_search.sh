#!/usr/bin/env bash
# =============================================================================
# blast_search.sh
# BLAST a query .fa file against all viable databases in a blast_db folder.
#
# Auto-detects query type (nucleotide or protein). When multiple BLAST
# programs are valid for that query type, -p must be specified.
#
#   nucleotide query → blastn  (vs _nucl dbs)   requires -p if prot dbs also present
#                    → blastx  (vs _prot dbs)
#   protein query    → blastp  (vs _prot dbs)   requires -p if nucl dbs also present
#                    → tblastn (vs _nucl dbs)
#
# Usage:
#   chmod +x blast_search.sh
#   ./blast_search.sh -q query.fa -d ./blast_db [-p blastn] [-i 80] [-o results.txt]
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
QUERY=""
DBDIR=""
OUTFILE=""
PROGRAM=""
EVALUE="1e-5"
MIN_PIDENT=""
THREADS=4
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'

usage() {
    echo "Usage: $0 -q <query.fa> -d <blast_db_dir> [-p <program>] [-i <min_pident>] [-o <out>] [-e <evalue>] [-t <threads>]"
    echo ""
    echo "  -q  Query FASTA file (nucleotide or protein)    (required)"
    echo "  -d  Path to blast_db directory                  (required)"
    echo "  -p  BLAST program: blastn, blastx, blastp, tblastn"
    echo "        Required when multiple programs are valid for the query type"
    echo "  -i  Minimum percent identity to report (0-100)  (default: all hits)"
    echo "  -o  Output file                                  (default: <query_basename>_blast_results.txt)"
    echo "  -e  E-value threshold                            (default: 1e-5)"
    echo "  -t  Number of threads                            (default: 4)"
    echo ""
    echo "Examples:"
    echo "  $0 -q my_gene.fa -d ./blast_db"
    echo "  $0 -q my_gene.fa -d ./blast_db -p blastn -i 90 -e 1e-10"
    exit 1
}

while getopts ":q:d:p:i:o:e:t:h" opt; do
    case $opt in
        q) QUERY="$OPTARG" ;;
        d) DBDIR="$OPTARG" ;;
        p) PROGRAM="$OPTARG" ;;
        i) MIN_PIDENT="$OPTARG" ;;
        o) OUTFILE="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: Option -$OPTARG requires an argument."; usage ;;
        \?) echo "ERROR: Unknown option -$OPTARG"; usage ;;
    esac
done

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ -z "$QUERY" ]]; then
    echo -e "${RED}ERROR:${NC} No query file specified. Use -q query.fa"
    usage
fi
if [[ ! -f "$QUERY" ]]; then
    echo -e "${RED}ERROR:${NC} Query file not found: $QUERY"
    exit 1
fi
if [[ -z "$DBDIR" ]]; then
    echo -e "${RED}ERROR:${NC} No database directory specified. Use -d ./blast_db"
    usage
fi
if [[ ! -d "$DBDIR" ]]; then
    echo -e "${RED}ERROR:${NC} Database directory not found: $DBDIR"
    exit 1
fi

# Validate -p if provided
if [[ -n "$PROGRAM" ]]; then
    case "$PROGRAM" in
        blastn|blastx|blastp|tblastn) ;;
        *) echo -e "${RED}ERROR:${NC} Invalid program '$PROGRAM'. Choose from: blastn blastx blastp tblastn"
           exit 1 ;;
    esac
fi

# Validate -i if provided
if [[ -n "$MIN_PIDENT" ]]; then
    if ! [[ "$MIN_PIDENT" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo -e "${RED}ERROR:${NC} -i must be a number between 0 and 100"
        exit 1
    fi
fi

# ── Dependency check ──────────────────────────────────────────────────────────
missing=()
for cmd in blastn blastp blastx tblastn python3 awk; do
    command -v "$cmd" &>/dev/null || missing+=("$cmd")
done
if [[ ${#missing[@]} -gt 0 ]]; then
    echo -e "${RED}ERROR:${NC} Missing required dependencies: ${missing[*]}"
    echo "       Install BLAST+ from: https://www.ncbi.nlm.nih.gov/books/NBK569861/"
    exit 1
fi

# ── Set default output filename ───────────────────────────────────────────────
if [[ -z "$OUTFILE" ]]; then
    query_base=$(basename "$QUERY" | sed 's/\.[^.]*$//')
    OUTFILE="${query_base}_blast_results.txt"
fi
SUMMARY_FILE="${OUTFILE%.txt}_hit_summary.tsv"
HITS_FASTA="${OUTFILE%.txt}_hit_sequences.fa"

# ── Detect query type (nucleotide or protein) ─────────────────────────────────
detect_query_type() {
    python3 << PYEOF
import re, sys

fasta = "$QUERY"
seq = []
with open(fasta) as f:
    for line in f:
        if line.startswith(">"):
            if seq:
                break
            continue
        seq.append(line.strip().upper())

seq = "".join(seq)
if not seq:
    print("unknown")
    sys.exit(0)

nuc_chars = len(re.findall(r'[ATCGNU]', seq))
total     = len(re.findall(r'[A-Z]',   seq))
nuc_frac  = nuc_chars / total if total > 0 else 0
print("nucleotide" if nuc_frac >= 0.85 else "protein")
PYEOF
}

QUERY_TYPE=$(detect_query_type)

if [[ "$QUERY_TYPE" == "unknown" ]]; then
    echo -e "${RED}ERROR:${NC} Could not determine query type from $QUERY"
    exit 1
fi

# ── Discover available databases ──────────────────────────────────────────────
NUCL_DBS=()
while IFS= read -r db; do
    [[ -n "$db" ]] && NUCL_DBS+=("$db")
done < <(find "$DBDIR" -name "*_nucl.nin" 2>/dev/null | sed 's/\.nin$//' | sort)

PROT_DBS=()
while IFS= read -r db; do
    [[ -n "$db" ]] && PROT_DBS+=("$db")
done < <(find "$DBDIR" -name "*_prot.pin" 2>/dev/null | sed 's/\.pin$//' | sort)

nucl_count=${#NUCL_DBS[@]}
prot_count=${#PROT_DBS[@]}

if [[ $nucl_count -eq 0 && $prot_count -eq 0 ]]; then
    echo -e "${RED}ERROR:${NC} No BLAST databases found in $DBDIR"
    echo "       Expected files ending in _nucl.nin or _prot.pin"
    exit 1
fi

# ── Determine valid programs for this query type + available dbs ──────────────
VALID_FOR_QUERY=()
if [[ "$QUERY_TYPE" == "nucleotide" ]]; then
    [[ $nucl_count -gt 0 ]] && VALID_FOR_QUERY+=("blastn")
    [[ $prot_count -gt 0 ]] && VALID_FOR_QUERY+=("blastx")
elif [[ "$QUERY_TYPE" == "protein" ]]; then
    [[ $prot_count -gt 0 ]] && VALID_FOR_QUERY+=("blastp")
    [[ $nucl_count -gt 0 ]] && VALID_FOR_QUERY+=("tblastn")
fi

valid_count=${#VALID_FOR_QUERY[@]}

# If -p was given, verify it's appropriate for this query type + available dbs
if [[ -n "$PROGRAM" ]]; then
    match=false
    for v in "${VALID_FOR_QUERY[@]}"; do
        [[ "$v" == "$PROGRAM" ]] && match=true && break
    done
    if [[ "$match" == false ]]; then
        echo -e "${RED}ERROR:${NC} '$PROGRAM' is not valid for a ${QUERY_TYPE} query against the available databases."
        echo "       Valid options: ${VALID_FOR_QUERY[*]}"
        exit 1
    fi
fi

# If multiple programs are valid and -p was not specified, require it
if [[ $valid_count -gt 1 && -z "$PROGRAM" ]]; then
    echo -e "${RED}ERROR:${NC} Multiple BLAST programs are valid for this ${QUERY_TYPE} query"
    echo "       and both nucl and prot databases are present."
    echo ""
    echo "       Please specify one with -p:"
    for v in "${VALID_FOR_QUERY[@]}"; do
        case "$v" in
            blastn)  echo "         -p blastn   nucleotide query vs nucleotide dbs" ;;
            blastx)  echo "         -p blastx   nucleotide query (translated) vs protein dbs" ;;
            blastp)  echo "         -p blastp   protein query vs protein dbs" ;;
            tblastn) echo "         -p tblastn  protein query vs nucleotide dbs (translated)" ;;
        esac
    done
    exit 1
fi

# If only one valid program, auto-select it
if [[ $valid_count -eq 1 && -z "$PROGRAM" ]]; then
    PROGRAM="${VALID_FOR_QUERY[0]}"
    echo -e "${YELLOW}[INFO]${NC} Auto-selected program: ${PROGRAM}"
fi

# Pick the right database set for the chosen program
case "$PROGRAM" in
    blastn|tblastn) TARGET_DBS=("${NUCL_DBS[@]}") ;;
    blastp|blastx)  TARGET_DBS=("${PROT_DBS[@]}") ;;
esac

# ── Print banner ──────────────────────────────────────────────────────────────
pident_display="${MIN_PIDENT:-(all hits)}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  BLAST Search"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Query      : $QUERY"
echo "  Query type : $QUERY_TYPE"
echo "  Program    : $PROGRAM"
echo "  DB dir     : $DBDIR"
echo "  Databases  : ${#TARGET_DBS[@]}"
echo "  E-value    : $EVALUE"
echo "  Min pident : $pident_display"
echo "  Threads    : $THREADS"
echo "  Output     : $OUTFILE"
echo "  Hit table  : $SUMMARY_FILE"
echo "  Hit seqs   : $HITS_FASTA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# ── Write output file header ──────────────────────────────────────────────────
{
    echo "# BLAST Search Results"
    echo "# Query:      $QUERY"
    echo "# Query type: $QUERY_TYPE"
    echo "# Program:    $PROGRAM"
    echo "# DB dir:     $DBDIR"
    echo "# E-value:    $EVALUE"
    echo "# Min pident: $pident_display"
    echo "# Generated:  $(date)"
    echo "# Columns:    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
} > "$OUTFILE"

# ── Initialise hit summary table ──────────────────────────────────────────────
echo -e "accession\thit" > "$SUMMARY_FILE"
# Initialise hit sequences FASTA
> "$HITS_FASTA"

# ── Extract hit sequences from BLAST results ──────────────────────────────────
# For protein programs (blastp/blastx): look up sseqid in the .faa file.
# For nucleotide programs (blastn/tblastn): extract sstart..send from the .fna,
# reverse-complementing if sstart > send (minus strand hit).
extract_hit_sequences() {
    local program="$1"   # blastn | blastp | blastx | tblastn
    local db="$2"        # path to db base (without _nucl/_prot suffix)
    local hits_file="$3" # filtered outfmt6 results
    local out_fasta="$4" # append sequences here

    python3 << PYEOF
import re, sys, os

program   = "$program"
db        = "$db"
hits_file = "$hits_file"
out_fasta = "$out_fasta"

# Read hit table: col 1=qseqid, col 1=sseqid, col 8=sstart, col 9=send
def normalize_id(seq_id):
    """Strip NCBI pipe-delimited wrappers like ref|WP_013027259.1| → WP_013027259.1"""
    if "|" in seq_id:
        parts = seq_id.split("|")
        # Format is type|accession| — return the accession part
        if len(parts) >= 2 and parts[-1] == "":
            return parts[-2]
        return parts[-1]
    return seq_id

hits = []
with open(hits_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 10:
            continue
        sseqid = normalize_id(cols[1])
        sstart = int(cols[8])
        send   = int(cols[9])
        pident = cols[2]
        evalue = cols[10]
        hits.append((sseqid, sstart, send, pident, evalue))

if not hits:
    sys.exit(0)

# ── Helper: load all sequences from a FASTA file into a dict ──────────────────
def load_fasta(path):
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                # Store full header but key by first word
                header = line[1:]
                current_id = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs

# ── Helper: reverse complement ────────────────────────────────────────────────
def revcomp(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

# ── Protein programs: look up sseqid in the .faa ─────────────────────────────
if program in ("blastp", "blastx"):
    faa = re.sub(r'_prot$', '', db) + ".faa"
    if not os.path.exists(faa):
        print(f"# WARNING: {faa} not found, cannot extract protein sequences", file=sys.stderr)
        sys.exit(0)
    seqs = load_fasta(faa)
    seen = set()
    with open(out_fasta, "a") as out:
        for (sseqid, sstart, send, pident, evalue) in hits:
            if sseqid in seen:
                continue
            seen.add(sseqid)
            seq = seqs.get(sseqid)
            if seq is None:
                print(f"# WARNING: {sseqid} not found in {faa}", file=sys.stderr)
                continue
            out.write(f">{sseqid} pident={pident} evalue={evalue}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

# ── Nucleotide programs: extract coordinate range from the .fna ───────────────
elif program in ("blastn", "tblastn"):
    fna = re.sub(r'_nucl$', '', db) + ".fna"
    if not os.path.exists(fna):
        print(f"# WARNING: {fna} not found, cannot extract nucleotide sequences", file=sys.stderr)
        sys.exit(0)
    seqs = load_fasta(fna)
    seen = set()
    with open(out_fasta, "a") as out:
        for (sseqid, sstart, send, pident, evalue) in hits:
            # Use (sseqid, sstart, send) as key to allow same contig hit at diff coords
            hit_key = (sseqid, sstart, send)
            if hit_key in seen:
                continue
            seen.add(hit_key)
            seq = seqs.get(sseqid)
            if seq is None:
                print(f"# WARNING: {sseqid} not found in {fna}", file=sys.stderr)
                continue
            # BLAST coords are 1-based; sstart > send means minus strand
            if sstart <= send:
                subseq = seq[sstart - 1:send]
                strand = "plus"
            else:
                subseq = revcomp(seq[send - 1:sstart])
                strand = "minus"
            out.write(f">{sseqid}:{sstart}-{send}({strand}) pident={pident} evalue={evalue}\n")
            for i in range(0, len(subseq), 60):
                out.write(subseq[i:i+60] + "\n")
PYEOF
}

# ── Run BLAST function ────────────────────────────────────────────────────────
run_blast() {
    local program="$1"
    local db="$2"
    local db_label
    db_label=$(basename "$db")
    local tmpout
    tmpout=$(mktemp /tmp/blast_XXXXXX.txt)

    echo -e "  ${CYAN}[RUN]${NC}   ${program} vs ${db_label}..."

    if "$program" \
        -query "$QUERY" \
        -db "$db" \
        -out "$tmpout" \
        -evalue "$EVALUE" \
        -num_threads "$THREADS" \
        -outfmt "$OUTFMT" \
        2>/dev/null; then

        # Apply percent identity filter if set
        local filtered
        filtered=$(mktemp /tmp/blast_filtered_XXXXXX.txt)
        if [[ -n "$MIN_PIDENT" ]]; then
            awk -v pid="$MIN_PIDENT" '$3 >= pid' "$tmpout" > "$filtered"
        else
            cp "$tmpout" "$filtered"
        fi

        local hits
        hits=$(grep -c $'^\S' "$filtered" 2>/dev/null || true)

        {
            echo ""
            echo "# ────────────────────────────────────────────────────────────"
            echo "# Program:  ${program}"
            echo "# Database: ${db_label}"
            echo "# Hits:     ${hits}"
            [[ -n "$MIN_PIDENT" ]] && echo "# Filter:   pident >= ${MIN_PIDENT}%"
            echo "# ────────────────────────────────────────────────────────────"
            if [[ $hits -gt 0 ]]; then
                cat "$filtered"
            else
                echo "# (no hits)"
            fi
        } >> "$OUTFILE"

        if [[ $hits -gt 0 ]]; then
            echo -e "  ${GREEN}[OK]${NC}    ${hits} hit(s)"
            echo -e "${db_label}\tyes" >> "$SUMMARY_FILE"
            extract_hit_sequences "$program" "$db" "$filtered" "$HITS_FASTA"
        else
            echo -e "  ${YELLOW}[OK]${NC}    no hits"
            echo -e "${db_label}\tno" >> "$SUMMARY_FILE"
        fi

        rm -f "$filtered"
    else
        echo -e "  ${RED}[FAIL]${NC}  ${program} vs ${db_label} failed"
        echo -e "${db_label}\terror" >> "$SUMMARY_FILE"
        {
            echo ""
            echo "# ────────────────────────────────────────────────────────────"
            echo "# Program:  ${program}"
            echo "# Database: ${db_label}"
            echo "# ERROR: BLAST failed for this database"
            echo "# ────────────────────────────────────────────────────────────"
        } >> "$OUTFILE"
    fi

    rm -f "$tmpout"
}

# ── Run against all target databases ─────────────────────────────────────────
total_runs=0
echo -e "${YELLOW}── ${PROGRAM} ──${NC}"
for db in "${TARGET_DBS[@]}"; do
    run_blast "$PROGRAM" "$db"
    ((total_runs++)) || true
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Program      : $PROGRAM"
echo "  Searches run : $total_runs"
echo "  Results file : $OUTFILE"
echo "  Hit table    : $SUMMARY_FILE"
echo "  Hit seqs     : $HITS_FASTA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
