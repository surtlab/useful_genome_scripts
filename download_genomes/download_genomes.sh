#!/usr/bin/env bash
# =============================================================================
# download_genbank.sh
# Download bacterial genome sequences from NCBI GenBank (.gbk format)
# using a text file of accession numbers.
#
# Supports:
#   - Standard accessions (NC_, NZ_, CP_, etc.)  → direct efetch
#   - WGS master records (e.g. JABEBS01)         → esearch → efetch (2-step)
#
# Usage:
#   chmod +x download_genbank.sh
#   ./download_genbank.sh -i accessions.txt -o ./sequences -e your@email.com
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
EMAIL="your.email@example.com"
OUTDIR="./sequences"
INPUT=""
RETRIES=5
DELAY=0.4
BUILD_BLAST=false

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

usage() {
    echo "Usage: $0 -i <accessions.txt> [-o <outdir>] [-e <email>] [-b]"
    echo ""
    echo "  -i  Path to accessions text file (required)"
    echo "  -o  Output directory              (default: ./sequences)"
    echo "  -e  Your email for NCBI Entrez    (default: ${EMAIL})"
    echo "  -b  Build BLAST databases after downloading"
    echo ""
    echo "Example:"
    echo "  $0 -i accessions.txt -o ./genomes -e jane@example.com -b"
    exit 1
}

while getopts ":i:o:e:bh" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        e) EMAIL="$OPTARG" ;;
        b) BUILD_BLAST=true ;;
        h) usage ;;
        :) echo "ERROR: Option -$OPTARG requires an argument."; usage ;;
        \?) echo "ERROR: Unknown option -$OPTARG"; usage ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo -e "${RED}ERROR:${NC} No input file specified. Use -i accessions.txt"
    usage
fi

if [[ ! -f "$INPUT" ]]; then
    echo -e "${RED}ERROR:${NC} Input file not found: $INPUT"
    exit 1
fi

if [[ "$EMAIL" == "your.email@example.com" ]]; then
    echo -e "${YELLOW}WARNING:${NC} Using placeholder email. NCBI requires a real email."
    echo "         Pass your email with: -e your.email@example.com"
    echo ""
fi

if command -v curl &>/dev/null; then
    DOWNLOADER="curl"
elif command -v wget &>/dev/null; then
    DOWNLOADER="wget"
else
    echo -e "${RED}ERROR:${NC} Neither curl nor wget found. Please install one."
    exit 1
fi

# ── Dependency check ──────────────────────────────────────────────────────────
missing=()
for cmd in grep sed tr mktemp head python3; do
    command -v "$cmd" &>/dev/null || missing+=("$cmd")
done
if [[ ${#missing[@]} -gt 0 ]]; then
    echo -e "${RED}ERROR:${NC} Missing required dependencies: ${missing[*]}"
    echo "       Please install them and re-run."
    exit 1
fi

if [[ "$BUILD_BLAST" == true ]] && ! command -v makeblastdb &>/dev/null; then
    echo -e "${RED}ERROR:${NC} -b was specified but makeblastdb not found."
    echo "       Install BLAST+ from: https://www.ncbi.nlm.nih.gov/books/NBK569861/"
    exit 1
fi

mkdir -p "$OUTDIR"

ESEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ── Detect if accession is a GCF/GCA assembly accession ──────────────────────
is_assembly() {
    [[ "$1" =~ ^GC[AF]_[0-9]+(\.[0-9]+)?$ ]]
}

# ── Resolve GCF/GCA → nuccore UIDs via targeted esearch ──────────────────────
resolve_assembly_uids() {
    local acc="$1"
    local term="${acc}%5BAssembly+Accession%5D"
    local search_url="${ESEARCH_URL}?db=nuccore&term=${term}&retmode=xml&retmax=1000&email=${EMAIL}&tool=download_genbank"
    local tmpxml
    tmpxml=$(mktemp /tmp/esearch_XXXXXX.xml)

    if ! http_get "$search_url" "$tmpxml"; then
        rm -f "$tmpxml"
        echo ""
        return 1
    fi

    local uids
    uids=$(grep -o '<Id>[0-9]*</Id>' "$tmpxml" | grep -o '[0-9]*' | tr '\n' ',' | sed 's/,$//' || true)
    local count
    count=$(grep -o '<Count>[0-9]*</Count>' "$tmpxml" | head -1 | grep -o '[0-9]*' || true)
    rm -f "$tmpxml"

    echo -e "  ${YELLOW}[INFO]${NC}  esearch returned ${count} nuccore UID(s) for ${acc}" >&2
    echo "$uids"
}

# ── Convert .gbk → nucleotide + protein FASTA ────────────────────────────────
gbk_to_fasta() {
    local gbk="$1"
    local base="$2"

    python3 - "$gbk" "$base" <<'PYEOF'
import sys, re

gbk_file = sys.argv[1]
base     = sys.argv[2]

# Split into individual records on //
records = []
current = []
with open(gbk_file) as f:
    for line in f:
        current.append(line)
        if line.strip() == "//":
            records.append("".join(current))
            current = []

nuc_out  = open(base + ".fna", "w")
prot_out = open(base + ".faa", "w")
seen_ids = set()  # track protein IDs to skip duplicates across contigs

for rec in records:
    # Accession and definition for nucleotide header
    acc  = ""
    defn = ""
    for line in rec.splitlines():
        if line.startswith("ACCESSION"):
            acc = line.split()[1]
        if line.startswith("DEFINITION"):
            defn = line[12:].strip()
    if not acc:
        continue

    # Nucleotide sequence from ORIGIN block
    in_origin = False
    nuc_seq   = []
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
        nuc_out.write(f">{acc} {defn}\n")
        for i in range(0, len(nuc_str), 60):
            nuc_out.write(nuc_str[i:i+60] + "\n")

    # CDS protein translations — use a line-by-line state machine
    # to handle multi-line /translation= values reliably
    in_translation = False
    prot_id  = None
    product  = None
    trans_lines = []
    in_features  = False

    for line in rec.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            # flush any open translation
            if in_translation and prot_id and trans_lines:
                seq = re.sub(r'[^A-Za-z]', '', "".join(trans_lines))
                if prot_id not in seen_ids:
                    seen_ids.add(prot_id)
                    hdr = prot_id + (f" {product}" if product else "")
                    prot_out.write(f">{hdr}\n")
                    for i in range(0, len(seq), 60):
                        prot_out.write(seq[i:i+60] + "\n")
            in_translation = False
            in_features = False
            continue

        if not in_features:
            continue

        # New top-level feature resets state
        if re.match(r'     \S', line):
            if in_translation and prot_id and trans_lines:
                seq = re.sub(r'[^A-Za-z]', '', "".join(trans_lines))
                if prot_id not in seen_ids:
                    seen_ids.add(prot_id)
                    hdr = prot_id + (f" {product}" if product else "")
                    prot_out.write(f">{hdr}\n")
                    for i in range(0, len(seq), 60):
                        prot_out.write(seq[i:i+60] + "\n")
            in_translation = False
            prot_id = None
            product = None
            trans_lines = []

        m_pid  = re.match(r'\s+/protein_id="([^"]+)"', line)
        m_prod = re.match(r'\s+/product="([^"]+)"', line)
        m_tr   = re.match(r'\s+/translation="([^"]*)"?', line)

        if m_pid:
            prot_id = m_pid.group(1)
        if m_prod:
            product = m_prod.group(1)
        if m_tr:
            in_translation = True
            trans_lines = [m_tr.group(1)]
            if line.rstrip().endswith('"'):
                # single-line translation
                in_translation = False
                if prot_id and trans_lines:
                    seq = re.sub(r'[^A-Za-z]', '', "".join(trans_lines))
                    if prot_id not in seen_ids:
                        seen_ids.add(prot_id)
                        hdr = prot_id + (f" {product}" if product else "")
                        prot_out.write(f">{hdr}\n")
                        for i in range(0, len(seq), 60):
                            prot_out.write(seq[i:i+60] + "\n")
                trans_lines = []
        elif in_translation:
            stripped = line.strip()
            if stripped.endswith('"'):
                trans_lines.append(stripped[:-1])
                in_translation = False
                if prot_id and trans_lines:
                    seq = re.sub(r'[^A-Za-z]', '', "".join(trans_lines))
                    if prot_id not in seen_ids:
                        seen_ids.add(prot_id)
                        hdr = prot_id + (f" {product}" if product else "")
                        prot_out.write(f">{hdr}\n")
                        for i in range(0, len(seq), 60):
                            prot_out.write(seq[i:i+60] + "\n")
                trans_lines = []
            else:
                trans_lines.append(stripped)

nuc_out.close()
prot_out.close()
PYEOF
}

# ── Build BLAST databases for one accession ───────────────────────────────────
# nucl and prot are tracked independently — a genome with no annotated CDSs
# (e.g. a draft WGS assembly) will simply skip the prot db without failing.
make_blast_db() {
    local acc="$1"
    local gbk="${OUTDIR}/${acc}.gbk"
    local blastdir="${OUTDIR}/blast_db"
    local dbbase="${blastdir}/${acc}"

    mkdir -p "$blastdir"
    gbk_to_fasta "$gbk" "$dbbase"

    # Nucleotide db
    if [[ -s "${dbbase}.fna" ]]; then
        if makeblastdb -in "${dbbase}.fna" -dbtype nucl -out "${dbbase}_nucl" \
                -title "$acc" -parse_seqids &>/dev/null; then
            echo -e "  ${GREEN}[BLAST]${NC} ${acc} → nucl db built"
        else
            echo -e "  ${RED}[BLAST]${NC} ${acc} — nucl db failed"
            return 1
        fi
    else
        echo -e "  ${YELLOW}[BLAST]${NC} ${acc} — no nucleotide sequences, skipping nucl db"
        return 1
    fi

    # Protein db — optional, skip silently if no translations found
    if [[ -s "${dbbase}.faa" ]]; then
        if makeblastdb -in "${dbbase}.faa" -dbtype prot -out "${dbbase}_prot" \
                -title "$acc" -parse_seqids &>/dev/null; then
            echo -e "  ${GREEN}[BLAST]${NC} ${acc} → prot db built"
        else
            echo -e "  ${YELLOW}[BLAST]${NC} ${acc} — prot db failed (non-fatal)"
        fi
    else
        echo -e "  ${YELLOW}[BLAST]${NC} ${acc} — no protein annotations found, skipping prot db"
    fi

    return 0
}

# ── Detect if accession is a WGS master record ────────────────────────────────
# WGS pattern: 4–6 uppercase letters + 2 digits (optionally .version), e.g. JABEBS01, ABCD01
is_wgs() {
    [[ "$1" =~ ^[A-Z]{4,6}[0-9]{2}(\.[0-9]+)?$ ]]
}

# ── Detect if accession is a GCF/GCA assembly accession ──────────────────────
is_assembly() {
    [[ "$1" =~ ^GC[AF]_[0-9]+(\.[0-9]+)?$ ]]
}

# ── Resolve GCF/GCA → nuccore UIDs via targeted esearch ──────────────────────
# Uses "GCF_XXXXXX[Assembly Accession]" field tag to return only sequences
# directly tagged with that assembly — no spurious linked records.
# Brackets are URL-encoded (%5B/%5D) to prevent curl misinterpreting them
# as range expressions.
resolve_assembly_uids() {
    local acc="$1"
    local term="${acc}%5BAssembly+Accession%5D"
    local search_url="${ESEARCH_URL}?db=nuccore&term=${term}&retmode=xml&retmax=1000&email=${EMAIL}&tool=download_genbank"
    local tmpxml
    tmpxml=$(mktemp /tmp/esearch_XXXXXX.xml)

    if ! http_get "$search_url" "$tmpxml"; then
        rm -f "$tmpxml"
        echo ""
        return 1
    fi

    local uids
    uids=$(grep -o '<Id>[0-9]*</Id>' "$tmpxml" | grep -o '[0-9]*' | tr '\n' ',' | sed 's/,$//' || true)
    local count
    count=$(grep -o '<Count>[0-9]*</Count>' "$tmpxml" | head -1 | grep -o '[0-9]*' || true)
    rm -f "$tmpxml"

    echo -e "  ${YELLOW}[INFO]${NC}  esearch returned ${count} nuccore UID(s) for ${acc}" >&2
    echo "$uids"
}

# ── Low-level HTTP fetch ──────────────────────────────────────────────────────
http_get() {
    local url="$1"
    local outfile="$2"
    local http_status

    if [[ "$DOWNLOADER" == "curl" ]]; then
        http_status=$(curl -sS --globoff --max-time 300 --connect-timeout 30 \
            -w "%{http_code}" -o "$outfile" "$url" 2>/tmp/dl_err)
    else
        wget -q --timeout=300 -O "$outfile" "$url" 2>/tmp/dl_err
        http_status=$?
    fi

    local err
    err=$(cat /tmp/dl_err 2>/dev/null || true)
    [[ -n "$err" ]] && echo -e "  ${YELLOW}[DEBUG]${NC} network msg: ${err}"

    # For curl, check HTTP status code
    if [[ "$DOWNLOADER" == "curl" && "$http_status" != "200" ]]; then
        echo -e "  ${YELLOW}[DEBUG]${NC} HTTP status: ${http_status}"
        return 1
    fi

    return 0
}

# ── esearch: resolve accession → all UIDs (comma-separated) ──────────────────
# For WGS masters, querying the prefix (e.g. JABEBS01) returns ALL contig UIDs,
# not just the master record — this is how the working Python/EDirect script does it.
# retmax=500 handles up to 500 contigs; raise if your assemblies are larger.
resolve_all_uids() {
    local acc="$1"
    local search_url="${ESEARCH_URL}?db=nuccore&term=${acc}&retmode=xml&retmax=500&email=${EMAIL}&tool=download_genbank"
    local tmpxml
    tmpxml=$(mktemp /tmp/esearch_XXXXXX.xml)

    if ! http_get "$search_url" "$tmpxml"; then
        rm -f "$tmpxml"
        echo ""
        return 1
    fi

    # Extract all <Id> values and join as comma-separated list for efetch
    local uids
    uids=$(grep -o '<Id>[0-9]*</Id>' "$tmpxml" | grep -o '[0-9]*' | tr '\n' ',' | sed 's/,$//' || true)
    local count
    count=$(grep -o '<Count>[0-9]*</Count>' "$tmpxml" | head -1 | grep -o '[0-9]*' || true)
    rm -f "$tmpxml"

    echo -e "  ${YELLOW}[INFO]${NC}  esearch returned ${count} UID(s)" >&2
    echo "$uids"
}

# ── Main download function ────────────────────────────────────────────────────
download_accession() {
    local acc="$1"
    local outfile="${OUTDIR}/${acc}.gbk"

    if [[ -f "$outfile" && -s "$outfile" ]]; then
        echo -e "  ${YELLOW}[SKIP]${NC}  ${acc} already exists"
        return 0
    fi

    local fetch_url
    if is_assembly "$acc"; then
        echo -e "  ${YELLOW}[INFO]${NC}  ${acc} is an assembly accession — querying nuccore by assembly tag..."
        local uids
        uids=$(resolve_assembly_uids "$acc")
        sleep "$DELAY"

        if [[ -z "$uids" ]]; then
            echo -e "  ${RED}[FAIL]${NC}  ${acc} — esearch returned no UIDs. Check the accession is correct."
            echo "$acc" >> "${OUTDIR}/failed_accessions.txt"
            return 1
        fi
        echo -e "  ${YELLOW}[INFO]${NC}  Fetching all sequences as combined .gbk..."
        fetch_url="${EFETCH_URL}?db=nuccore&id=${uids}&rettype=gbwithparts&retmode=text&email=${EMAIL}&tool=download_genbank"
    elif is_wgs "$acc"; then
        # WGS: esearch returns ALL contig UIDs for the prefix query.
        # We pass them all to efetch as a comma-separated list → one combined .gbk
        echo -e "  ${YELLOW}[INFO]${NC}  ${acc} is a WGS master — resolving all contig UIDs via esearch..."
        local uids
        uids=$(resolve_all_uids "$acc")
        sleep "$DELAY"   # be polite between the two API calls

        if [[ -z "$uids" ]]; then
            echo -e "  ${RED}[FAIL]${NC}  ${acc} — esearch returned no UIDs. Check the accession is correct."
            echo "$acc" >> "${OUTDIR}/failed_accessions.txt"
            return 1
        fi
        echo -e "  ${YELLOW}[INFO]${NC}  Fetching all contigs as combined .gbk..."
        # Pass all UIDs at once; efetch returns one GenBank record per contig, concatenated
        fetch_url="${EFETCH_URL}?db=nuccore&id=${uids}&rettype=gbwithparts&retmode=text&email=${EMAIL}&tool=download_genbank"
    else
        # Standard accession: direct efetch
        fetch_url="${EFETCH_URL}?db=nuccore&id=${acc}&rettype=gbwithparts&retmode=text&email=${EMAIL}&tool=download_genbank"
    fi

    local attempt=0
    local success=false

    while [[ $attempt -lt $RETRIES ]]; do
        ((attempt++))

        # FIX: capture http_get return code explicitly rather than letting it silently pass
        if http_get "$fetch_url" "$outfile"; then
            if [[ -s "$outfile" ]] && grep -q "^LOCUS" "$outfile" 2>/dev/null; then
                echo -e "  ${GREEN}[OK]${NC}    ${acc}  →  ${outfile}"
                success=true
                break
            fi
        fi

        # If we reach here, the download didn't produce a valid GenBank file
        if [[ -f "$outfile" ]]; then
            local preview
            preview=$(head -3 "$outfile" 2>/dev/null || true)
            echo -e "  ${YELLOW}[DEBUG]${NC} Server response: ${preview}"
        fi
        rm -f "$outfile"
        echo -e "  ${YELLOW}[WARN]${NC}  ${acc} attempt ${attempt}/${RETRIES} failed"
        local backoff=$(( attempt * attempt * 5 + RANDOM % 5 ))
        echo -e "  ${YELLOW}[INFO]${NC}  Retrying in ${backoff}s..."
        sleep "$backoff"
    done

    if [[ "$success" == false ]]; then
        echo -e "  ${RED}[FAIL]${NC}  ${acc} could not be downloaded after ${RETRIES} attempts"
        echo "$acc" >> "${OUTDIR}/failed_accessions.txt"
        return 1
    fi
}

# ── Read accessions ───────────────────────────────────────────────────────────
ACCESSIONS=()
while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    [[ "$line" =~ ^[[:space:]]*$ ]] && continue
    acc=$(echo "$line" | tr -d '[:space:]')
    ACCESSIONS+=("$acc")
done < "$INPUT"

total=${#ACCESSIONS[@]}
if [[ $total -eq 0 ]]; then
    echo -e "${RED}ERROR:${NC} No accessions found in $INPUT"
    exit 1
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  GenBank Sequence Downloader (.gbk)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Input file : $INPUT"
echo "  Output dir : $OUTDIR"
echo "  Accessions : $total"
echo "  Downloader : $DOWNLOADER"
if [[ "$BUILD_BLAST" == true ]]; then
echo "  BLAST dbs  : yes → ${OUTDIR}/blast_db"
else
echo "  BLAST dbs  : no  (use -b to enable)"
fi
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

rm -f "${OUTDIR}/failed_accessions.txt"
success_count=0
fail_count=0
blast_count=0
blast_fail_count=0

for acc in "${ACCESSIONS[@]}"; do
    if download_accession "$acc"; then
        ((success_count++)) || true
        if [[ "$BUILD_BLAST" == true ]]; then
            if make_blast_db "$acc"; then
                ((blast_count++)) || true
            else
                ((blast_fail_count++)) || true
            fi
        fi
    else
        ((fail_count++)) || true
    fi
    sleep "$DELAY"
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "  ${GREEN}Downloaded${NC} : ${success_count}"
echo -e "  ${RED}Failed${NC}     : ${fail_count}"
if [[ "$BUILD_BLAST" == true ]]; then
    echo -e "  ${GREEN}BLAST dbs${NC}  : ${blast_count}"
    if [[ $blast_fail_count -gt 0 ]]; then
        echo -e "  ${RED}BLAST fail${NC} : ${blast_fail_count}"
    fi
    echo "  BLAST db dir: ${OUTDIR}/blast_db"
fi
if [[ $fail_count -gt 0 ]]; then
    echo ""
    echo "  Failed accessions saved to: ${OUTDIR}/failed_accessions.txt"
fi
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
