#!/usr/bin/env python3
"""
Extract nucleotide sequences from one or more GenBank files.

Modes
-----
  Gene mode   (-g)      : search by gene name / locus_tag / product annotation
  BLAST mode  (--blast) : extract subject regions from an outfmt-6 BLAST table

Input (-i) can be:
  • a single .gbk / .gb file
  • a directory  → all *.gbk / *.gb files inside are processed

Output
------
  Single file  : all hits pooled into one FASTA (default when -i is a file)
  Per-genome   : one FASTA per input GenBank (--per-genome, auto-enabled for dirs)
"""

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument(
    "-i", "--input",
    required=True,
    help="GenBank file (.gbk/.gb) OR a directory containing GenBank files",
)
parser.add_argument(
    "-o", "--outfile",
    default=None,
    help=(
        "Output FASTA file when pooling all hits (single-file mode). "
        "Ignored when --per-genome is active. "
        "Default: <gene>.fna or blast_regions.fna"
    ),
)
parser.add_argument(
    "--per-genome",
    action="store_true",
    default=False,
    help=(
        "Write a separate FASTA for each input GenBank file "
        "(auto-enabled when -i is a directory). "
        "Files are written to --outdir."
    ),
)
parser.add_argument(
    "--outdir",
    default="extracted_seqs",
    metavar="DIR",
    help="Output directory for per-genome FASTA files (default: extracted_seqs/)",
)

# --- Gene mode ---
gene_group = parser.add_argument_group("Gene name mode")
gene_group.add_argument(
    "-g", "--gene",
    default=None,
    help="Gene name to search for (case-insensitive); matched against "
         "/gene, /locus_tag, and /product qualifiers",
)
gene_group.add_argument(
    "--feature-type",
    default=None,
    help="Restrict to a specific feature type, e.g. CDS, gene, rRNA (default: any)",
)

# --- BLAST mode ---
blast_group = parser.add_argument_group("BLAST region mode")
blast_group.add_argument(
    "--blast",
    default=None,
    metavar="BLAST_FILE",
    help="BLAST -outfmt 6 tabular file; extracts sbjct (subject) regions",
)
blast_group.add_argument(
    "--padding",
    type=int,
    default=0,
    metavar="BP",
    help="Extra bases to include on each side of each BLAST hit (default: 0)",
)
blast_group.add_argument(
    "--min-pident",
    type=float,
    default=0.0,
    metavar="FLOAT",
    help="Minimum percent identity to keep a BLAST hit (default: 0)",
)
blast_group.add_argument(
    "--min-length",
    type=int,
    default=0,
    metavar="INT",
    help="Minimum alignment length to keep a BLAST hit (default: 0)",
)

args = parser.parse_args()

# ---------------------------------------------------------------------------
# Validate mode
# ---------------------------------------------------------------------------
if args.gene is None and args.blast is None:
    sys.exit("ERROR: Supply either -g/--gene (gene mode) or --blast (BLAST region mode).")
if args.gene and args.blast:
    sys.exit("ERROR: --gene and --blast are mutually exclusive.")

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    sys.exit("ERROR: Biopython is required.\nInstall it with:  pip install biopython")

# ---------------------------------------------------------------------------
# Resolve input files
# ---------------------------------------------------------------------------
GBK_SUFFIXES = {".gbk", ".gb", ".genbank"}

input_path = Path(args.input)
if not input_path.exists():
    sys.exit(f"ERROR: Input path not found: {input_path}")

if input_path.is_dir():
    gbk_files = sorted(
        p for p in input_path.iterdir()
        if p.suffix.lower() in GBK_SUFFIXES
    )
    if not gbk_files:
        sys.exit(f"ERROR: No .gbk/.gb files found in directory: {input_path}")
    print(f"Found {len(gbk_files)} GenBank file(s) in '{input_path}'")
    # Directories always use per-genome output
    args.per_genome = True
else:
    gbk_files = [input_path]

# ---------------------------------------------------------------------------
# Set up output
# ---------------------------------------------------------------------------
if args.per_genome:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Per-genome mode: output directory -> {outdir}/")
else:
    if args.outfile:
        pooled_outfile = Path(args.outfile)
    elif args.gene:
        pooled_outfile = Path(f"{args.gene}.fna")
    else:
        pooled_outfile = Path("blast_regions.fna")

# ---------------------------------------------------------------------------
# BLAST file parsing (done once, shared across all GenBank files)
# ---------------------------------------------------------------------------
QSEQID, SSEQID, PIDENT, LENGTH = 0, 1, 2, 3
SSTART, SEND = 8, 9

blast_by_sseqid = {}
if args.blast:
    from collections import defaultdict
    blast_path = Path(args.blast)
    if not blast_path.exists():
        sys.exit(f"ERROR: BLAST file not found: {blast_path}")

    print(f"Reading BLAST file: {blast_path}")
    blast_hits = []
    with open(blast_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                print(f"  WARNING: Skipping malformed line: {line}")
                continue
            if float(cols[PIDENT]) < args.min_pident:
                continue
            if int(cols[LENGTH]) < args.min_length:
                continue
            blast_hits.append(cols)

    if not blast_hits:
        sys.exit("No BLAST hits passed the filters.")

    print(f"  {len(blast_hits)} hit(s) retained after filtering")

    # Build lookup: sseqid -> [list of hit rows]
    blast_by_sseqid = defaultdict(list)
    for cols in blast_hits:
        blast_by_sseqid[cols[SSEQID]].append(cols)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def extract_gene_hits(gbk_path, gene_query, feature_type_filter):
    """Return list of SeqRecord hits for gene mode."""
    from Bio.Seq import UndefinedSequenceError
    hits = []
    seen_coords = set()  # deduplicate by (record_id, start, end, strand)

    for record in SeqIO.parse(str(gbk_path), "genbank"):

        # Check upfront whether this record has a defined sequence
        try:
            _ = str(record.seq)
            record_has_seq = True
        except UndefinedSequenceError:
            record_has_seq = False

        if not record_has_seq:
            print(f"  WARNING: Record '{record.id}' has no sequence data — skipping")
            continue

        for feature in record.features:
            if feature_type_filter and feature.type.lower() != feature_type_filter.lower():
                continue

            matched = any(
                gene_query in v.lower()
                for qual in ["gene", "locus_tag", "product"]
                for v in feature.qualifiers.get(qual, [])
            )
            if not matched:
                continue

            strand    = "+" if feature.location.strand >= 0 else "-"
            coords    = f"{int(feature.location.start)+1}..{int(feature.location.end)}"
            dedup_key = (record.id, coords, strand)
            if dedup_key in seen_coords:
                continue
            seen_coords.add(dedup_key)

            try:
                nt_seq = feature.extract(record.seq)
            except Exception as e:
                print(f"  WARNING: Could not extract sequence: {e}")
                continue

            gene_val  = feature.qualifiers.get("gene",      ["?"])[0]
            locus_val = feature.qualifiers.get("locus_tag", ["?"])[0]
            product   = feature.qualifiers.get("product",   ["?"])[0]

            seq_id = f"{record.id}|{gene_val}|{locus_val}"
            description = (
                f"product={product} "
                f"location={coords}({strand}) "
                f"type={feature.type} "
                f"source={gbk_path.name}"
            )
            hit = SeqRecord(nt_seq, id=seq_id, description=description)
            hits.append(hit)
            print(f"    [HIT] {seq_id}  {coords}({strand})  {product}")
    return hits


def extract_blast_hits(gbk_path, blast_by_sseqid, padding):
    """Return list of SeqRecord hits for BLAST region mode."""
    hits = []

    # Index all records in this GenBank for fast contig lookup
    gbk_index = SeqIO.to_dict(SeqIO.parse(str(gbk_path), "genbank"))
    if not gbk_index:
        return hits

    # Build a secondary lookup keyed by version-stripped ID (e.g. "NZ_JABEBQ010000003"
    # for a record whose full ID is "NZ_JABEBQ010000003.1"). This lets BLAST hits
    # that lack the version suffix still resolve correctly.
    stripped_index = {rid.rsplit(".", 1)[0]: rec
                      for rid, rec in gbk_index.items()}

    seen_regions = set()  # deduplicate by (sseqid, start_0, end_0, strand)

    for sseqid, hit_rows in blast_by_sseqid.items():
        if sseqid in gbk_index:
            record = gbk_index[sseqid]
        elif sseqid in stripped_index:
            record = stripped_index[sseqid]
        else:
            continue  # contig not in this genome — fine for multi-genome runs
        contig_len = len(record.seq)

        for cols in hit_rows:
            qseqid  = cols[QSEQID]
            pident  = float(cols[PIDENT])
            aln_len = int(cols[LENGTH])
            sstart  = int(cols[SSTART])
            send    = int(cols[SEND])

            start_1based = min(sstart, send)
            end_1based   = max(sstart, send)
            strand       = "+" if sstart <= send else "-"

            start_0 = max(0,          start_1based - 1 - padding)
            end_0   = min(contig_len, end_1based        + padding)

            dedup_key = (sseqid, start_0, end_0, strand)
            if dedup_key in seen_regions:
                print(f"    [SKIP] Duplicate region {sseqid} {start_0+1}..{end_0}({strand}) — skipping")
                continue
            seen_regions.add(dedup_key)

            nt_seq = record.seq[start_0:end_0]
            if strand == "-":
                nt_seq = nt_seq.reverse_complement()

            coords = f"{start_0+1}..{end_0}"
            seq_id = f"{sseqid}|{coords}({strand})|query={qseqid}"
            description = (
                f"pident={pident} "
                f"aln_len={aln_len} "
                f"padding={padding} "
                f"source={gbk_path.name}"
            )
            hit = SeqRecord(nt_seq, id=seq_id, description=description)
            hits.append(hit)
            print(f"    [HIT] {seq_id}  pident={pident}  aln={aln_len}bp")

    return hits


# ---------------------------------------------------------------------------
# Main loop — iterate over all GenBank files
# ---------------------------------------------------------------------------
all_hits   = []
total_hits = 0

for gbk_path in gbk_files:
    print(f"\nProcessing: {gbk_path.name}")

    if args.gene:
        file_hits = extract_gene_hits(gbk_path, args.gene.strip().lower(), args.feature_type)
    else:
        file_hits = extract_blast_hits(gbk_path, blast_by_sseqid, args.padding)

    if not file_hits:
        print(f"  (no hits)")
    else:
        print(f"  {len(file_hits)} hit(s) found")
        total_hits += len(file_hits)

    if args.per_genome:
        if file_hits:
            out_path = Path(args.outdir) / f"{gbk_path.stem}.fna"
            SeqIO.write(file_hits, str(out_path), "fasta")
            print(f"  -> Written to {out_path}")
    else:
        all_hits.extend(file_hits)

# ---------------------------------------------------------------------------
# Write pooled output (single-file mode only)
# ---------------------------------------------------------------------------
if not args.per_genome:
    if not all_hits:
        if args.gene:
            sys.exit(
                f"\nNo features found matching gene '{args.gene}'.\n"
                "Tip: check spelling, or try --feature-type CDS"
            )
        else:
            sys.exit(
                "\nNo sequences extracted — check that contig names in the BLAST "
                "output match record IDs in the GenBank file."
            )
    SeqIO.write(all_hits, str(pooled_outfile), "fasta")
    print(f"\n{len(all_hits)} sequence(s) written to {pooled_outfile}")
else:
    print(f"\nDone. {total_hits} total hit(s) across {len(gbk_files)} file(s).")
    if total_hits == 0:
        if args.gene:
            print(f"Tip: no features matched '{args.gene}'. Check spelling or try --feature-type CDS.")
        else:
            print("Tip: check that contig names in the BLAST output match record IDs in your GenBank files.")
