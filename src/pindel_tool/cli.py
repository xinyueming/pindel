#!/usr/bin/env python3
"""CLI entry point for pindel_tool.

Subcommands:
    pindel      - Run Pindel to detect structural variants
    pindel2vcf  - Convert Pindel output to standard VCF format
    anno        - ANNOVAR-based gene annotation
    filter      - Filter and format annotated VCF output

Usage:
    pindel-tool <subcommand> [options]
    pindel-tool --help
"""

import argparse
import os
import subprocess
import sys

try:
    from importlib.resources import files as _pkg_files
except ImportError:
    from importlib_resources import files as _pkg_files

from pindel_tool import __version__


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tool_path(*parts):
    """Resolve a packaged resource path."""
    return str(_pkg_files("pindel_tool").joinpath(*parts))


def _run(cmd, label):
    """Run a subprocess, streaming stdout/stderr, and return its exit code."""
    print(f"[{label}] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"[{label}] Failed with exit code {result.returncode}", file=sys.stderr)
    return result.returncode


def _check_file(path, label="file"):
    """Check if a file/directory exists. Returns 1 with error message if not."""
    if not os.path.exists(path):
        print(f"Error: {label} not found: {path}", file=sys.stderr)
        return 1
    return 0


def _check_command(cmd):
    """Check if a command is available."""
    try:
        subprocess.run([cmd, "--version"], capture_output=True)
        return True
    except FileNotFoundError:
        return False


def _get_read_length(bam_path):
    """Get read length from samtools stats."""
    result = subprocess.run(
        ["samtools", "stats", bam_path],
        capture_output=True, text=True,
    )
    for line in result.stdout.splitlines():
        if line.startswith("RL"):
            parts = line.split()
            if len(parts) >= 2:
                return int(parts[1])
    raise RuntimeError("Could not determine read length from BAM")


def _get_insert_size(bam_path):
    """Get median and mean insert size from Picard CollectInsertSizeMetrics."""
    import tempfile

    picard_cmds = [
        ["picard", "CollectInsertSizeMetrics"],
        ["java", "-jar", _tool_path("picard.jar"), "CollectInsertSizeMetrics"],
    ]

    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
        metrics_file = f.name

    try:
        for cmd in picard_cmds:
            try:
                full_cmd = cmd + [
                    f"I={bam_path}",
                    f"O={metrics_file}",
                    f"H={metrics_file}.hist",
                ]
                result = subprocess.run(full_cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    break
            except FileNotFoundError:
                continue
        else:
            raise RuntimeError(
                "Picard not found. Install with: conda install -c bioconda picard"
            )

        with open(metrics_file) as f:
            for line in f:
                if line.startswith("MEDIAN_INSERT_SIZE"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    try:
                        median = int(float(parts[0]))
                        mean = int(float(parts[1]))
                        return median, mean
                    except ValueError:
                        continue
        raise RuntimeError("Could not determine insert size from BAM")
    finally:
        import os

        for f in [metrics_file, metrics_file + ".hist"]:
            if os.path.exists(f):
                os.unlink(f)


def _generate_config(bam_path, sample_name, output_path):
    """Generate pindel config file from BAM."""
    rl = _get_read_length(bam_path)
    median, mean = _get_insert_size(bam_path)
    insert_size = max(median, rl)
    with open(output_path, "w") as f:
        f.write(f"{bam_path}\t{insert_size}\t{sample_name}\n")
    return output_path


def cmd_config(args):
    """Generate pindel config file from BAM."""
    if not _check_command("samtools"):
        print("Error: samtools not found. Install with: conda install -c bioconda samtools", file=sys.stderr)
        return 1
    rc = _check_file(args.bam, "BAM file")
    if rc:
        return rc
    try:
        config_path = _generate_config(args.bam, args.sample, args.output)
        print(f"[config] Wrote config to {config_path}")
    except (RuntimeError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


# ---------------------------------------------------------------------------
# Subcommand implementations
# ---------------------------------------------------------------------------

def cmd_pindel(args):
    """Call the pindel executable to detect structural variants."""
    rc = _check_file(args.fasta, "Reference FASTA")
    if rc:
        return rc

    # Auto-generate config from BAM if --bam provided
    if args.bam:
        if not args.sample:
            print("Error: --sample required with --bam", file=sys.stderr)
            return 1
        if not _check_command("samtools"):
            print("Error: samtools not found. Install with: conda install -c bioconda samtools", file=sys.stderr)
            return 1
        rc = _check_file(args.bam, "BAM file")
        if rc:
            return rc
        args.config_file = _generate_config(args.bam, args.sample, args.output_prefix + ".config")
    elif not args.config_file:
        print("Error: --config-file or --bam required", file=sys.stderr)
        return 1

    rc = _check_file(args.config_file, "Config file")
    if rc:
        return rc
    cmd = [
        _tool_path("pindel"),
        "-f", args.fasta,
        "-i", args.config_file,
        "-o", args.output_prefix,
        "-c", args.chromosome,
        "-T", str(args.number_of_threads),
        "-a", str(args.additional_mismatch),
        "-M", str(args.minimum_support_for_event),
    ]
    return _run(cmd, "pindel")


def cmd_pindel2vcf(args):
    """Call pindel2vcf to convert Pindel output to VCF format."""
    rc = _check_file(args.reference, "Reference FASTA")
    if rc:
        return rc
    cmd = [
        _tool_path("pindel2vcf"),
        "-P", args.pindel_output_root,
        "-r", args.reference,
        "-R", args.reference_name,
        "-d", args.reference_date,
        "-e", str(args.min_supporting_reads),
        "-he", str(args.het_cutoff),
        "-ho", str(args.hom_cutoff),
    ]
    if args.vcf:
        cmd.extend(["-v", args.vcf])
    return _run(cmd, "pindel2vcf")


def cmd_anno(args):
    """Run ANNOVAR table_annovar.pl for gene annotation."""
    rc = _check_file(args.input_vcf, "Input VCF")
    if rc:
        return rc
    rc = _check_file(args.db_path, "ANNOVAR database directory")
    if rc:
        return rc
    out = args.out if args.out else os.path.splitext(os.path.basename(args.input_vcf))[0]
    cmd = [
        "perl",
        _tool_path("annovar_scripts", "table_annovar.pl"),
        args.input_vcf,
        args.db_path,
        "--buildver", args.buildver,
        "--protocol", args.protocol,
        "--operation", args.operation,
        "--nastring", args.nastring,
        "--outfile", out,
        "--vcfinput",
    ]
    return _run(cmd, "anno")


def cmd_filter(args):
    """Filter annotated VCF by gene/transcript and output as TSV table."""
    rc = _check_file(args.input_vcf, "Input VCF")
    if rc:
        return rc
    header_fields = [
        "CHROM", "POS", "ID", "REF", "ALT", "Gene", "Transcript",
        "SVTYPE", "SVLEN", "Insertion", "CDS", "AA",
        "GT", "AD", "VD", "DP", "AF", "AR", "Sample",
    ]

    # Parse filter targets
    target_genes = set(args.gene.split(",")) if args.gene else None
    target_transcripts = set(_strip_version(t) for t in args.transcript.split(",")) if args.transcript else None
    gene_transcript_pairs = {}
    if args.gene_transcript_pair:
        for pair_str in args.gene_transcript_pair:
            for pair in pair_str.split(","):
                pair = pair.strip()
                if not pair:
                    continue
                g, t = pair.split(":", 1)
                gene_transcript_pairs.setdefault(g, set()).add(_strip_version(t))

    # Parse VCF
    records = _parse_vcf_records(
        args.input_vcf,
        target_genes,
        target_transcripts,
        gene_transcript_pairs,
    )

    # Output
    lines = ["\t".join(header_fields)]
    for rec in records:
        lines.append("\t".join(str(rec.get(f, ".")) for f in header_fields))

    out_text = "\n".join(lines) + "\n"
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(out_text)
        print(f"[filter] Wrote {len(records)} record(s) to {args.output}")
    else:
        sys.stdout.write(out_text)
    return 0


# ---------------------------------------------------------------------------
# VCF parsing helpers
# ---------------------------------------------------------------------------

def _parse_info(info_str):
    """Parse VCF INFO field into a dict. Values stay as raw strings."""
    info = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            info[k] = v
        else:
            info[token] = True
    return info


def _parse_aachange(aachange_str):
    """Parse AAChange.refGeneWithVer into list of dicts.

    Format: GENE:TRANSCRIPT:exonN:cDNA:AA, comma-separated.
    Encoded semicolons (\\x3b) decode to ','.
    """
    aachange_str = aachange_str.replace("\\x3b", ",")
    entries = []
    for entry in aachange_str.split(","):
        parts = entry.split(":")
        if len(parts) >= 5:
            entries.append({
                "gene": parts[0],
                "transcript": parts[1],
                "exon": parts[2],
                "cds": parts[3],
                "aa": parts[4],
            })
        elif len(parts) == 2:
            entries.append({
                "gene": parts[0],
                "transcript": parts[1],
                "exon": ".",
                "cds": ".",
                "aa": ".",
            })
    return entries


def _strip_version(transcript):
    """Remove version suffix from transcript ID (e.g. NM_004119.3 → NM_004119)."""
    if not transcript:
        return transcript
    return transcript.split(".")[0] if "." in transcript else transcript


def _matches_filter(gene, transcript, target_genes, target_transcripts, gene_transcript_pairs):
    """Check if a gene/transcript combo matches any of the filter criteria."""
    if target_genes is None and target_transcripts is None and not gene_transcript_pairs:
        return True

    transcript_base = _strip_version(transcript)

    if gene_transcript_pairs:
        if gene in gene_transcript_pairs and transcript_base in gene_transcript_pairs[gene]:
            return True

    if target_genes is not None and gene in target_genes:
        if not target_transcripts:
            return True

    if target_transcripts is not None and transcript_base in target_transcripts:
        if not target_genes:
            return True

    return False


def _parse_vcf_records(vcf_path, target_genes, target_transcripts, gene_transcript_pairs):
    """Read VCF, parse records, filter, and return list of row dicts."""
    records = []
    sample_name = None
    format_keys = None

    with open(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                cols = line.lstrip("#").split("\t")
                if len(cols) > 9:
                    sample_name = cols[9]
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info_str = parts[7]

            info = _parse_info(info_str)
            gene = info.get("Gene.refGeneWithVer", ".")
            aachange_raw = info.get("AAChange.refGeneWithVer", "")
            svtype = info.get("SVTYPE", ".")
            svlen = info.get("SVLEN", ".")

            if len(parts) > 8:
                format_keys = parts[8].split(":")
            sample_data = {}
            if len(parts) > 9 and format_keys:
                sample_vals = parts[9].split(":")
                sample_data = dict(zip(format_keys, sample_vals))

            gt = sample_data.get("GT", ".")
            ad_raw = sample_data.get("AD", ".")

            if ad_raw != ".":
                ad_parts = [int(x) for x in ad_raw.split(",")]
                if len(ad_parts) >= 2:
                    vd = ad_parts[1]
                    dp = sum(ad_parts)
                    af = round(vd / dp, 4) if dp > 0 else 0
                    ar = round(vd / ad_parts[0], 4) if ad_parts[0] > 0 else "."
                    ad = ad_raw
                else:
                    vd = "."
                    dp = "."
                    af = "."
                    ar = "."
                    ad = ad_raw
            else:
                vd = dp = af = ar = "."
                ad = "."

            aachanges = _parse_aachange(aachange_raw) if aachange_raw else []

            if not aachanges:
                if _matches_filter(gene, None, target_genes, target_transcripts, gene_transcript_pairs):
                    insertion = len(alt) - len(ref) if svtype == "INS" else "."
                    records.append({
                        "CHROM": chrom, "POS": pos, "ID": vid or ".",
                        "REF": ref, "ALT": alt, "Gene": gene,
                        "Transcript": ".", "SVTYPE": svtype, "SVLEN": svlen,
                        "Insertion": insertion, "CDS": ".", "AA": ".",
                        "GT": gt, "AD": ad, "VD": vd, "DP": dp, "AF": af, "AR": ar,
                        "Sample": sample_name or ".",
                    })
                continue

            for aa_entry in aachanges:
                a_gene = aa_entry["gene"]
                a_transcript = aa_entry["transcript"]
                a_cds = aa_entry["cds"]
                a_aa = aa_entry["aa"]

                if not _matches_filter(a_gene, a_transcript, target_genes, target_transcripts, gene_transcript_pairs):
                    continue

                insertion = len(alt) - len(ref) if svtype == "INS" else "."

                records.append({
                    "CHROM": chrom, "POS": pos, "ID": vid or ".",
                    "REF": ref, "ALT": alt, "Gene": a_gene,
                    "Transcript": a_transcript, "SVTYPE": svtype,
                    "SVLEN": svlen, "Insertion": insertion,
                    "CDS": a_cds, "AA": a_aa,
                    "GT": gt, "AD": ad, "VD": vd, "DP": dp, "AF": af, "AR": ar,
                    "Sample": sample_name or ".",
                })

    return records


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser():
    """Build and return the top-level argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="pindel-tool",
        description="CLI wrapper for Pindel structural variant detection pipeline.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # -- pindel -----------------------------------------------------------
    p_pindel = subparsers.add_parser(
        "pindel",
        help="Run Pindel to detect structural variants from BAM files",
    )
    p_pindel.add_argument("-f", "--fasta", required=True, help="Reference genome FASTA file")
    p_pindel.add_argument("-i", "--config-file", default=None, help="BAM config file (optional if --bam provided)")
    p_pindel.add_argument("-b", "--bam", default=None, help="BAM file path (auto-generates config)")
    p_pindel.add_argument("-s", "--sample", default=None, help="Sample name (required with --bam)")
    p_pindel.add_argument("-o", "--output-prefix", required=True, help="Output file prefix")
    p_pindel.add_argument("-c", "--chromosome", default="ALL", help="Chromosome region (default: ALL)")
    p_pindel.add_argument("-T", "--number-of-threads", type=int, default=1, help="Number of threads (default: 1)")
    p_pindel.add_argument("-a", "--additional-mismatch", type=int, default=1, help="Additional mismatch tolerance (default: 1)")
    p_pindel.add_argument("-M", "--minimum-support-for-event", type=int, default=1, help="Minimum supporting reads for an event (default: 1)")
    p_pindel.set_defaults(func=cmd_pindel)

    # -- config -----------------------------------------------------------
    p_config = subparsers.add_parser(
        "config",
        help="Generate pindel config file from BAM",
    )
    p_config.add_argument("-b", "--bam", required=True, help="BAM file path")
    p_config.add_argument("-s", "--sample", required=True, help="Sample name")
    p_config.add_argument("-o", "--output", default="pindel.config", help="Output config file (default: pindel.config)")
    p_config.set_defaults(func=cmd_config)

    # -- pindel2vcf -------------------------------------------------------
    p_pindel2vcf = subparsers.add_parser(
        "pindel2vcf",
        help="Convert Pindel output to standard VCF format",
    )
    p_pindel2vcf.add_argument("-P", "--pindel-output-root", required=True, help="Pindel output file prefix")
    p_pindel2vcf.add_argument("-r", "--reference", required=True, help="Reference genome FASTA file")
    p_pindel2vcf.add_argument("-R", "--reference-name", required=True, help="Reference name (e.g. hg38)")
    p_pindel2vcf.add_argument("-d", "--reference-date", required=True, help="Reference version date")
    p_pindel2vcf.add_argument("-v", "--vcf", default=None, help="Output VCF filename")
    p_pindel2vcf.add_argument("-e", "--min-supporting-reads", type=int, default=1, help="Minimum supporting reads (default: 1)")
    p_pindel2vcf.add_argument("-he", "--het-cutoff", type=float, default=0.2, help="Heterozygous cutoff (default: 0.2)")
    p_pindel2vcf.add_argument("-ho", "--hom-cutoff", type=float, default=0.8, help="Homozygous cutoff (default: 0.8)")
    p_pindel2vcf.set_defaults(func=cmd_pindel2vcf)

    # -- anno -------------------------------------------------------------
    p_anno = subparsers.add_parser(
        "anno",
        help="Run ANNOVAR for gene annotation on VCF files",
    )
    p_anno.add_argument("input_vcf", help="Input VCF file (annotated or unannotated)")
    p_anno.add_argument(
        "db_path",
        nargs="?",
        default="/mnt/nas/zhangrs/3.database/hg38/humandb/",
        help="ANNOVAR database path (default: /mnt/nas/zhangrs/3.database/hg38/humandb/)",
    )
    p_anno.add_argument("--buildver", default="hg38", help="Genome build version (default: hg38)")
    p_anno.add_argument("--protocol", default="refGeneWithVer", help="Annotation protocol (default: refGeneWithVer)")
    p_anno.add_argument("--operation", default="g", help="Operation type (default: g)")
    p_anno.add_argument("--nastring", default=".", help="String for missing values (default: .)")
    p_anno.add_argument("--out", default=None, help="Output file prefix (default: based on input filename)")
    p_anno.set_defaults(func=cmd_anno)

    # -- filter -----------------------------------------------------------
    p_filter = subparsers.add_parser(
        "filter",
        help="Filter annotated VCF by gene/transcript and output as table",
    )
    p_filter.add_argument("input_vcf", help="Annotated VCF file (*_multianno.vcf)")
    p_filter.add_argument("--gene", default=None, help="Gene name(s), comma-separated. Multiple filter types use OR logic (e.g. FLT3,KMT2A)")
    p_filter.add_argument("--transcript", default=None, help="Transcript ID(s), comma-separated. Multiple filter types use OR logic")
    p_filter.add_argument("--gene-transcript-pair", action="append", default=None, help="Gene:transcript pair (e.g. FLT3:NM_004119.3, repeatable). Multiple filter types use OR logic")
    p_filter.add_argument("-o", "--output", default=None, help="Output file path (default: stdout)")
    p_filter.set_defaults(func=cmd_filter)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = build_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    rc = args.func(args)
    sys.exit(rc)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)
