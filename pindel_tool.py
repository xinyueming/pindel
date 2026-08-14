#!/usr/bin/env python3
"""pindel_tool.py - CLI wrapper for Pindel structural variant detection pipeline.

Subcommands:
    pindel      - Run Pindel to detect structural variants
    pindel2vcf  - Convert Pindel output to standard VCF format
    anno        - ANNOVAR-based gene annotation
    filter      - Filter and format annotated VCF output

Usage:
    python pindel_tool.py <subcommand> [options]
    python pindel_tool.py --help
"""

import argparse
import os
import subprocess
import sys

__version__ = "0.1.0"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tool_path(*parts):
    """Resolve a path relative to this script's directory."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), *parts)


def _run(cmd, label):
    """Run a subprocess, streaming stdout/stderr, and return its exit code."""
    print(f"[{label}] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"[{label}] Failed with exit code {result.returncode}", file=sys.stderr)
    return result.returncode


# ---------------------------------------------------------------------------
# Subcommand implementations
# ---------------------------------------------------------------------------

def cmd_pindel(args):
    """Call the pindel executable to detect structural variants."""
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
    """Filter annotated VCF by gene/transcript and output table format."""
    # TODO: implement VCF parsing and filtering
    print("[filter] Not yet implemented")
    return 0


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser():
    """Build and return the top-level argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="pindel_tool",
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
    p_pindel.add_argument("-i", "--config-file", required=True, help="BAM config file (bam_path insert_size sample_tag)")
    p_pindel.add_argument("-o", "--output-prefix", required=True, help="Output file prefix")
    p_pindel.add_argument("-c", "--chromosome", default="ALL", help="Chromosome region (default: ALL)")
    p_pindel.add_argument("-T", "--number-of-threads", type=int, default=1, help="Number of threads (default: 1)")
    p_pindel.add_argument("-a", "--additional-mismatch", type=int, default=1, help="Additional mismatch tolerance (default: 1)")
    p_pindel.add_argument("-M", "--minimum-support-for-event", type=int, default=1, help="Minimum supporting reads for an event (default: 1)")
    p_pindel.set_defaults(func=cmd_pindel)

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
    p_filter.add_argument("--gene", default=None, help="Gene name(s), comma-separated (e.g. FLT3,KMT2A)")
    p_filter.add_argument("--transcript", default=None, help="Transcript ID(s), comma-separated")
    p_filter.add_argument("--gene-transcript-pair", action="append", default=None, help="Gene:transcript pair (e.g. FLT3:NM_004119.3, repeatable)")
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
    main()
