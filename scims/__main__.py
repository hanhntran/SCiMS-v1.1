"""
SCiMS: Sex Calling in Metagenomic Sequencing

A tool for inferring host chromosomal sex from metagenomic sequencing data.

Subcommands:
    scims call        - Classify sex from .idxstats files or BAM files

Author: Hanh Tran
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

from ._version import __version__
from .utils import read_metadata, find_sample_id_column, extract_sample_id
from .process_input_file import process_input_file

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

BANNER = """
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|  
    _|      _|         _|    _|_|  _|_|   _|        
    _|_|_|  _|         _|    _|  _|  _|   _|_|_|    
        _|  _|         _|    _|      _|       _|  
    _|_|_|   _|_|_|  _|_|_|  _|      _|   _|_|_|    
================================================="""


def setup_logger(output_dir: str, log_to_file: bool = False) -> logging.Logger:
    """Configure and return the SCiMS logger."""
    logger = logging.getLogger("scims")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    if logger.hasHandlers():
        logger.handlers.clear()

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_to_file and output_dir:
        os.makedirs(output_dir, exist_ok=True)
        log_path = os.path.join(output_dir, "scims.log")
        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info(f"Log file created at: {log_path}")

    logger.info(f"\n================================================={BANNER}")
    logger.info("SCiMS: Sex Calling in Metagenomic Sequencing")
    logger.info(f"Version: {__version__}")
    logger.info("=================================================")

    return logger

# ---------------------------------------------------------------------------
# Subcommand: call
# ---------------------------------------------------------------------------

def add_call_parser(subparsers):
    """Add the 'call' subcommand parser."""
    p = subparsers.add_parser(
        "call",
        help="Classify host chromosomal sex from BAM or .idxstats files",
        description=(
            "Classify host chromosomal sex from alignment data. "
            "Uses an analytical multinomial likelihood model by default (no training data needed). "
            "Accepts BAM files directly (requires samtools) or pre-computed .idxstats files."
        ),
    )

    # Input: BAM or idxstats (user provides one)
    input_grp = p.add_argument_group("Input (provide one)")
    input_grp.add_argument("--bam", help="Path to a single BAM file")
    input_grp.add_argument("--bam_folder", help="Path to folder of BAM files (batch mode)")
    input_grp.add_argument("--idxstats_file", help="Path to a single .idxstats file")
    input_grp.add_argument("--idxstats_folder", help="Path to folder of .idxstats files (batch mode)")

    p.add_argument("--scaffolds", dest="scaffold_ids_file", required=True, help="Path to scaffolds.txt")
    p.add_argument("--homogametic_id", dest="homogametic_scaffold", required=True, help="Scaffold ID for homogametic sex chromosome (X or Z)")
    p.add_argument("--heterogametic_id", dest="heterogametic_scaffold", required=True, help="Scaffold ID for heterogametic sex chromosome (Y or W)")
    p.add_argument("--ZW", dest="ZW", action="store_true", help="Use ZW sex determination system (default: XY)")
    p.add_argument("--threshold", type=float, default=0.95, help="Posterior probability threshold (default: 0.95)")
    p.add_argument("--mismapping_rate", type=float, default=1e-4, help="Expected mismapping rate for reads on zero-ploidy chromosomes (default: 1e-4)")
    p.add_argument("--output_dir", required=True, help="Output directory")
    p.add_argument("--metadata", help="Path to metadata file (optional)")
    p.add_argument("--id_column", help="Sample ID column name in metadata")
    p.add_argument("--log", action="store_true", help="Write log file to output directory")
    p.set_defaults(func=run_call)


def _bam_to_idxstats(bam_path: str, output_dir: str, logger) -> str:
    """
    Convert a BAM file to an idxstats file using samtools.
    Indexes the BAM first if needed. Returns the path to the .idxstats file.
    """
    import shutil
    import subprocess

    samtools = shutil.which("samtools")
    if samtools is None:
        raise RuntimeError(
            "samtools not found on PATH. Install it (e.g., 'conda install -c bioconda samtools') "
            "to use BAM files as input, or provide pre-computed .idxstats files instead."
        )

    # Check if BAM index exists, create if not
    bai_path = bam_path + ".bai"
    csi_path = bam_path + ".csi"
    if not os.path.exists(bai_path) and not os.path.exists(csi_path):
        logger.info(f"  Indexing {os.path.basename(bam_path)}...")
        result = subprocess.run(
            [samtools, "index", bam_path],
            capture_output=True,
        )
        if result.returncode != 0:
            stderr = result.stderr.decode("utf-8", errors="replace")
            raise RuntimeError(
                f"samtools index failed for {bam_path}. "
                f"Is the BAM sorted? Error: {stderr}"
            )

    # Generate idxstats
    base_name = os.path.splitext(os.path.basename(bam_path))[0]
    idxstats_path = os.path.join(output_dir, f"{base_name}.idxstats")

    result = subprocess.run(
        [samtools, "idxstats", bam_path],
        capture_output=True,
    )
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace")
        raise RuntimeError(f"samtools idxstats failed for {bam_path}: {stderr}")

    with open(idxstats_path, "wb") as f:
        f.write(result.stdout)

    return idxstats_path


def _resolve_input(args, output_dir, logger) -> tuple:
    """
    Resolve the input source to a list of idxstats file paths.
    Accepts BAM files, BAM folders, idxstats files, or idxstats folders.
    Returns (list_of_input_paths, is_batch).
    """
    if args.bam:
        logger.info(f"[call] Converting BAM to idxstats: {os.path.basename(args.bam)}")
        idxstats_path = _bam_to_idxstats(args.bam, output_dir, logger)
        return [idxstats_path], False

    elif args.bam_folder:
        bam_files = sorted([
            os.path.join(args.bam_folder, f)
            for f in os.listdir(args.bam_folder)
            if f.endswith(".bam")
        ])
        if not bam_files:
            logger.error("No .bam files found in the provided folder.")
            sys.exit(1)
        logger.info(f"[call] Converting {len(bam_files)} BAM files to idxstats...")
        idxstats_files = []
        for bam in bam_files:
            idxstats_files.append(_bam_to_idxstats(bam, output_dir, logger))
        return idxstats_files, True

    elif args.idxstats_file:
        return [args.idxstats_file], False

    elif args.idxstats_folder:
        folder_files = sorted([
            os.path.join(args.idxstats_folder, f)
            for f in os.listdir(args.idxstats_folder)
            if f.endswith(".idxstats")
        ])
        if not folder_files:
            logger.error("No .idxstats files found in the provided folder.")
            sys.exit(1)
        return folder_files, True

    else:
        logger.error(
            "Provide one of: --bam, --bam_folder, --idxstats_file, or --idxstats_folder."
        )
        sys.exit(1)


def run_call(args):
    """Execute the 'call' subcommand."""
    logger = setup_logger(args.output_dir, log_to_file=args.log)
    os.makedirs(args.output_dir, exist_ok=True)

    if args.metadata and not args.id_column:
        logger.error("When providing a metadata file, you must also specify --id_column.")
        sys.exit(1)

    # Resolve inputs (BAM to idxstats conversion happens here if needed)
    try:
        input_files, is_batch = _resolve_input(args, args.output_dir, logger)
    except RuntimeError as e:
        logger.error(str(e))
        sys.exit(1)

    # Load scaffold IDs
    try:
        with open(args.scaffold_ids_file, "r") as f:
            scaffold_ids = [line.strip() for line in f if line.strip()]
    except Exception as e:
        logger.error(f"Failed to read scaffold IDs: {e}")
        sys.exit(1)

    logger.info("[call] Using analytical multinomial likelihood classifier")

    # Process all input files
    all_results = []
    for input_file in input_files:
        result = process_input_file(
            input_file=input_file,
            scaffold_ids=scaffold_ids,
            args=args,
        )
        all_results.append(result)
        out_dict = {
            "SCiMS_ID": result.get("SCiMS_ID"),
            "SCiMS_predicted_sex": result.get("SCiMS_sex"),
            "SCiMS_male_post_prob": result.get("SCiMS_male_post_prob"),
            "SCiMS_female_post_prob": result.get("SCiMS_female_post_prob"),
        }
        base_name = extract_sample_id(os.path.basename(input_file))
        output_file = os.path.join(args.output_dir, f"{base_name}_results.txt")
        pd.DataFrame([out_dict]).to_csv(output_file, sep="\t", index=False)
        logger.info(f"Results written to {output_file}")

    # Merge with metadata if provided
    if args.metadata:
        try:
            results_df = pd.DataFrame(all_results)
            metadata = read_metadata(args.metadata)
            sample_id_col = find_sample_id_column(metadata, args.id_column)
            merged_df = pd.merge(
                metadata, results_df,
                left_on=sample_id_col, right_on="SCiMS_ID", how="left",
            )
            merged_df.drop(columns=["SCiMS_ID"], inplace=True, errors="ignore")
            metadata_basename = os.path.basename(args.metadata).split(".")[0]
            metadata_file = os.path.join(args.output_dir, f"{metadata_basename}_scims_updated.txt")
            merged_df.to_csv(metadata_file, sep="\t", index=False)
            logger.info(f"Updated metadata written to {metadata_file}")
        except Exception as e:
            logger.error(f"Error updating metadata: {e}")
            sys.exit(1)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="scims",
        description="SCiMS: Sex Calling in Metagenomic Sequencing - "
                    "Infer host chromosomal sex from metagenomic data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Classify directly from BAM files:\n"
            "  scims call --bam sample.bam --scaffolds scaffolds.txt \\\n"
            "       --homogametic_id NC_000023.11 --heterogametic_id NC_000024.10 --output_dir out/\n\n"
            "  # Batch mode from a folder of BAMs:\n"
            "  scims call --bam_folder /path/to/bams/ --scaffolds scaffolds.txt \\\n"
            "       --homogametic_id NC_000023.11 --heterogametic_id NC_000024.10 --output_dir out/\n\n"
            "  # From pre-computed idxstats files:\n"
            "  scims call --idxstats_file sample.idxstats --scaffolds scaffolds.txt \\\n"
            "       --homogametic_id NC_000023.11 --heterogametic_id NC_000024.10 --output_dir out/\n\n"
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    add_call_parser(subparsers)

    # Handle no arguments
    argv = sys.argv[1:]

    if not argv:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    main()