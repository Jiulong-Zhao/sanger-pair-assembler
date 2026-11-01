"""Command line entry point and reproducible result files."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

from . import __version__
from .core import assemble, read_sequence, reverse_complement, to_fasta

CONFLICT_COLUMNS = [
    "consensus_pos_1based", "F_trimmed_pos_1based", "R_revcomp_pos_1based",
    "F_raw_pos_1based", "R_raw_pos_1based", "F_base", "R_base",
    "F_phred", "R_phred", "consensus", "rule",
]


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="sanger-merge",
        description="Assemble forward and original-orientation reverse Sanger reads, reference-free.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--forward",
                        help="Forward AB1/ABI/FASTA/FASTQ/SEQ/TXT")
    parser.add_argument("-r", "--reverse",
                        help="Reverse read in ORIGINAL sequencer orientation")
    parser.add_argument("-o", "--output", default=None,
                        help="Single-pair consensus FASTA path")
    parser.add_argument("-i", "--input-dir", default=None,
                        help="Batch mode: folder of multiple F/R sequence files")
    parser.add_argument("--output-dir", default=None,
                        help="Batch result directory (default: <input-dir>/assembled_results)")
    parser.add_argument("--recursive", action="store_true",
                        help="Batch mode: include files in nested subfolders")
    parser.add_argument("--trim-cutoff", type=float, default=0.05,
                        help="Modified Mott error-probability cutoff (AB1 only); "
                             "e.g. 0.05 ~ Q13, 0.01 = Q20")
    parser.add_argument("--no-trim", action="store_true",
                        help="Do not trim AB1 read ends")
    parser.add_argument("--min-overlap", type=int, default=40,
                        help="Minimum aligned overlap columns")
    parser.add_argument("--min-identity", type=float, default=0.90,
                        help="Minimum fraction of exact A/C/G/T matches in the overlap")
    parser.add_argument("--min-q", type=int, default=20,
                        help="Minimum Phred score to resolve a conflicting base by quality")
    parser.add_argument("--quality-delta", type=int, default=5,
                        help="Minimum Phred advantage to resolve a conflicting base")
    parser.add_argument("--max-r-prefix-clip", type=int, default=15,
                        help="Max unaligned bases before R* overlap (orientation check)")
    parser.add_argument("--max-both-suffix", type=int, default=15,
                        help="Max unaligned bases in BOTH read suffixes")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing output files")
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    args = parser.parse_args(argv)
    if not 0 < args.trim_cutoff < 1:
        parser.error("--trim-cutoff must be strictly between 0 and 1")
    if args.min_overlap < 1 or not 0 < args.min_identity <= 1:
        parser.error("--min-overlap >= 1 and --min-identity in (0, 1] are required")
    if args.min_q < 0 or args.quality_delta < 0:
        parser.error("--min-q and --quality-delta cannot be negative")
    if args.max_r_prefix_clip < 0 or args.max_both_suffix < 0:
        parser.error("Terminal clip limits cannot be negative")
    if args.input_dir is not None:
        if args.forward or args.reverse or args.output:
            parser.error("Use --input-dir for batch mode OR -f/-r/-o for one pair, not both")
        if args.output_dir is None:
            args.output_dir = str(Path(args.input_dir).expanduser() / 'assembled_results')
    else:
        if not args.forward or not args.reverse:
            parser.error("Single mode requires both -f and -r; batch mode uses --input-dir")
        if args.output_dir is not None or args.recursive:
            parser.error("--output-dir and --recursive are only valid with --input-dir")
        if args.output is None:
            args.output = 'merged.fasta'
        if Path(args.forward).expanduser().resolve() == Path(args.reverse).expanduser().resolve():
            parser.error("Forward and reverse inputs must be different files")
    return args


def output_paths(destination: str) -> Dict[str, Path]:
    out = Path(destination).expanduser()
    if out.suffix.lower() not in (".fasta", ".fa", ".fna", ".fas"):
        raise ValueError("Output filename must end in .fasta, .fa, .fna or .fas")
    stem = out.with_suffix("")
    return {
        "fasta": out,
        "forward": Path(str(stem) + ".F_trimmed.fasta"),
        "reverse_rc": Path(str(stem) + ".R_revcomp.fasta"),
        "alignment": Path(str(stem) + ".alignment.txt"),
        "disagreements": Path(str(stem) + ".disagreements.tsv"),
        "summary": Path(str(stem) + ".summary.json"),
    }


def write_alignment(path: Path, alignment) -> None:
    f, marks, r = alignment
    with path.open("w", encoding="utf-8") as handle:
        handle.write("F vs reverse-complemented R; aligned overlap only\n")
        handle.write("| = identical character; . = mismatch; blank = gap\n\n")
        for start in range(0, len(f), 80):
            handle.write("F  {}\n".format(f[start:start + 80]))
            handle.write("   {}\n".format(marks[start:start + 80]))
            handle.write("R* {}\n\n".format(r[start:start + 80]))


def run(args: argparse.Namespace) -> int:
    paths = output_paths(args.output)
    inputs = {Path(args.forward).expanduser().resolve(), Path(args.reverse).expanduser().resolve()}
    if any(path.expanduser().resolve() in inputs for path in paths.values()):
        raise ValueError('Output path must not overwrite an input sequence file')
    existing = [str(p) for p in paths.values() if p.exists()]
    if existing and not args.force:
        raise ValueError("Output already exists: {} (use --force to overwrite)".format(
            ", ".join(existing)))
    if args.force:
        # Prevent a previous consensus surviving when a later run fails QC.
        for path in paths.values():
            if path.is_file():
                path.unlink()

    forward = read_sequence(args.forward, trim=not args.no_trim, cutoff=args.trim_cutoff)
    reverse = read_sequence(args.reverse, trim=not args.no_trim, cutoff=args.trim_cutoff)
    reverse_rc = reverse_complement(reverse)
    print("F: {} nt (raw {}, removed 5'={}, 3'={})".format(
        forward.retained_length, forward.raw_length,
        forward.trim_5prime, forward.trim_3prime))
    print("R: {} nt (raw {}, removed 5'={}, 3'={}); reverse-complemented".format(
        reverse.retained_length, reverse.raw_length,
        reverse.trim_5prime, reverse.trim_3prime))

    summary = {
        "software": "sanger-pair-assembler",
        "version": __version__,
        "status": "failed",
        "settings": {
            "trim_cutoff": args.trim_cutoff,
            "no_trim": args.no_trim,
            "min_overlap": args.min_overlap,
            "min_identity": args.min_identity,
            "min_q": args.min_q,
            "quality_delta": args.quality_delta,
            "max_r_prefix_clip": args.max_r_prefix_clip,
            "max_both_suffix": args.max_both_suffix,
        },
        "forward": forward.report(),
        "reverse": reverse.report(),
        "warnings": [
            "No automatic primer removal is performed.",
            "Gaps in the overlap must be checked against chromatograms before publication.",
        ],
    }
    paths["fasta"].parent.mkdir(parents=True, exist_ok=True)
    paths["forward"].write_text(to_fasta("F_trimmed", forward.sequence), encoding="utf-8")
    paths["reverse_rc"].write_text(
        to_fasta("R_reverse_complement", reverse_rc.sequence), encoding="utf-8")
    try:
        result = assemble(forward, reverse_rc, args.min_overlap, args.min_identity,
                          args.min_q, args.quality_delta,
                          args.max_r_prefix_clip, args.max_both_suffix)
    except ValueError as exc:
        summary["error"] = str(exc)
        paths["summary"].write_text(json.dumps(summary, indent=2, ensure_ascii=False)
                                    + "\n", encoding="utf-8")
        print("ERROR: no reliable assembly: {}".format(exc), file=sys.stderr)
        print("Saved separate oriented reads and failure summary.", file=sys.stderr)
        return 2

    conflicts = result.pop("conflicts")
    alignment = result.pop("alignment")
    summary["status"] = "assembled"
    summary["assembly"] = {
        "consensus_length": len(result["sequence"]),
        "overlap_columns": result["overlap_columns"],
        "identity": result["identity"],
        "matches": result["matches"],
        "mismatches": result["mismatches"],
        "internal_gap_bases": result["gaps"],
        "ambiguous_overlap_columns": result["ambiguous_overlap_columns"],
        "difference_rows": len(conflicts),
        "unaligned_R_prefix_bases": result["unaligned_R_prefix_bases"],
        "F_suffix_bases": result["F_suffix_bases"],
        "R_suffix_bases": result["R_suffix_bases"],
        "tail_source": result["tail_source"],
        "alignment_score": result["alignment_score"],
    }
    paths["fasta"].write_text(to_fasta("F_R_consensus", result["sequence"]),
                              encoding="utf-8")
    write_alignment(paths["alignment"], alignment)
    with paths["disagreements"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONFLICT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(conflicts)
    paths["summary"].write_text(json.dumps(summary, indent=2, ensure_ascii=False)
                                + "\n", encoding="utf-8")
    print("Merged: {} nt".format(len(result["sequence"])))
    print("Overlap: {} columns; identity: {:.2%}".format(
        result["overlap_columns"], result["identity"]))
    print("Mismatches: {}; internal gap bases: {}; ambiguous columns: {}".format(
        result["mismatches"], result["gaps"], result["ambiguous_overlap_columns"]))
    print("Outputs:")
    for key, path in paths.items():
        print("  {}: {}".format(key, path))
    if conflicts:
        print("NOTE: inspect {} disagreement/gap rows against raw chromatograms".format(
            len(conflicts)), file=sys.stderr)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    try:
        args = parse_args(argv)
        if args.input_dir is not None:
            from .batch import run_batch
            return run_batch(args)
        return run(args)
    except (ValueError, OSError, IndexError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2
