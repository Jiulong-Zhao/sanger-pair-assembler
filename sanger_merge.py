#!/usr/bin/env python3
"""Sanger Pair Assembler: standalone, batch + single-pair version.

Generated from src/sanger_assembler/{core,cli,batch}.py.
Use only Biopython as a dependency. MIT license: see LICENSE in the project.
"""
from __future__ import annotations
__version__ = "0.2.1"


"""Quality-aware, reference-free assembly of two Sanger sequencing reads.

The default trimmer reproduces Biopython's abi-trim (Modified Mott) slicing
rule at a *configurable* error-probability cutoff. All coordinates in internal
functions are 0-based; human-readable reports use 1-based coordinates.
"""

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IUPAC_SETS = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("CG"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"), "B": set("CGT"),
    "D": set("AGT"), "H": set("ACT"), "V": set("ACG"),
    "N": set("ACGT"),
}
SETS_IUPAC = {frozenset(v): k for k, v in IUPAC_SETS.items()}


@dataclass
class Read:
    sequence: str
    qualities: Optional[List[int]]
    input_file: str
    raw_length: int
    trim_5prime: int = 0
    trim_3prime: int = 0
    trimmed: bool = False

    @property
    def retained_length(self) -> int:
        return len(self.sequence)

    def report(self) -> Dict[str, object]:
        return {
            "input_file": self.input_file,
            "raw_length": self.raw_length,
            "trimmed_5prime_bases": self.trim_5prime,
            "trimmed_3prime_bases": self.trim_3prime,
            "retained_length": self.retained_length,
            "retained_raw_coordinates_1based": [
                self.trim_5prime + 1,
                self.raw_length - self.trim_3prime,
            ],
            "quality_available": self.qualities is not None,
            "quality_trim_applied": self.trimmed,
        }


def mott_trim_bounds(qualities: List[int], cutoff: float = 0.05) -> Tuple[int, int]:
    """Return [start, end) following Biopython's legacy abi-trim convention.

    The default 0.05 results match SeqIO.read(path, 'abi-trim') on AB1 reads.
    Based on segment scores `cutoff - 10**(-Q/10)`; this does NOT remove
    every base below a given Q threshold. Like abi-trim, the first base is
    ignored when accumulating scores; the max-score index is exclusive.
    """
    if not 0 < cutoff < 1:
        raise ValueError("Trim cutoff must be a probability strictly between 0 and 1")
    n = len(qualities)
    if n <= 20:
        return 0, n
    scores = [cutoff - 10.0 ** (-q / 10.0) for q in qualities]
    cumulative = [0.0]
    start = None
    for i in range(1, n):
        value = cumulative[-1] + scores[i]
        if value < 0:
            cumulative.append(0.0)
        else:
            cumulative.append(value)
            if start is None:
                start = i
    if start is None:
        raise ValueError("No high-quality stretch found at the chosen trim cutoff")
    end = cumulative.index(max(cumulative))
    if end <= start:
        raise ValueError("Quality trimming left no usable sequence")
    return start, end


def read_sequence(filename: str, trim: bool = True, cutoff: float = 0.05) -> Read:
    path = Path(filename).expanduser()
    if not path.is_file():
        raise ValueError("Input file not found: {}".format(path))
    ext = path.suffix.lower()
    if ext in (".ab1", ".abi"):
        # Always read the full ABI file so trimming boundaries are reportable.
        record = SeqIO.read(str(path), "abi")
    elif ext in (".fa", ".fasta", ".fna", ".fas"):
        record = SeqIO.read(str(path), "fasta")
    elif ext in (".fq", ".fastq"):
        record = SeqIO.read(str(path), "fastq")
    elif ext in (".seq", ".txt"):
        content = path.read_text(encoding="utf-8-sig").strip()
        if content.startswith(">"):
            record = SeqIO.read(StringIO(content), "fasta")
        else:
            record = SeqRecord(Seq("".join(content.split())), id=path.stem)
    else:
        raise ValueError("Unsupported extension: {}".format(ext))

    original = str(record.seq).upper()
    if not original:
        raise ValueError("Empty sequence: {}".format(path))
    invalid = set(original).difference(IUPAC_SETS)
    if invalid:
        raise ValueError("Invalid DNA characters in {}: {}".format(
            path, "".join(sorted(invalid))))
    quals = record.letter_annotations.get("phred_quality")
    q = list(quals) if quals is not None else None
    if q is not None and len(q) != len(original):
        raise ValueError("Phred quality length disagrees with read length")
    left, right = 0, len(original)
    applied = trim and ext in (".ab1", ".abi")
    if applied:
        if q is None:
            raise ValueError("AB1 file does not include Phred quality scores")
        left, right = mott_trim_bounds(q, cutoff)
    trimmed = original[left:right]
    if not trimmed:
        raise ValueError("No usable sequence after trimming: {}".format(path))
    return Read(trimmed, q[left:right] if q is not None else None, str(path),
                len(original), left, len(original) - right, applied)


def reverse_complement(read: Read) -> Read:
    """Reverse-complement a read, reversing qualities in the same order.

    Metadata remain in original ABI orientation, not R* orientation.
    """
    return Read(str(Seq(read.sequence).reverse_complement()),
                list(reversed(read.qualities)) if read.qualities is not None else None,
                read.input_file, read.raw_length,
                read.trim_5prime, read.trim_3prime, read.trimmed)


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -4.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0
    # End gap penalties removed to permit a dovetail (semi-global) alignment.
    # This property is available in Biopython 1.83 and 1.86.
    aligner.end_gap_score = 0.0
    return aligner


def consensus_base(a: str, b: str, qa: Optional[int], qb: Optional[int],
                   min_q: int = 20, quality_delta: int = 5) -> Tuple[str, str]:
    """Conservative base-call resolution; never invent a base at a conflict."""
    if a == b:
        return a, "agree"
    if a == "N":
        return b, "F_is_N"
    if b == "N":
        return a, "R_is_N"
    if qa is not None and qb is not None:
        if qa >= min_q and qa - qb >= quality_delta:
            return a, "higher_F_quality"
        if qb >= min_q and qb - qa >= quality_delta:
            return b, "higher_R_quality"
    sa, sb = IUPAC_SETS[a], IUPAC_SETS[b]
    if sa.issubset(sb):
        return a, "IUPAC_compatible"
    if sb.issubset(sa):
        return b, "IUPAC_compatible"
    return SETS_IUPAC.get(frozenset(sa | sb), "N"), "ambiguous"


def phred_at(qualities: Optional[List[int]], index: int) -> Optional[int]:
    return qualities[index] if qualities is not None else None


def assemble(forward: Read, reverse_rc: Read, min_overlap: int = 40,
             min_identity: float = 0.90, min_q: int = 20,
             quality_delta: int = 5, max_r_prefix_clip: int = 15,
             max_both_suffix: int = 15) -> Dict[str, object]:
    """Align and merge F / reverse-complement(R); fail closed on weak overlaps."""
    f, r = forward.sequence, reverse_rc.sequence
    align = make_aligner().align(f, r)[0]
    f_indices, r_indices = align.indices
    paired = [k for k, (i, j) in enumerate(zip(f_indices, r_indices))
              if i >= 0 and j >= 0]
    if not paired:
        raise ValueError("No aligned bases between forward and reverse reads")
    first, last = paired[0], paired[-1]
    first_f, first_r = int(f_indices[first]), int(r_indices[first])
    last_f, last_r = int(f_indices[last]), int(r_indices[last])
    if first_r > max_r_prefix_clip:
        raise ValueError("Invalid read order: {} unaligned bases at beginning "
                         "of R*. Check F/R assignments and sequence direction."
                         .format(first_r))
    f_tail, r_tail = len(f) - last_f - 1, len(r) - last_r - 1
    if f_tail > max_both_suffix and r_tail > max_both_suffix:
        raise ValueError("Internal match, not a trustworthy terminal overlap "
                         "(F tail {}, R* tail {})".format(f_tail, r_tail))

    matches = mismatches = gaps = ambiguous = 0
    for k in range(first, last + 1):
        i, j = int(f_indices[k]), int(r_indices[k])
        if i < 0 or j < 0:
            gaps += 1
        elif f[i] in "ACGT" and r[j] in "ACGT":
            if f[i] == r[j]:
                matches += 1
            else:
                mismatches += 1
        else:
            ambiguous += 1
    overlap_columns = last - first + 1
    identity = matches / overlap_columns if overlap_columns else 0.0
    if overlap_columns < min_overlap or identity < min_identity:
        raise ValueError("Overlap QC failed: {} columns (minimum {}), "
                         "identity {:.1%} (minimum {:.1%})".format(
                             overlap_columns, min_overlap, identity, min_identity))

    assembled = list(f[:first_f])
    conflicts = []
    aligned_f, aligned_r, align_marks = [], [], []
    for k in range(first, last + 1):
        i, j = int(f_indices[k]), int(r_indices[k])
        bf = f[i] if i >= 0 else "-"
        br = r[j] if j >= 0 else "-"
        fq = phred_at(forward.qualities, i) if i >= 0 else None
        rq = phred_at(reverse_rc.qualities, j) if j >= 0 else None
        if i < 0:
            base, rule = br, "gap_in_F_review_trace"
        elif j < 0:
            base, rule = bf, "gap_in_R_review_trace"
        else:
            base, rule = consensus_base(bf, br, fq, rq, min_q, quality_delta)
        assembled.append(base)
        aligned_f.append(bf)
        aligned_r.append(br)
        align_marks.append("|" if bf == br else "." if i >= 0 and j >= 0 else " ")
        if bf != br:
            # R* coordinate j maps to original R 1-based coordinate
            # (raw_length - original 3'-trim - j).
            conflicts.append({
                "consensus_pos_1based": len(assembled),
                "F_trimmed_pos_1based": i + 1 if i >= 0 else "",
                "R_revcomp_pos_1based": j + 1 if j >= 0 else "",
                "F_raw_pos_1based": forward.trim_5prime + i + 1 if i >= 0 else "",
                "R_raw_pos_1based": reverse_rc.raw_length - reverse_rc.trim_3prime - j
                if j >= 0 else "",
                "F_base": bf, "R_base": br,
                "F_phred": fq if fq is not None else "",
                "R_phred": rq if rq is not None else "",
                "consensus": base, "rule": rule,
            })

    if r_tail >= f_tail:
        assembled.extend(r[last_r + 1:])
        tail_source = "reverse_complemented_R"
    else:
        assembled.extend(f[last_f + 1:])
        tail_source = "F"
    return {
        "sequence": "".join(assembled),
        "overlap_columns": overlap_columns,
        "identity": identity,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "ambiguous_overlap_columns": ambiguous,
        "conflicts": conflicts,
        "alignment": ("".join(aligned_f), "".join(align_marks), "".join(aligned_r)),
        "unaligned_R_prefix_bases": first_r,
        "F_suffix_bases": f_tail,
        "R_suffix_bases": r_tail,
        "tail_source": tail_source,
        "alignment_score": float(align.score),
    }


def to_fasta(name: str, sequence: str, width: int = 80) -> str:
    return ">{}\n".format(name) + "\n".join(
        sequence[k:k + width] for k in range(0, len(sequence), width)) + "\n"


"""Command line entry point and reproducible result files."""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional


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
            return run_batch(args)
        return run(args)
    except (ValueError, OSError, IndexError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 2


"""Batch pairing and assembly with fail-closed sample matching.

Recognizes typical Sanger filename suffixes like sample_F, sample_R,
sample_27F, sample_1492R and sample_(S1)_[27F]. Files that cannot be
unambiguously paired are NEVER guessed or placed in combined outputs.
"""

import contextlib
import csv
import io
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO


EXTENSIONS = {'.ab1', '.abi', '.fasta', '.fa', '.fna', '.fas', '.fastq', '.fq', '.seq', '.txt'}
# Recognize only a terminal, delimited label (not an arbitrary F/R inside a name).
LABEL_PATTERN = re.compile(
    r'(?i)(?:\[(27f|1492r|forward|reverse|f|r)\]'
    r'|\((27f|1492r|forward|reverse|f|r)\)'
    r'|(?:^|[._\-\s])(27f|1492r|forward|reverse|f|r))$'
)
DUPLICATE_SUFFIX = re.compile(r'\s*\(\d+\)$')  # Finder's "(1)" copies
WINDOWS_RESERVED = re.compile(r'(?i)^(?:con|prn|aux|nul|com[1-9]|lpt[1-9])(?:\.|$)')
SAMPLE_IN_PARENS = re.compile(r'\(([^()]+)\)')
SUMMARY_COLUMNS = [
    'sample_name', 'status', 'forward_file', 'reverse_file',
    'merged_length', 'overlap_columns', 'identity', 'difference_rows',
    'F_raw_length', 'F_trimmed_5prime', 'F_trimmed_3prime', 'F_retained_length',
    'R_raw_length', 'R_trimmed_5prime', 'R_trimmed_3prime', 'R_retained_length',
    'message',
]


def parse_pair_label(filename: str) -> Optional[Tuple[str, str]]:
    """Return (sample_name, F|R), or None for unrecognized filenames.

    Sequencing providers sometimes prepend different numbers to F and R.
    If an explicit '(S1)' token is present, it is preferred as sample name.
    Otherwise, only a leading numerical run index such as '0001_' is removed.
    """
    stem = Path(filename).stem
    while DUPLICATE_SUFFIX.search(stem):
        stem = DUPLICATE_SUFFIX.sub('', stem)
    hit = LABEL_PATTERN.search(stem)
    if not hit:
        return None
    label = next(group for group in hit.groups() if group is not None).lower()
    direction = 'F' if label in ('27f', 'f', 'forward') else 'R'
    prefix = stem[:hit.start()].rstrip('._- ')
    # Some labs include both direction and primer: sample_F_27F / sample_R_1492R.
    # Remove a redundant direction marker, but refuse contradictory markers.
    redundant = LABEL_PATTERN.search(prefix)
    if redundant:
        prior = next(group for group in redundant.groups() if group is not None).lower()
        prior_direction = 'F' if prior in ('27f', 'f', 'forward') else 'R'
        if prior_direction != direction:
            return None
        prefix = prefix[:redundant.start()].rstrip('._- ')
    candidates = [value.strip() for value in SAMPLE_IN_PARENS.findall(prefix)
                  if value.strip() and not value.strip().isdigit()]
    if candidates:
        name = candidates[-1]
    else:
        # Provider export labels may use distinct serial numbers for F and R.
        name = re.sub(r'^\d{3,}[._-]', '', prefix).strip('._- ')
    if not name or name in ('.', '..'):
        return None
    return name, direction


def safe_dirname(name: str) -> str:
    """Keep Unicode sample names but avoid filesystem path separators/control chars."""
    safe = re.sub(r'[^\w.\-]+', '_', name, flags=re.UNICODE).strip('._-')
    if safe in ('', '.', '..'):
        raise ValueError('Sample name cannot be used for a directory: {!r}'.format(name))
    # Windows device names are reserved even when followed by a suffix.
    if WINDOWS_RESERVED.match(safe):
        safe = '_' + safe
    return safe


def discover(input_dir: Path, output_dir: Path, recursive: bool = False):
    files = input_dir.rglob('*') if recursive else input_dir.iterdir()
    group = defaultdict(lambda: {'F': [], 'R': []})
    display_names = {}
    unrecognized = []
    for file in sorted(files, key=lambda p: str(p).casefold()):
        if not file.is_file() or file.suffix.lower() not in EXTENSIONS:
            continue
        # Do not reingest results when output is nested inside the input folder.
        if file == output_dir or output_dir in file.parents:
            continue
        parsed = parse_pair_label(file.name)
        if parsed is None:
            unrecognized.append((str(file), 'Filename lacks a supported F/R suffix'))
            continue
        name, direction = parsed
        key = name.casefold()
        display_names.setdefault(key, name)
        group[key][direction].append(file)
    return [(display_names[key], group[key]) for key in sorted(group)], unrecognized


def summary_row(name: str, direction_files: Dict[str, List[Path]]) -> Dict[str, object]:
    return {
        'sample_name': name, 'status': '',
        'forward_file': '; '.join(map(str, direction_files['F'])),
        'reverse_file': '; '.join(map(str, direction_files['R'])),
        'merged_length': '', 'overlap_columns': '', 'identity': '',
        'difference_rows': '', 'F_raw_length': '', 'F_trimmed_5prime': '',
        'F_trimmed_3prime': '', 'F_retained_length': '', 'R_raw_length': '',
        'R_trimmed_5prime': '', 'R_trimmed_3prime': '', 'R_retained_length': '',
        'message': '',
    }


def write_tsv(path: Path, headers: List[str], rows: List[Dict]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def run_batch(args) -> int:
    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    if not input_dir.is_dir():
        raise ValueError('Input folder not found: {}'.format(input_dir))
    if input_dir == output_dir:
        raise ValueError('--output-dir must differ from --input-dir')
    aggregate_fasta = output_dir / 'all_merged.fasta'
    aggregate_tsv = output_dir / 'all_merged.tsv'
    batch_summary = output_dir / 'batch_summary.tsv'
    issues_tsv = output_dir / 'unrecognized_files.tsv'
    result_paths = (aggregate_fasta, aggregate_tsv, batch_summary, issues_tsv)
    existing = [str(p) for p in result_paths if p.exists()]
    if existing and not args.force:
        raise ValueError('Batch output already exists: {} (use --force to overwrite)'.format(
            ', '.join(existing)))

    groups, unrecognized = discover(input_dir, output_dir, args.recursive)
    if not groups:
        raise ValueError('No files with recognizable paired F/R labels found in {}'.format(
            input_dir))

    # Ensure different source sample names cannot collide in the output tree.
    dirnames = {}
    for name, _ in groups:
        pathpart = safe_dirname(name)
        key = pathpart.casefold()
        if key in dirnames:
            raise ValueError('Sample names have colliding output folder names: {!r} and {!r}'.format(
                name, dirnames[key]))
        dirnames[key] = name

    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    successes = []
    failed = 0
    for idx, (name, pair) in enumerate(groups, 1):
        row = summary_row(name, pair)
        if len(pair['F']) != 1 or len(pair['R']) != 1:
            row['status'] = 'unpaired' if not pair['F'] or not pair['R'] else 'ambiguous'
            row['message'] = 'Found {} forward and {} reverse files; need exactly one of each'.format(
                len(pair['F']), len(pair['R']))
            rows.append(row)
            failed += 1
            print('[{}/{}] {}: {} ({})'.format(idx, len(groups), name, row['status'], row['message']))
            continue

        dirname = safe_dirname(name)
        destination = output_dir / 'samples' / dirname / (dirname + '_merged.fasta')
        local = type(args)(**vars(args))
        local.forward, local.reverse, local.output = str(pair['F'][0]), str(pair['R'][0]), str(destination)
        log_out, log_err = io.StringIO(), io.StringIO()
        try:
            # Keep batch mode legible while preserving full per-sample files.
            with contextlib.redirect_stdout(log_out), contextlib.redirect_stderr(log_err):
                code = run(local)
            summary_path = output_paths(str(destination))['summary']
            report = json.loads(summary_path.read_text(encoding='utf-8')) if summary_path.exists() else {}
            for side, label in [('forward', 'F'), ('reverse', 'R')]:
                stats = report.get(side, {})
                row[label + '_raw_length'] = stats.get('raw_length', '')
                row[label + '_trimmed_5prime'] = stats.get('trimmed_5prime_bases', '')
                row[label + '_trimmed_3prime'] = stats.get('trimmed_3prime_bases', '')
                row[label + '_retained_length'] = stats.get('retained_length', '')
            if code != 0:
                row['status'] = 'failed'
                row['message'] = report.get('error', log_err.getvalue().strip() or 'Assembly failed')
            else:
                assembly = report['assembly']
                sequence = str(SeqIO.read(str(destination), 'fasta').seq)
                row['status'] = 'assembled'
                row['merged_length'] = len(sequence)
                row['overlap_columns'] = assembly['overlap_columns']
                row['identity'] = '{:.6f}'.format(assembly['identity'])
                row['difference_rows'] = assembly['difference_rows']
                successes.append((name, sequence))
        except (ValueError, OSError, IndexError, KeyError) as exc:
            row['status'] = 'failed'
            row['message'] = str(exc)
        if row['status'] != 'assembled':
            failed += 1
        rows.append(row)
        print('[{}/{}] {}: {}{}'.format(idx, len(groups), name, row['status'],
              ' ({} bp)'.format(row['merged_length']) if row['status'] == 'assembled'
              else ' ({})'.format(row['message'])))

    aggregate_fasta.write_text(''.join(to_fasta(safe_dirname(name), sequence) for name, sequence in successes),
                               encoding='utf-8')
    # Exactly two requested columns, each sequence on one TSV row.
    write_tsv(aggregate_tsv, ['样本名称', '拼接后序列'], [
        {'样本名称': name, '拼接后序列': sequence} for name, sequence in successes
    ])
    write_tsv(batch_summary, SUMMARY_COLUMNS, rows)
    write_tsv(issues_tsv, ['file', 'reason'], [
        {'file': file, 'reason': reason} for file, reason in unrecognized
    ])
    print('\nBatch result: {} assembled, {} failed/unpaired, {} unrecognized input files'.format(
        len(successes), failed, len(unrecognized)))
    print('Combined FASTA: {}'.format(aggregate_fasta))
    print('Two-column TSV: {}'.format(aggregate_tsv))
    print('Batch QC report: {}'.format(batch_summary))
    if unrecognized:
        print('Unrecognized files: {}'.format(issues_tsv))
    if not successes:
        return 2
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
