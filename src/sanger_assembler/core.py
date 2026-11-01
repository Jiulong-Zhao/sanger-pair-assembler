"""Quality-aware, reference-free assembly of two Sanger sequencing reads.

The default trimmer reproduces Biopython's abi-trim (Modified Mott) slicing
rule at a *configurable* error-probability cutoff. All coordinates in internal
functions are 0-based; human-readable reports use 1-based coordinates.
"""
from __future__ import annotations

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
