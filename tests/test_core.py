import random
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from sanger_assembler.core import (
    Read, assemble, consensus_base, mott_trim_bounds,
    read_sequence, reverse_complement,
)

RNG = random.Random(127)
REFERENCE = "".join(RNG.choices("ACGT", k=1000))


def pair_from_forward_coordinates(f, r_forward_strand, q=None):
    # Create the reverse input in the *original* sequencer orientation.
    original_r = str(Seq(r_forward_strand).reverse_complement())
    return (Read(f, q, "F.fasta", len(f)),
            reverse_complement(Read(original_r, q, "R.fasta", len(original_r))))


def test_perfect_dovetail():
    f, rc = pair_from_forward_coordinates(REFERENCE[:650], REFERENCE[350:])
    got = assemble(f, rc)
    assert got["sequence"] == REFERENCE
    assert got["overlap_columns"] == 300
    assert got["identity"] == 1
    assert got["conflicts"] == []


def test_reverse_contained():
    f, rc = pair_from_forward_coordinates(REFERENCE[:650], REFERENCE[350:550])
    assert assemble(f, rc)["sequence"] == REFERENCE[:650]


def test_reverse_extends_past_f():
    f, rc = pair_from_forward_coordinates(REFERENCE[:450], REFERENCE[100:])
    assert assemble(f, rc)["sequence"] == REFERENCE


def test_too_short_overlap_refused():
    f, rc = pair_from_forward_coordinates(REFERENCE[:400], REFERENCE[380:])
    with pytest.raises(ValueError, match="Overlap QC failed"):
        assemble(f, rc)


def test_unrelated_refused():
    other = "".join(RNG.choices("ACGT", k=600))
    f, rc = pair_from_forward_coordinates(REFERENCE[:450], other)
    with pytest.raises(ValueError):
        assemble(f, rc)


def test_wrong_direction_refused():
    f, rc = pair_from_forward_coordinates(REFERENCE[450:], REFERENCE[:550])
    with pytest.raises(ValueError):
        assemble(f, rc)


def test_mismatch_without_quality_is_iupac():
    mutant = list(REFERENCE[350:])
    mutant[120] = "A" if mutant[120] != "A" else "C"
    f, rc = pair_from_forward_coordinates(REFERENCE[:650], "".join(mutant))
    result = assemble(f, rc)
    assert result["mismatches"] == 1
    assert result["conflicts"][0]["rule"] == "ambiguous"
    assert result["sequence"][470] in "RYSWKMBDHVN"


def test_high_quality_preferred():
    assert consensus_base("A", "G", 40, 12) == ("A", "higher_F_quality")
    assert consensus_base("A", "G", 12, 40) == ("G", "higher_R_quality")
    assert consensus_base("A", "G", 22, 21) == ("R", "ambiguous")


def test_indel_reported_with_raw_coordinate():
    mutated = REFERENCE[350:480] + "G" + REFERENCE[480:]
    f, rc = pair_from_forward_coordinates(REFERENCE[:650], mutated)
    result = assemble(f, rc)
    assert result["gaps"] >= 1
    assert any("gap" in item["rule"] for item in result["conflicts"])
    assert all("R_raw_pos_1based" in c for c in result["conflicts"])


def test_trim_nonpositive_region_error():
    with pytest.raises(ValueError, match="No high-quality"):
        mott_trim_bounds([0] * 200, 0.01)


def test_trim_returns_bounds_and_respects_cutoff():
    quality = [2] * 20 + [35] * 200 + [2] * 25
    start1, end1 = mott_trim_bounds(quality, 0.05)
    start2, end2 = mott_trim_bounds(quality, 0.01)
    assert 0 <= start1 < end1 <= len(quality)
    assert 0 <= start2 < end2 <= len(quality)
    assert start1 >= 20 and start2 >= 20


def test_fasta_file_and_orientation(tmp_path: Path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">test\nacgtryn\n", encoding="utf-8")
    r = read_sequence(str(fasta))
    assert r.sequence == "ACGTRYN"
    assert r.qualities is None
    assert reverse_complement(r).sequence == "NRYACGT"
    assert r.report()["trimmed_5prime_bases"] == 0


def test_plain_seq_file(tmp_path: Path):
    seq = tmp_path / "test.seq"
    seq.write_text("ACA TG\nCAT\n", encoding="utf-8")
    assert read_sequence(str(seq)).sequence == "ACATGCAT"


def test_fastq_quality_preserved(tmp_path: Path):
    p = tmp_path / "x.fastq"
    p.write_text("@x\nACGT\n+\nIIII\n", encoding="utf-8")
    r = read_sequence(str(p))
    assert r.qualities == [40, 40, 40, 40]
    assert reverse_complement(r).qualities == [40, 40, 40, 40]
    assert not r.trimmed


def test_abi_trim_default_matches_biopython_if_fixture_available():
    # Optional private integration: never ship real AB1 data to the public repo.
    import os
    file = os.environ.get("SANGER_PRIVATE_AB1")
    if not file:
        pytest.skip("Set SANGER_PRIVATE_AB1 to a locally held private AB1 file")
    trimmed = read_sequence(file, trim=True, cutoff=0.05)
    expected = SeqIO.read(file, "abi-trim")
    assert trimmed.sequence == str(expected.seq).upper()
