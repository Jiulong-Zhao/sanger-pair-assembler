import json
import random
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

RNG = random.Random(8675309)
REF = "".join(RNG.choices("ACGT", k=300))


def fixture_pair(tmp_path):
    f = tmp_path / "forward.fasta"
    r = tmp_path / "reverse.fasta"
    f.write_text(">F\n" + REF[:190] + "\n", encoding="utf-8")
    r.write_text(">R\n" + str(Seq(REF[110:]).reverse_complement()) + "\n",
                 encoding="utf-8")
    return f, r


def cli(*args):
    return subprocess.run([sys.executable, "-m", "sanger_assembler", *map(str, args)],
                          capture_output=True, text=True)


def test_cli_creates_all_outputs(tmp_path):
    f, r = fixture_pair(tmp_path)
    output = tmp_path / "merged.fasta"
    run = cli("-f", f, "-r", r, "-o", output)
    assert run.returncode == 0, run.stderr
    assert str(SeqIO.read(output, "fasta").seq) == REF
    for ext in (".F_trimmed.fasta", ".R_revcomp.fasta", ".alignment.txt",
                ".disagreements.tsv", ".summary.json"):
        assert (tmp_path / ("merged" + ext)).exists()
    summary = json.loads((tmp_path / "merged.summary.json").read_text())
    assert summary["status"] == "assembled"
    assert summary["assembly"]["consensus_length"] == 300
    assert summary["assembly"]["overlap_columns"] == 80
    assert not summary["forward"]["quality_available"]


def test_refuse_clobber_without_force(tmp_path):
    f, r = fixture_pair(tmp_path)
    output = tmp_path / "merged.fasta"
    assert cli("-f", f, "-r", r, "-o", output).returncode == 0
    bad = cli("-f", f, "-r", r, "-o", output)
    assert bad.returncode == 2
    assert "--force" in bad.stderr
    assert cli("-f", f, "-r", r, "-o", output, "--force").returncode == 0


def test_no_fake_join_and_failure_summary(tmp_path):
    f, r = fixture_pair(tmp_path)
    output = tmp_path / "refused.fasta"
    failed = cli("-f", f, "-r", r, "-o", output, "--min-overlap", "1000")
    assert failed.returncode == 2
    assert not output.exists()
    summary = json.loads((tmp_path / "refused.summary.json").read_text())
    assert summary["status"] == "failed"
    assert (tmp_path / "refused.F_trimmed.fasta").exists()


def test_help_and_bad_params():
    assert cli("--help").returncode == 0
    assert cli("--version").returncode == 0
    assert cli("-f", "a", "-r", "b", "--trim-cutoff", "2").returncode == 2
