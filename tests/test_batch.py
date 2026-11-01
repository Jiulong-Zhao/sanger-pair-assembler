"""Synthetic-only integration tests: no real user sequencing files in repo."""
import csv
import json
import random
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from sanger_assembler.batch import parse_pair_label, safe_dirname

RNG = random.Random(98791)
REF1 = ''.join(RNG.choices('ACGT', k=300))
REF2 = ''.join(RNG.choices('ACGT', k=320))


def run_cli(*args):
    return subprocess.run([sys.executable, '-m', 'sanger_assembler', *map(str, args)],
                          text=True, capture_output=True)


def add_pair(folder, f_name, r_name, ref):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f_name).write_text('>F\n' + ref[:190] + '\n', encoding='utf-8')
    (folder / r_name).write_text('>R\n' + str(Seq(ref[110:]).reverse_complement()) + '\n',
                                 encoding='utf-8')


def read_rows(filename):
    with filename.open(encoding='utf-8', newline='') as file:
        return list(csv.DictReader(file, delimiter='\t'))


def test_real_world_filename_patterns():
    assert parse_pair_label('0001_31626102900435_(S1)_[27F](1).ab1') == ('S1', 'F')
    assert parse_pair_label('0002_31626102900435_(S1)_[1492R](1).ab1') == ('S1', 'R')
    assert parse_pair_label('S2_F.ab1') == ('S2', 'F')
    assert parse_pair_label('S2_R.ab1') == ('S2', 'R')
    assert parse_pair_label('isolate-27F.fasta') == ('isolate', 'F')
    assert parse_pair_label('isolate-1492R.fasta') == ('isolate', 'R')
    assert parse_pair_label('isolate_forward.fasta') == ('isolate', 'F')
    assert parse_pair_label('isolate_reverse.fasta') == ('isolate', 'R')
    assert parse_pair_label('sample_F_27F.ab1') == ('sample', 'F')
    assert parse_pair_label('sample_R_1492R.ab1') == ('sample', 'R')
    assert parse_pair_label('sample_R_27F.ab1') is None  # inconsistent direction labels
    assert parse_pair_label('differential.fasta') is None
    assert parse_pair_label('nope_R_final.ab1') is None
    assert safe_dirname('S 1/@@') == 'S_1'


def test_batch_success_fasta_tsv_and_per_sample_reports(tmp_path):
    source = tmp_path / 'input'
    add_pair(source, '0001_316_(S1)_[27F](1).fasta',
             '0002_316_(S1)_[1492R](1).fasta', REF1)
    add_pair(source, 'S2_F.fasta', 'S2_R.fasta', REF2)
    out = tmp_path / 'results'
    cmd = run_cli('--input-dir', source, '--output-dir', out)
    assert cmd.returncode == 0, cmd.stderr + cmd.stdout
    records = list(SeqIO.parse(out / 'all_merged.fasta', 'fasta'))
    assert {r.id: str(r.seq) for r in records} == {'S1': REF1, 'S2': REF2}
    with (out / 'all_merged.tsv').open(encoding='utf-8') as handle:
        assert handle.readline().strip() == '样本名称\t拼接后序列'
    assert read_rows(out / 'all_merged.tsv') == [
        {'样本名称': 'S1', '拼接后序列': REF1},
        {'样本名称': 'S2', '拼接后序列': REF2},
    ]
    rows = read_rows(out / 'batch_summary.tsv')
    assert len(rows) == 2
    assert [row['status'] for row in rows] == ['assembled', 'assembled']
    assert [row['merged_length'] for row in rows] == ['300', '320']
    assert all(row['overlap_columns'] == '80' for row in rows)
    for name in ('S1', 'S2'):
        destination = out / 'samples' / name / (name + '_merged')
        for ext in ('.fasta', '.F_trimmed.fasta', '.R_revcomp.fasta', '.alignment.txt',
                    '.disagreements.tsv', '.summary.json'):
            assert (destination.parent / (destination.name + ext)).exists()
        assert json.loads((destination.parent / (destination.name + '.summary.json')).read_text())['status'] == 'assembled'
    assert read_rows(out / 'unrecognized_files.tsv') == []
    # Do not silently replace prior results.
    rerun = run_cli('--input-dir', source, '--output-dir', out)
    assert rerun.returncode == 2 and '--force' in rerun.stderr
    assert run_cli('--input-dir', source, '--output-dir', out, '--force').returncode == 0


def test_missing_duplicate_and_unrecognized_are_not_merged(tmp_path):
    source = tmp_path / 'raw'
    add_pair(source, 'S1_F.fasta', 'S1_R.fasta', REF1)
    (source / 'S2_F.fasta').write_text('>F\n' + REF1 + '\n', encoding='utf-8')
    add_pair(source, 'S3_F.fasta', 'S3_R.fasta', REF2)
    (source / 'S3_27F.fasta').write_text('>F\n' + REF2 + '\n', encoding='utf-8')
    (source / 'README.txt').write_text('hello\n', encoding='utf-8')
    out = tmp_path / 'result'
    cmd = run_cli('-i', source, '--output-dir', out)
    assert cmd.returncode == 1, cmd.stderr + cmd.stdout
    assert [r.id for r in SeqIO.parse(out / 'all_merged.fasta', 'fasta')] == ['S1']
    assert len(read_rows(out / 'all_merged.tsv')) == 1
    status = {r['sample_name']: r['status'] for r in read_rows(out / 'batch_summary.tsv')}
    assert status == {'S1': 'assembled', 'S2': 'unpaired', 'S3': 'ambiguous'}
    unknown = read_rows(out / 'unrecognized_files.tsv')
    assert len(unknown) == 1 and 'README.txt' in unknown[0]['file']
    assert not (out / 'samples' / 'S2').exists()


def test_failed_assembly_excluded_from_aggregates(tmp_path):
    source = tmp_path / 'source'
    add_pair(source, 'S1_F.fasta', 'S1_R.fasta', REF1)
    add_pair(source, 'S2_F.fasta', 'S2_R.fasta', REF2)
    out = tmp_path / 'output'
    cmd = run_cli('-i', source, '--output-dir', out, '--min-overlap', 100)
    assert cmd.returncode == 2  # zero success
    assert (out / 'all_merged.fasta').read_text() == ''
    assert read_rows(out / 'all_merged.tsv') == []
    assert all(row['status'] == 'failed' for row in read_rows(out / 'batch_summary.tsv'))
    assert (out / 'samples' / 'S1' / 'S1_merged.summary.json').exists()
    assert not (out / 'samples' / 'S1' / 'S1_merged.fasta').exists()
    assert run_cli('-i', source, '--output-dir', out, '--force').returncode == 0
    # Rerunning with --force but a failing threshold must remove previous consensus.
    assert run_cli('-i', source, '--output-dir', out, '--force', '--min-overlap', '1000').returncode == 2
    assert not (out / 'samples' / 'S1' / 'S1_merged.fasta').exists()


def test_recursive_skips_output_folder(tmp_path):
    source = tmp_path / 'raw'
    add_pair(source / 'S1', 'S1_F.fasta', 'S1_R.fasta', REF1)
    out = source / 'assembled_results'
    assert run_cli('-i', source, '--recursive').returncode == 0
    assert run_cli('-i', source, '--recursive', '--force').returncode == 0
    assert [r.id for r in SeqIO.parse(out / 'all_merged.fasta', 'fasta')] == ['S1']
    assert run_cli('-i', source).returncode == 2  # without --recursive can't see nested pair


def test_bad_batch_cli_params(tmp_path):
    assert run_cli('--input-dir', tmp_path, '-f', 'a', '-r', 'b').returncode == 2
    assert run_cli('-i', tmp_path, '--output-dir', tmp_path).returncode == 2
    assert run_cli('--output-dir', tmp_path).returncode == 2


def test_standalone_file_batch_run(tmp_path):
    source = tmp_path / 'input'
    add_pair(source, 'sample_F_27F.fasta', 'sample_R_1492R.fasta', REF1)
    standalone = Path(__file__).resolve().parents[1] / 'sanger_merge.py'
    cmd = subprocess.run([sys.executable, str(standalone), '--input-dir', str(source)],
                         capture_output=True, text=True)
    assert cmd.returncode == 0, cmd.stdout + cmd.stderr
    result = source / 'assembled_results' / 'all_merged.fasta'
    assert str(SeqIO.read(result, 'fasta').seq) == REF1
