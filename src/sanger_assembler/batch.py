"""Batch pairing and assembly with fail-closed sample matching.

Recognizes typical Sanger filename suffixes like sample_F, sample_R,
sample_27F, sample_1492R and sample_(S1)_[27F]. Files that cannot be
unambiguously paired are NEVER guessed or placed in combined outputs.
"""
from __future__ import annotations

import contextlib
import csv
import io
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from .cli import output_paths, run
from .core import to_fasta

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
