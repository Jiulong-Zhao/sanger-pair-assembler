# Changelog

## 0.2.0 — Batch mode and combined exports (2026-11-01)

- New `--input-dir`, `--output-dir`, and `--recursive` options for folder-based batch assembly.
- Recognize lab/provider F/R naming conventions, including `(S1)_[27F]` vs `(S1)_[1492R]` with differing run prefixes.
- Fail-closed pairing for missing/duplicate read files; list unrecognized files.
- Combined `all_merged.fasta` and two-column UTF-8 `all_merged.tsv` (样本名称, 拼接后序列).
- Per-sample reports retained under `samples/<sample>/`, and new `batch_summary.tsv`.
- Standalone `sanger_merge.py` generated from the maintained Python package.
- New synthetic integration tests for batch pairing, partial failure, recursion and safe reruns.

## 0.1.0 — Initial public-ready package (2026-11-01)

- Package entry point `sanger-merge` and `python -m sanger_assembler`.
- AB1/ABI, FASTA, FASTQ, plain `.seq` / `.txt` input.
- Reverse-complemented R read and conservative reference-free F/R assembly.
- Configurable Biopython-compatible Modified Mott AB1 trim cutoff.
- JSON trimming/QC summary and TSV disagreement report with raw AB1 coordinates.
- Guardrails against implausible orientation, internal matches, or weak overlap.
- Chinese/English documentation, synthetic examples, automated tests, MIT license.
