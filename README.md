# Sanger Pair Assembler

[中文说明](README.zh-CN.md) | English

**A small, reference-free command-line tool for assembling a forward (F) and reverse (R) Sanger sequencing read into a consensus sequence.** Designed for routine bidirectional 16S rRNA, ITS, PCR amplicon, and other overlapping Sanger read pairs, including `.ab1` chromatogram files.

> **Scientific status:** Early-stage open-source software (v0.2.0). Assemblies and unresolved indels/mismatches require chromatogram review before interpretation, publication, or database submission. It does not assess mixed templates or automatically remove PCR primer sequences.

## Features

- Accept `.ab1` / `.abi`, `.fasta` / `.fa` / `.fna` / `.fas`, `.fastq` / `.fq`, `.seq`, `.txt`.
- Uses the **original sequencer orientation** of R; reverse-complements it automatically.
- Trims AB1 ends with a **configurable, Biopython-compatible Modified Mott algorithm** (default error probability `0.05`). Reports **5′ and 3′ bases trimmed from each original read**.
- Semi-global alignment with adjustable minimum overlap and identity; **does not concatenate reads when no reliable overlap is found**.
- Uses per-base Phred qualities (if available) to resolve clear disagreements, otherwise writes IUPAC ambiguity codes. Flags all mismatches and indels for review.
- Produces consensus FASTA, both oriented reads, overlap alignment, discrepancy TSV, and JSON QC summary.
- **Batch assembly of a folder** with safe filename-based F/R pairing, per-sample reports, one combined FASTA, and a two-column TSV of merged sequences.
- Works on macOS (Apple Silicon / Intel), Linux, and Windows with compatible Python/Biopython versions. No reference sequence required.

## Install (macOS)

Python **3.8+** is supported at the package level, with `biopython>=1.83,<2` installed according to the Python version. Note that Python 3.8 is end-of-life; use a newer maintained Python if practical.

### Recommended: clean virtual environment

```bash
cd ~/Downloads/sanger-pair-assembler
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

### Existing Conda environment (e.g. named `python3`)

```bash
conda activate python3
cd ~/Downloads/sanger-pair-assembler
python -m pip install -e .
python -c "from Bio import SeqIO; print('Biopython is ready')"
```

**Mac tip:** `python` and `python3` can point to *different interpreters*, even with Conda activated. Install with `python -m pip` and run with that same `python`.

## Quick start

```bash
sanger-merge \
  -f '0001_sample_(S1)_[27F].ab1' \
  -r '0002_sample_(S1)_[1492R].ab1' \
  -o S1_merged.fasta
```

Alternatively: `python -m sanger_assembler -f ... -r ... -o ...`

> **Important:** Supply the raw R read from the sequencing provider. **Do not manually reverse-complement R** before running the command. Quote file paths containing spaces, `(`, `)`, `[` or `]`.

Try the **fully synthetic** example:

```bash
sanger-merge -f examples/toy_F.fasta -r examples/toy_R.fasta -o toy_merged.fasta
```

## Batch assembly (new in v0.2.0)

Use the installed CLI to process a folder of sequencing files:

```bash
sanger-merge --input-dir /path/to/reads --output-dir /path/to/assembled_results
# Optional: include subfolders, or adjust the AB1 trimming threshold
sanger-merge -i /path/to/reads --output-dir /path/to/q20_results --recursive --trim-cutoff 0.01
```

Alternatively, the repository's standalone **`sanger_merge.py`** needs only Biopython (no package installation):

```bash
python sanger_merge.py --input-dir /path/to/reads --output-dir /path/to/results
```

The filename suffix must clearly identify the read direction: `_F` / `_R`, `_27F` / `_1492R`, `_forward` / `_reverse`, or the equivalent bracketed/parenthesized token. Common pairings include:

| Forward | Reverse | Sample name |
|---|---|---|
| `S1_F.ab1` | `S1_R.ab1` | `S1` |
| `S2_27F.ab1` | `S2_1492R.ab1` | `S2` |
| `sample_F_27F.ab1` | `sample_R_1492R.ab1` | `sample` |
| `S3_forward.fasta` | `S3_reverse.fasta` | `S3` |
| `0001_316_(S1)_[27F](1).ab1` | `0002_316_(S1)_[1492R](1).ab1` | `S1` |

For sequencing-provider exports, an explicit `(S1)` sample token takes precedence over different numeric run prefixes. Unrecognized names are recorded instead of guessed. Missing or duplicate direction files are excluded from assembly, along with QC failures. Batch mode searches the immediate input folder by default; use `--recursive` for subfolders. Existing outputs require explicit `--force` to overwrite.

### Batch outputs

```text
assembled_results/
  all_merged.fasta           # All successfully merged sequences, one >sample_name per read
  all_merged.tsv             # Exactly 2 columns: 样本名称, 拼接后序列 (UTF-8)
  batch_summary.tsv          # QC metrics, 5'/3' trimming, success/failure reasons
  unrecognized_files.tsv     # Sequencing files lacking recognized F/R labels
  samples/
    S1/
      S1_merged.fasta
      S1_merged.F_trimmed.fasta
      S1_merged.R_revcomp.fasta
      S1_merged.alignment.txt
      S1_merged.disagreements.tsv
      S1_merged.summary.json
    S2/ ...
```

`all_merged.tsv` contains **only** two fields (`样本名称`, `拼接后序列`), with each complete consensus on a single row. The corresponding FASTA headers are filesystem-safe sample IDs. Only assembled samples are included; inspect `batch_summary.tsv` for failed/unpaired samples. Exit codes: 0 = all recognized pairs assembled, 1 = partial success, 2 = no assembled sequences or invalid input. Unrecognized files are reported and do not change exit status on their own.

---

## Trimming: Q13, Q20, Q30 are not equivalent to hard filters

The default `--trim-cutoff 0.05` corresponds to a base-score error probability boundary of ~Q13.01. `0.01` corresponds to Q20; `0.001` corresponds to Q30. **Modified Mott considers the cumulative quality of a stretch of bases**. It does not guarantee every retained base has Q ≥13, Q20, or Q30.

```bash
# Default / permissive, often helpful for retaining an overlap
sanger-merge -f F.ab1 -r R.ab1 -o q13.fasta --trim-cutoff 0.05

# More stringent: compare assembly stability
sanger-merge -f F.ab1 -r R.ab1 -o q20.fasta --trim-cutoff 0.01

# Strict: may remove the entire overlap
sanger-merge -f F.ab1 -r R.ab1 -o q30.fasta --trim-cutoff 0.001

# Completely disable automatic AB1 end trimming
sanger-merge -f F.ab1 -r R.ab1 -o untrimmed.fasta --no-trim
```

The **default trimming convention reproduces the slice boundaries in Biopython's `abi-trim` implementation** at `0.05`, while making the cutoff adjustable. Trimming is applied only to AB1/ABI files; FASTA has no quality information, and FASTQ qualities are retained but not automatically end-trimmed.

There is **no universally mandatory Sanger end-trimming cutoff**. Phred Q20 is a widely used per-base *quality benchmark*, not a rule that all sequences must be terminally trimmed at Q20. Document your actual trimming and conflict-resolution settings.

## Output files

For `-o S1_merged.fasta`:

| File | Meaning |
|---|---|
| `S1_merged.fasta` | Assembled consensus in F→R direction |
| `S1_merged.F_trimmed.fasta` | F after optional end trimming |
| `S1_merged.R_revcomp.fasta` | R after trimming **and** reverse complement |
| `S1_merged.alignment.txt` | Aligned overlap with match markers |
| `S1_merged.disagreements.tsv` | Conflict and gap positions, Phred scores, and decision rules |
| `S1_merged.summary.json` | Reproducible parameters, original/retained lengths, 5′/3′ cuts and QC |

Coordinates in TSV are **1-based**. `R_raw_pos_1based` refers to the original reverse chromatogram orientation; `R_revcomp_pos_1based` refers to the trimmed reverse-complemented read. A one-sided gap is **not independently confirmed** by the other read and must be reviewed against the peaks.

The program prevents accidental overwrites. To replace previous results, add `--force`. If QC fails, it writes separate oriented reads and a failure summary but **not** a fake consensus.

## All options

```text
--trim-cutoff FLOAT       Mott error-probability cutoff for AB1 (default: 0.05)
--no-trim                 Disable AB1 trimming
--min-overlap INT         Minimum aligned overlap columns (default: 40)
--min-identity FLOAT      Exact A/C/G/T matches / overlap columns (default: 0.90)
--min-q INT               Conflict-resolution min Phred quality (default: 20)
--quality-delta INT       Min Phred advantage for a base decision (default: 5)
--max-r-prefix-clip INT   Max unaligned bases at R* start (default: 15)
--max-both-suffix INT     Max unaligned bases at BOTH suffixes (default: 15)
--input-dir DIR           Assemble all F/R pairs in a folder
--output-dir DIR          Batch output destination
--recursive               Include nested subfolders in batch mode
--force                   Overwrite existing outputs
```

See `sanger-merge --help` for the full argument list.

## Scope and important limitations

1. Each assembled unit is **one F/R pair** with genuine overlap. Batch mode repeats this process for multiple samples; it is not a general contig assembler or variant caller.
2. The read geometry assumes a forward amplicon read and the reverse read from the other end. Non-overlapping reads are rejected; no placeholder `N` bridge is inserted.
3. A quality-trimmed 16S consensus may **still contain primer-derived bases**; identify and remove primers separately as appropriate.
4. A single indel may result from a peak-calling/alignment artifact or genuine sequence difference; all gap sites are written to the TSV and require manual examination.
5. Mixed templates, repeats, homopolymers, low-quality single-read ends and primer artifacts require human review. The percentage identity is not proof of biological accuracy.
6. FASTA lacks Phred qualities. At unresolved substitutions, IUPAC ambiguity codes are used; ambiguous bases in final data may require further review.
7. For clinical, regulatory, or high-stakes use, validate against appropriate controls and an independently verified workflow.

## Development and tests

```bash
python -m pip install -e . pytest
python -m pytest -q
```

A small toy dataset under `examples/` is **synthetic**; no real sample data are included. An optional local test supports `SANGER_PRIVATE_AB1=/path/to/private.ab1 python -m pytest -q` to check default trimming against Biopython's `abi-trim`, without publishing the file.

## Scientific references

- Crossley BM et al. (2020). [Guidelines for Sanger sequencing and molecular assay monitoring](https://doi.org/10.1177/1040638720905833). *J Vet Diagn Invest.*
- Ewing B, Green P (1998). [Base-calling of automated sequencer traces using phred. II. Error probabilities](https://doi.org/10.1101/gr.8.3.186). *Genome Research.*
- [Biopython documentation](https://biopython.org/docs/latest/): `Bio.SeqIO` for ABI and `Bio.Align.PairwiseAligner`.

## License and contribution

MIT License. See [LICENSE](LICENSE) and [CONTRIBUTING.md](CONTRIBUTING.md). Community contributions and issue reports with synthetic data are welcome.
