# Contributing

Thanks for improving Sanger Pair Assembler! Bug reports, documentation updates,
compatibility fixes, and tests are welcome.

1. Open an issue describing the use case, ideally with versions, command, and **synthetic** reproduction data.
2. Create a branch and make a small, focused change.
3. Run `python -m pip install -e . pytest` and `python -m pytest -q`.
4. Open a pull request explaining the behavior before and after.

**Privacy:** do not attach real patient, clinical, unpublished, or other
sensitive AB1/chromatogram data to a public issue without permission. Create a
synthetic reproduction instead.

**Scientific correctness:** do not silently change the default trimming scheme,
consensus behavior, or alignment thresholds. Add tests and explain the consequences
for previously assembled datasets.
