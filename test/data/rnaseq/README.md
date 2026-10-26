# Synthetic RNA-seq fixtures

These files are synthetic, with no network or personal data dependency.
`genome.fa` is a 30 kb sequence drawn from A/C/G/T using Python's
`random.Random(718)`. Good reads are 75 bases long with Phred 40 qualities.
There are 30 single-end reads, plus two paired libraries in `one/` and `two/`
with identical basenames but distinct paths and sequences. The paired reads
come from inward-facing positions on the same synthetic contig.

`mixed/` contains a passing read/pair and deliberate failures for length,
ambiguous bases and low quality. It verifies that a pair is discarded when
either mate fails. Expected retained counts are encoded in
`test_rnaseq_qc.nf`. These tiny inputs test plumbing and file correctness;
they do not validate biological annotation quality.
