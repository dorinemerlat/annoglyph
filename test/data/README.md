# Protein evidence for the yeast test profiles

`saccharomyces_uniprot_2000.fasta` contains 2,000 reviewed, nonredundant
Saccharomyces proteins (975,272 amino acids; 1,008,572 bytes). It is intended
for the Saccharomycetales genomes in `test`, `test_full1` and `test_full2`.
The final gene annotations are not expected to match a run with a complete
27-million-protein evidence bank; this is a bounded integration dataset.

Selection preserves complete proteins, accepts lengths from 90 to 5,000 amino
acids, removes exact sequence duplicates, and retains all 16 eligible reviewed
S. paradoxus proteins. The remaining proteins are selected deterministically
by a hash of their accession. FASTA identifiers are UniProt accessions without
spaces, avoiding GALBA header rewriting. The adjacent TSV records the original
UniProt header, taxonomy, length and SHA-256 for every sequence.

Source: the existing downloaded UniProt Saccharomyces FASTA:
`test/test_full2/cache/genes_annotation/structural_annotation/prepare_protein_evidence/saccharomyces_uniprot.fa`.
Its evidence metadata identifies the clade as Saccharomyces, NCBI taxid 4930.
No OrthoDB download is needed for these test profiles.

To reproduce the fixture:

```bash
python test/validation/prepare_test_proteins.py \
  test/test_full2/cache/genes_annotation/structural_annotation/prepare_protein_evidence/saccharomyces_uniprot.fa \
  test/data/saccharomyces_uniprot_2000.fasta
```

`conf/test_proteins.config` supplies the fixture to all three yeast profiles.
Structural and functional gene predictions are cached below
`genes_annotation/test_proteins_2000/`; existing RNA-seq and repeat caches remain
reusable. Automatic evidence generation, if enabled instead of the fixture,
is capped at 2,000 sequences in these test profiles. Production defaults remain
unchanged.

In the real `test_full2` alignment stage, these proteins produced 1,994 miniprot
alignment models for S. cerevisiae and 1,963 for S. paradoxus, with 1,894 and
1,836 candidate training transcripts respectively. These are evidence
alignments, not a count of final validated gene predictions.
