# Test assembly statistics

Computed directly from the installed genomic FASTAs by `prepare_test_sets.py`.
Assemblies are used in full in both test modes. N50 uses all FASTA records;
records may be contigs or scaffolds. GC excludes non-ACGT bases.

[Download the TSV](assembly_stats.tsv)

| Group | Species | Accession | Sequences | Size (Mb) | N50 (bp) | L50 | GC (%) | Ambiguous bases |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hanseniaspora | *Hanseniaspora guilliermondii* | GCA_050947605.1 | 7 | 9.303 | 1,551,527 | 3 | 31.224 | 0 |
| hanseniaspora | *Hanseniaspora menglaensis* | GCA_040256855.1 | 8 | 9.558 | 1,490,982 | 3 | 30.338 | 0 |
| hanseniaspora | *Hanseniaspora opuntiae* | GCA_050947625.1 | 8 | 8.910 | 1,500,692 | 3 | 35.222 | 0 |
| hanseniaspora | *Hanseniaspora osmophila* | GCA_050947655.1 | 8 | 11.407 | 1,658,769 | 3 | 37.279 | 0 |
| hanseniaspora | *Hanseniaspora vineae* | GCA_050947685.1 | 8 | 11.700 | 1,753,121 | 3 | 37.633 | 0 |
| pieris | *Pieris brassicae* | GCF_905147105.1 | 401 | 292.326 | 21,549,451 | 7 | 34.129 | 5,637 |
| pieris | *Pieris napi* | GCF_905475465.1 | 42 | 319.187 | 13,068,865 | 12 | 33.954 | 4,900 |
| pieris | *Pieris rapae* | GCF_905147795.1 | 40 | 256.370 | 10,651,053 | 12 | 33.225 | 1,000 |
| utricularia | *Utricularia aurea* | GCA_056322005.1 | 38 | 180.315 | 8,730,432 | 9 | 37.365 | 2,100 |
| utricularia | *Utricularia gibba* | GCA_982292365.1 | 966 | 118.554 | 5,593,449 | 9 | 40.666 | 0 |
| utricularia | *Utricularia reniformis* | GCA_009725065.1 | 1,825 | 304.539 | 488,004 | 183 | 36.679 | 5,787,106 |
