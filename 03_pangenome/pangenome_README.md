# *Oryza* Genome Type-Level Pangenome Generation

Workflow of commands to generate genome type-level pangenomes (i.e. AA, BB, CC, DD) of the *Oryza* species.

PGGB and PANACUS software are used to process genome data, perform pangenome analysis and generate visualizations.

## Table of Contents

1. [Requirements](#requirements)
2. [Pipeline Steps](#pipeline-steps)
   1. [PGGB](#1-pggb)
      - [For AA Genome Types](#for-aa-genome-types)
      - [For BB, CC, and DD Genome Types](#for-bb-cc-and-dd-genome-types)
   2. [PANACUS](#2-panacus)
      - [Histgrowth Analysis](#histgrowth-analysis)
      - [Pangenome Visualization](#pangenome-visualization)
3. [References](#references)
4. [Authors](#authors)

## Requirements

The following software are required for this workflow. You can find installation instructions on their respective GitHub pages:

- **PGGB**:
  GitHub: (https://github.com/derika/pggb)
- **PANACUS**:
  GitHub: (https://github.com/michaellyon/panacus)

---

## Pipeline Steps

For each genome type (i.e. AA, BB, CC, DD), generate individual multifasta files for each chromosome.

PGGB is used here to analyze each chromosome individually, and PANACUS to generate pangenome growth statistics and core size estimation.


### 1. PGGB

PGGB was run for each chromosome of each subgenome using the following command:

-i: input file (e.g. CHR##.fasta: input file for each chromosome (replace ## with chromosome number));

-p: minimum identity percentage in the mapping step;

-s: length of the mapped and aligned segment;

-n: number of threads;

--n-mappings: number of mappings per segment of each genome;

-k: removes any match shorter than N bp from the initial alignment;

-o: output directory.


#### For AA genome types:

```bash
pggb -i CHR##.fasta -p 90 -s 15000 -n 13 --n-mappings 1 -k 7 -o out_CHR##_genome_type
```

#### For BB, CC, and DD genome types:

```bash
pggb -i CHR##.fasta -p 80 -s 15000 -n 13 --n-mappings 1 -k 7 -o out_CHR##_genome_type
```

### 2. PANACUS

PANACUS processes the output of PGGB to calculate coverage and pangenome growth and visualize results.

#### Histgrowth Analysis:

Run *panacus histgrowth* to calculate coverage and pangenome growth:

-t: number of threads;

-c: unit (base-pair);

-l: coverage thresholds;

-q: quorum thresholds;

-S: input file (e.g. CHR##.smooth.final.gfa: the output file from PGGB (replace ## with chromosome number)).

``` bash
panacus histgrowth -t 40 -c bp -l 1,1,1 -q 0,0.1,1 -S CHR##.smooth.final.gfa > output.tsv
```

#### Pangenome Visualization:

The resulting output.tsv file is used to visualize the coverage histogram and pangenome growth curve:

```bash
panacus-visualize -e output.tsv > output.pdf
```

---

## References

- **PGGB**: Garrison E et al., *bioRxiv.*, 2023. https://pubmed.ncbi.nlm.nih.gov/37066137/
- **Panacus**: Parmigiani L et al., *Bioinformatics*, 2024. [DOI: 10.1093/bioinformatics/btae720](https://doi.org/10.1093/bioinformatics/btae720)

## Authors

- Andrea Zuccolo - King Abdullah University of Science and Technology
- Alice Fornasiero - King Abdullah University of Science and Technology
