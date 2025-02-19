# Oryza Genome Type-Level Pangenome Generation

This repository provides a pipeline for generating the pangenome of different genome types (AA, BB, CC, DD) of the *Oryza* species.
It uses the PGGB and PANACUS tools to process the genome data, perform pangenome analysis, and generate visualizations.

## Overview

For each genome type (AA, BB, CC, DD), separate multifasta files for each individual chromosome are generated. 
PGGB is used to analyze each chromosome individually, and PANACUS is then used to further process the pangenome and visualize.

## Requirements

The following software are required for this workflow. You can find installation instructions on their respective GitHub pages:

- **PGGB**:
  GitHub: (https://github.com/derika/pggb)
- **PANACUS**:
  GitHub: (https://github.com/michaellyon/panacus)

---

## Pipeline Steps

### 1. PGGB Analysis

PGGB was run for each chromosome of each genome type to generate genome-type level pangenomes using the following commands:

CHR##.fasta: Input file for each chromosome (replace ## with chromosome number).
-p: Percentage of matching reads to be considered.
-s: Sliding window size.
-n: Number of threads to use.
--n-mappings: Number of mappings to perform.
-k: k-mer size for alignment.
-o: Output prefix for result files.

#### For AA genome types:

```bash
pggb -i CHR##.fasta -p 90 -s 15000 -n 13 --n-mappings 1 -k 7 -o out_CHR##_genome_type
```

#### For BB, CC, and DD genome types:

```bash
pggb -i CHR##.fasta -p 80 -s 15000 -n 13 --n-mappings 1 -k 7 -o out_CHR##_genome_type
```

---

### 2. PANACUS Processing

PANACUS processes the output to generate pangenome statistics and visualization.

CHR##.smooth.final.gfa: Input file from PGGB analysis (replace ## with chromosome number).
-t 40: Number of threads.
-c bp: Units for histogram (base pairs).
-l 1,1,1: Length of sliding windows for analysis.
-q 0,0.1,1: Quantiles for data processing.
-S: Input file format.

#### Histgrowth Analysis:

``` bash
panacus histgrowth -t 40 -c bp -l 1,1,1 -q 0,0.1,1 -S CHR##.smooth.final.gfa > output.tsv
```

---

### 3. Pangenome Visualization

The resulting output.tsv file can be used to generate a visualization of the pangenome using the following command:
This will generate a PDF file (output.pdf) that visualizes the pangenome results.

```bash
panacus-visualize -e output.tsv > output.pdf
```

---

## References

If you use this pipeline in your research, please cite the respective software:

- **PGGB**: Garrison E et al., *bioRxiv.*, 2023. [https://doi.org/10.1101/2023.04.05.535718](https://www.biorxiv.org/content/10.1101/2023.04.05.535718v1)
- **Panacus**: Parmigiani L et al., *Bioinformatics*, 2024. [DOI: 10.1093/bioinformatics/btae720](https://doi.org/10.1093/bioinformatics/btae720)

---

## Author
Andrea Zuccolo - King Abdullah University of Science and Technology
