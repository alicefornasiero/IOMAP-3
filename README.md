# Oryza_Pangenome

# Genome Assembly Workflow

Workflow of commands for assembling the genomes of the wild *Oryza* species using different software for long-read sequencing data. 

The software used are:
- **Mecat**, **Canu**, **Arrow**, **Pilon** for PacBio CLR (Circular Long Reads) data
- **HiFiAsm** for PacBio HiFi data
- **Bionano optical map** hybrid scaffolding

## Table of Contents
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
   - [MECAT](#mecat)
   - [Flye](#flye)
   - [Canu](#canu)
4. [Configuration Files](#configuration-files)
5. [Notes](#notes)

---

## Requirements
- **MECAT** (for CLR assembly)
- **Canu** (for CLR assembly)
- **Flye** (for CLR assembly)
- **Arrow** (for polishing CLR reads)
- **Pilon** (for polishing CLR reads)
- **HiFiAsm** (for HiFi reads)

Make sure to install these tools as required for your specific sequencing data.

---

## Requirements
The following tools are required for this workflow. You can find installation instructions on their respective GitHub pages:

MECAT (for CLR assembly)
GitHub: https://github.com/xiaochuanle/MECAT

Canu (for CLR assembly)
GitHub: https://github.com/marbl/canu

Arrow (for polishing CLR assembly)
GitHub: https://github.com/pacificbio/arrow

Pilon (for polishing CLR assembly)
GitHub: https://github.com/broadinstitute/pilon

Flye (for de novo assembly of long reads)
GitHub: https://github.com/fenderglass/Flye

HiFiAsm (for HiFi data assembly)
GitHub: https://github.com/PacificBiosciences/HiFiAsm

---

## Usage

### MECAT
For CLR data, use MECAT for assembling and correcting the reads. You can run MECAT with the following commands:

```bash
# Configuring MECAT
mecat.pl config oryza_spp_config_file.txt

# Correcting reads
mecat.pl correct oryza_spp_config_file.txt

# Assembling the genome
mecat.pl assemble oryza_spp_config_file.txt
```

#### MECAT Configuration File (oryza_spp_config_file.txt)
```txt
PROJECT=Oryza_alta
RAWREADS=/ibex/scratch/projects/c2016/andrea_O_alta/Oryza_alta.fasta
GENOME_SIZE=1008000000
THREADS=32
MIN_READ_LENGTH=500
CNS_OVLP_OPTIONS=""
CNS_OPTIONS="-r 0.6 -a 1000 -c 4 -l 2000"
CNS_OUTPUT_COVERAGE=30
TRIM_OVLP_OPTIONS="-B"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
CLEANUP=0
GRID_NODE=0
USE_GRID=false
```

---

### Flye
For Flye assembly, you can use the following command:

```bash
flye --pacbio-raw Oryza_alta.fasta --out-dir FLYE_test --genome-size 1008m --threads 32 --asm-coverage 35
```

---

### Canu
To use Canu for genome assembly, use the following command:

```bash
canu -d canu4 -p canu4 genomeSize=1008m -pacbio-raw ../Oryza_alta.fasta usegrid=1 gridOptions="--time=5-00:00:00 --partition=batch --mem-per-cpu=16g" gridOptionsJobName=Oalta-using-grid
```

---

## Configuration Files
You need to adjust and prepare the configuration files as required by each tool. The example configuration for MECAT is already included above.

---

## Notes
- The pipeline currently supports CLR data assembly using MECAT, Canu, and Flye. 
- HiFi data support via **HiFiAsm** will be added in future updates.
- Hybrid scaffolding may be introduced later to improve the assembly accuracy and contiguity.

---
