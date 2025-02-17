# Oryza_Pangenome
Sure! Here's a GitHub README based on the provided information:

---

# Oryza Alta Genome Assembly Workflow

This repository contains a set of scripts and configurations for assembling the genome of *Oryza alta* using a combination of tools and workflows for long-read sequencing data analysis. The tools included in this pipeline are:

- **Mecat**, **Canu**, **Arrow**, **Pilon** for CLR (Circular Long Reads) data
- **HiFiAsm** for HiFi data
- Hybrid scaffolding may be included eventually.

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
- **MECAT** (for CLR assembly and correction)
- **Canu** (for CLR assembly)
- **Flye** (for de novo assembly of long reads)
- **Arrow** (for polishing)
- **Pilon** (for polishing)
- **HiFiAsm** (for HiFi long reads)

Make sure to install these tools as required for your specific sequencing data.

---

## Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/oryza-alta-genome-assembly.git
   cd oryza-alta-genome-assembly
   ```

2. Ensure all dependencies (MECAT, Canu, Arrow, Pilon, Flye, etc.) are installed on your system. Instructions for each tool can be found on their respective GitHub pages or documentation.

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

Feel free to open issues or submit pull requests for any improvements or bug fixes!

