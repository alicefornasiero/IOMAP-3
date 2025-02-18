# Oryza_Pangenome

# Genome Assembly Workflow

Workflow of commands for assembling the genomes of the wild *Oryza* species using different software for PacBio long-read sequencing data. 

The software used are:
- **Mecat**, **Canu**, **Arrow**, **Pilon** for PacBio CLR (Circular Long Reads) data
- **pbmm2**, **Arrow** and **Pilon** for CLR assembly polishing
- **HiFiAsm** for PacBio HiFi data
- **Bionano optical map** hybrid scaffolding

## Table of Contents
1. [Requirements](#requirements)
2. [Usage](#usage)
   - [MECAT](#mecat)
   - [Canu](#canu)
   - [Flye](#flye)
3. [Notes](#notes)

---

## Requirements
The following tools are required for this workflow. You can find installation instructions on their respective GitHub pages:

**1. PacBio CLR assemblying**
**MECAT** (for de novo assembly of CLR reads)
GitHub: https://github.com/xiaochuanle/MECAT

**Canu** (for de novo assembly of CLR reads)
GitHub: https://github.com/marbl/canu

**Flye** (for de novo assembly of CLR reads)
GitHub: https://github.com/fenderglass/Flye

**2. PacBio CLR assembly polishing**
**pbmm2**: For aligning PacBio subreads to the reference contig assembly.
  - GitHub: [https://github.com/PacificBiosciences/pbmm2](https://github.com/PacificBiosciences/pbmm2)
  
**samtools**: For merging and indexing BAM files.
  - GitHub: [https://github.com/samtools/samtools](https://github.com/samtools/samtools)
  
**pbindex**: For indexing PacBio BAM files.
  - GitHub: [https://github.com/PacificBiosciences/pbindex](https://github.com/PacificBiosciences/pbindex)
  
**Arrow**: For polishing genome assemblies using PacBio reads.
  - GitHub: [https://github.com/PacificBiosciences/gcpp](https://github.com/PacificBiosciences/gcpp)

**Pilon** (for polishing CLR read assembly)
GitHub: https://github.com/broadinstitute/pilon

**3. PacBio HiFi assembly**
**HiFiAsm** (for de novo assembly of HiFi reads)
GitHub: https://github.com/PacificBiosciences/HiFiAsm

---

## Usage

**1. PacBio CLR assemblying**

### MECAT
You can run MECAT with the following commands:

```bash
# Configuring MECAT
mecat.pl config oryza_spp_config_file.txt

# Correcting reads
mecat.pl correct oryza_spp_config_file.txt

# Assembling the genome
mecat.pl assemble oryza_spp_config_file.txt
```

#### MECAT Configuration file example
You need to generate and modify the configuration file (e.g. oryza_spp_config_file.txt):

```txt
PROJECT=Oryza_species
RAWREADS=/path/to/directory/Oryza_species.fasta
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
### Canu
Use the following command to run Canu genome assembly:

```bash
time canu -d output -p Oryza_species genomeSize=1008m -pacbio-raw Oryza_species.fasta gnuplot=undef usegrid=1 gridOptions="--time=10-00:00:00 --partition=batch --mem-per-cpu=16g" gridOptionsJobName=Oryza-using-grid
```

---
### Flye
Use the following command to run Flye genome assembly:

```bash
flye --pacbio-raw Oryza_species.fasta --out-dir output --genome-size 1008000000 --threads n_threads --asm-coverage 35
```

---

**2. PacBio CLR assembly polishing**

## Polishing
The pipeline utilizes tools such as **pbmm2**, **samtools**, **pbindex**, and **Arrow** for CLR (Circular Long Reads) assembly polishing.

### Workflow Steps

1. **Align Reads**  
   Align PacBio subreads to the contig assembly using **pbmm2**.
   ```bash
pbmm2 align ${INPUT} ${SAMPLE} ${OUTDIR}/${PREFIX}_${SUFFIX}.sort.bam --sort -j 24 -J 8 ${TMPDIR}"
   ```

2. **Merge Sorted BAM Files**
   Merge the sorted BAM files generated in the previous step.
   ```bash
  samtools merge ${OUTDIR}/${PREFIX}.merged.bam ${OUTDIR}/*.sort.bam"
   ```

3. **Index Merged BAM File Using pbindex**  
   Create an index for the merged BAM file using **pbindex**.
   ```bash
   pbindex ${OUTDIR}/${PREFIX}.merged.bam
   ```

4. **Index Merged BAM File Using samtools**  
   Index the merged BAM file using **samtools**.
   ```bash
  samtools index ${OUTDIR}/${PREFIX}.merged.bam ${OUTDIR}/${PREFIX}.merged.bam.bai
   ```

5. **Index Contig Assembly**
   Index the reference contig assembly using **samtools**.
   ```bash
   samtools faidx ${INPUT}
   ```

6. **Polish Assembly with Arrow**
   Use **Arrow** to polish the genome assembly.
   ```bash
   gcpp --algorithm=arrow -j ${CORES} -r ${INPUT} -o ${OUTDIR}/${PREFIX}.arrow.fasta ${OUTDIR}/${PREFIX}.merged.bam
   ```

---
