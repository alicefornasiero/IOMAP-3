#!/bin/bash
#### Run Mecat2 assembler on PacBio CLR reads ####
#    Alice Fornasiero
#    last modified: February 2025
##################################################

## Software Modules
module load minimap2/2.15
module load samtools/1.8
module load smrtlink/8.0

## Data Variables 
export RAWDATA=/ibex/scratch/projects/c2079/O_malampuzhaensis/3_C01
export INPUT=/ibex/scratch/projects/c2079/analysis/O_malampuzhaensis/mecat2/assembly/4-fsa/contigs.fasta
export OUTDIR=/ibex/scratch/projects/c2079/analysis/O_malampuzhaensis/mecat2/polishing
export TMPDIR=${OUTDIR}/tmp
export PREFIX=O_malampuzhaensis_mecat2
mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR}


## Workflow steps 
for SAMPLE in `ls $RAWDATA/*.subreads.bam`;
do
    SUFFIX=`basename $SAMPLE .subreads.bam`
		
    ## Run pmmb2 to align reads on the contig assembly
    MEM="450gb"
    CORES=32
    JOB1_NAME="pbmm2_align"
    JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME}.${SPECIES} --time=3-00:00:00 --output=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
    JOB1_CMD="pbmm2 align ${INPUT} ${SAMPLE} ${OUTDIR}/${PREFIX}_${SUFFIX}.sort.bam --sort -j 24 -J 8 ${TMPDIR}";
    ${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}"

done

## Merge the sorted .bam files
MEM="450gb"
CORES=1
JOB2_NAME="samtools merge"
JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB2_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB2_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB2_CMD="samtools merge ${OUTDIR}/${PREFIX}.merged.bam ${OUTDIR}/*.sort.bam";
${JOB2_TYPE} --parsable --wrap="${JOB2_CMD}"

## Create index using pbindex
MEM="450gb"
CORES=1
JOB3_NAME="index"
JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME}.${SPECIES} --time=20:00:00 --output=${OUTDIR}/logs/${JOB3_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB3_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB3_CMD="pbindex ${OUTDIR}/${PREFIX}.merged.bam";
${JOB3_TYPE} --parsable --wrap="${JOB3_CMD}"

## Create index using samtools
MEM="450gb"
CORES=1
JOB4_NAME="index"
JOB4_TYPE="sbatch --partition=batch --job-name=${JOB4_NAME}.${SPECIES} --time=20:00:00 --output=${OUTDIR}/logs/${JOB4_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB4_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB4_CMD="samtools index ${OUTDIR}/${PREFIX}.merged.bam ${OUTDIR}/${PREFIX}.merged.bam.bai";
${JOB4_TYPE} --parsable --wrap="${JOB4_CMD}"

## Create index of the contig assembly using samtools
MEM="450gb"
CORES=1
JOB5_NAME="index"
JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME}.${SPECIES} --time=2:00:00 --output=${OUTDIR}/logs/${JOB5_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB5_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB5_CMD="samtools faidx ${INPUT}";
${JOB5_TYPE} --parsable --wrap="${JOB5_CMD}"

## Run Arrow software
MEM="450gb"
CORES=1
JOB6_NAME="run_arrow"
JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB6_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB6_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB6_CMD="gcpp --algorithm=arrow -j ${CORES} -r ${INPUT} -o ${OUTDIR}/${PREFIX}.arrow.fasta ${OUTDIR}/${PREFIX}.merged.bam";
${JOB6_TYPE} --parsable --wrap="${JOB6_CMD}"
