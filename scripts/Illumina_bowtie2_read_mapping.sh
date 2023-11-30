#!/bin/bash

#########################
########HEADER###########
#SBATCH --qos=env6932
#SBATCH --time=72:00:00
#SBATCH --mem=8gb 
#########################


#########################
### Things that you need to define in the command line ###
# $1 = Accession
# $2 = Reference Genome With Path
# $3 = Output Directory
#########################


#########################
### If you haven't indexed your reference genome you need to do that first.
# bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
#########################


#########################
### Load Modules
module load bowtie2
module load samtools
module load picard
#########################



#########################
##### MAIN SCRIPT #####
### Unzip files
gunzip  SRA_DLs/$1.*.fq.gz # Assumes that you've already downloaded and processed reads with SeqPrep

### Map Reads and filter low quality mappings.
bowtie2 -x $2 -U SRA_DLs/$1.M.fq | samtools view -q 30 -u - | samtools sort -o $3/$1.SE.sort.bam -
samtools rmdup -s $3/$1.SE.sort.bam $3/$1.SE.rmdup.bam

bowtie2 -x $2 -1 SRA_DLs/$1.F.fq -2 SRA_DLs/$1.R.fq | samtools view -q 30 -u - | samtools sort -o $3/$1.PE.sort.bam -
samtools rmdup $3/$1.PE.sort.bam $3/$1.PE.rmdup.bam

### Combine SE and PE reads into a single file
samtools merge $3/$1.merged.unsort.bam  $3/$1.SE.rmdup.bam $3/$1.PE.rmdup.bam
samtools sort -o $3/$1.final.bam $3/$1.merged.unsort.bam 

### Index Bam File
samtools index $3/$1.final.bam

### Read Group Information
picard AddOrReplaceReadGroups \
I=BAM_Files/$1.final.bam \
O=BAM_Files/$2.$1.rg.bam \
RGID=$2 \
RGLB=$1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=$2

samtools index $1.rg.bam
#########################




#########################
### Cleanup
gzip SRA_DLs/$1.*.fq
rm $3/$1.merged.unsort.bam
rm $3/$1.SE.*.bam
rm $3/$1.PE.*.bam
rm BAM_Files/$1.final.bam*
#########################

