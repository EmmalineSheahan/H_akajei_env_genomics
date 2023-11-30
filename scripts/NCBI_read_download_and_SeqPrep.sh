#!/bin/bash

#########################
########HEADER###########
#SBATCH --qos=env6932
#SBATCH --time=72:00:00
#SBATCH --mem=8gb
#########################


#########################
##### Things you need to define in the command line
# $1 = SRA Accession
#########################


#########################
### Load Modules
module load sra
#########################



#########################
##### MAIN SCRIPT #####
#########################
### Create Directories
mkdir SRA_DLs
mkdir BAM_Files

### Download files
fastq-dump -A $1 --split-files -O SRA_DLs

### Merge Reads
/home/james.cahill/gitLibs/SeqPrep/SeqPrep -f SRA_DLs/$1_1.fastq -r SRA_DLs/$1_2.fastq  -1 SRA_DLs/$1.F.fq.gz -2 SRA_DLs/$1.R.fq.gz -s SRA_DLs/$1.M.fq.gz
#########################


#########################
### Cleanup
rm SRA_DLs/$1_1.fastq
rm SRA_DLs/$1_2.fastq
#########################
