#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --workdir=../Work
#SBATCH --mail-type=ALL
#SBATCH --time=00:10:00
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --mail-user=JGibbons1@mail.usf.edu
#SBATCH --job-name=FASTQC
#SBATCH --output=fastqc.out

##Close any open program modules and load the modules you need
module purge
module load apps/fastqc/0.11.5

##Run the fastqc command. Use a wildcard to specify want all of your fastqs as input. Specify output location
fastqc ../FASTQs/*.fastq.gz --outdir ../FASTQ_QC


