#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --workdir=../Work
#SBATCH --mail-type=ALL
#SBATCH --time=00:30:00
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --mail-user=jgibbons1@usf.edu
#SBATCH --job-name=run_featureCounts
#SBATCH --output=run_featureCounts_exon.out

##Close any open program modules and load the modules you needmodule purge
module purge
module load apps/subread/1.6.0

#Variables shared between samples
REF_ANNOTATION=../Reference_Data/PlasmoDB-41_Pfalciparum3D7.gff
NUMBER_OF_PROCESSORS=1
FEATURE_TYPE=exon
ATTRIBUTE_TYPE=Parent

#Infiles
SAMPLE1_BAM=vehicle_rep1.bam
SAMPLE2_BAM=vehicle_rep2.bam

SAMPLE3_BAM=drug_rep1.bam
SAMPLE4_BAM=drug_rep2.bam

#Outfile
OUTFILE=vehicle_drug_feature_counts_exons.txt

##featureCounts is used to get "raw" counts of reads to genes
##-T specifies the number of processors
##-t specifies the feature type (name in the GFF that cooresponds to
##the features (usualy genes) you want to count reads mapped to
##-a specifies the genome reference file (gff)
## -p means we are using paired-end reads
##Final input is the alignment (bam) files 
featureCounts -T $NUMBER_OF_PROCESSORS -t $FEATURE_TYPE -g $ATTRIBUTE_TYPE -a $REF_ANNOTATION -o $OUTFILE -p $SAMPLE1_BAM $SAMPLE2_BAM $SAMPLE3_BAM $SAMPLE4_BAM
