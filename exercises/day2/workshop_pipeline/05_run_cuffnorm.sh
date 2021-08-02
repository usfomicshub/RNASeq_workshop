#!/bin/bash
#SBATCH --ntasks=3
#SBATCH --workdir=../Work
#SBATCH --mail-type=ALL
#SBATCH --time=00:30:00
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --mail-user=your_email@usf.edu
#SBATCH --job-name=run_cuffnorm
#SBATCH --output=run_cuffnorm.out

##Close any open program modules and load the modules you needmodule purge
module purge
module load apps/cufflinks/2.2.1

##Variables shared between samples
REF_ANNOTATION=../Reference_Data/PlasmoDB-41_Pfalciparum3D7.gff
NUMBER_OF_PROCESSORS=3
##Sample specific stuff

SAMPLE1_BAM=vehicle_rep1.bam
SAMPLE2_BAM=vehicle_rep2.bam

SAMPLE3_BAM=drug_rep1.bam
SAMPLE4_BAM=drug_rep2.bam

CONTROL_LABEL=Vehicle
EXPERIMENT_LABEL=Drug

OUTPUT=Cuffnorm_Output
NORMALIZTION_METHOD=classic-fpkm

##Cufnorm can be used to create normalized expression measurements
##-p specifies the number of processors
##-o specifies the name of the output directory
##--library-norm-method is the type of normalization you want performed
##-L specifies a comma seperated list of names for the different groups
##The final input is the sequence alignments (bams) cooresponding to the samples
##Samples from the same group are seperated by a comma, groups are 
##seperated by a space
cuffnorm -p $NUMBER_OF_PROCESSORS -o $OUTPUT --library-norm-method $NORMALIZTION_METHOD $REF_ANNOTATION -L $CONTROL_LABEL,$EXPERIMENT_LABEL $SAMPLE1_BAM,$SAMPLE2_BAM $SAMPLE3_BAM,$SAMPLE4_BAM
