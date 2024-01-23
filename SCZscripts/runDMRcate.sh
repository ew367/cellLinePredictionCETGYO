#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=12:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.m.walker@exeter.ac.uk # email me at job completion
#SBATCH --error=DMRcate.err # error file
#SBATCH --output=DMRcate.log # output file
#SBATCH --job-name=DMRcate


# print start date and time
echo Job started on:
date -u
    
# needs to be executed from the scripts folder
echo "Changing Folder to: " $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

# load config file provided on command line when submitting job
#echo "Loading config file for project: " $1
#export PROJECT=$1
#source ./DNAm/config/config.txt 

module load R

# run null CT model on all celltypes
Rscript DMRcateNullCTmodel.r /lustre/projects/Research_Project-MRC190311/DNAm/MRC 'Double-'
#Rscript DMRcateNullCTmodel.r /lustre/projects/Research_Project-MRC190311/DNAm/MRC 'NeuN+'
Rscript DMRcateNullCTmodel.r /lustre/projects/Research_Project-MRC190311/DNAm/MRC 'Sox10+'


# run model LM on NeuN+ and Double-ve
Rscript DMRcateModelLM.r /lustre/projects/Research_Project-MRC190311/DNAm/MRC 'Double-'
Rscript DMRcateModelLM.r /lustre/projects/Research_Project-MRC190311/DNAm/MRC 'NeuN+'

#Rscript DNAm/analysis/EWAS/crr.r ${DATADIR}

# print end date and time
echo Job finished on:
date -u
    