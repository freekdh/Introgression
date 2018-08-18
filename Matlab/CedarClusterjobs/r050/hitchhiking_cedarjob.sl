#!/bin/bash -l
#SBATCH --job-name=matlab_r050
#SBATCH --account=def-spotto #adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=2:00:00     #adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  #adjust this if you are using PCT
#SBATCH --mem=8000         #adjust this according to your the memory requirement per node you need
#SBATCH --mail-user=dehaas@zoology.ubc.ca #adjust this to match your email address
#SBATCH --mail-type=ALL

#Load the appropriate matlab module
module load matlab
#Remove -singleCompThread if you are using PCT 
matlab -nodisplay -singleCompThread -r "hitchhiking"
