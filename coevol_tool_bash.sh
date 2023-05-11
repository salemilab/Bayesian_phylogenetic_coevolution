#!/bin/bash
#SBATCH --job-name= coevol_tool   
#SBATCH --mail-type= END,FAIL  
#SBATCH --mail-user=a user@email.com 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb  
#SBATCH --time=7-00:00:00   
#SBATCH --output=coevol_tool_job.%j.out   
#Use modules to load the environment for R
module load R

#Run R script 
Rscript coevol_tool.R