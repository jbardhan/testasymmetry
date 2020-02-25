#!/bin/bash
#SBATCH --job-name=cosmo
#SBATCH --output=cB_HB.out
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --time=4-24:00:00
#SBATCH --exclusive
#SBATCH --mem=0
module load matlab/R2019a 
module load anaconda3/3.7
module unload anaconda2/2018.12
matlab -nosplash -nodesktop -r "parpool('local',12);runCosmo_Bondii_HB;exit"
