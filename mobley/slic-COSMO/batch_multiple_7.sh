#!/bin/bash
#SBATCH --job-name=m147
#SBATCH --output=output.out
#SBATCH --partition=gpu_test
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --exclusive
#SBATCH --mem=120Gb
module load matlab/R2019a 
module load anaconda3/3.7
module unload anaconda2/2018.12
cd test_multiple_14/test_random_7
matlab -nosplash -nodesktop -r "parpool('local',12);runbenzene;runhexadecane;runchloroform;runcyclohexane;runcarbontet;runhexane;runtoluene;runxylene;runwater;exit"
