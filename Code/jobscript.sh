#!/bin/bash
#SBATCH --account=def-jgrons
#SBATCH --time=00:04:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --job-name="sims"

module load nixpgs/16.09
modele load gcc/8.3.0
module load r/4.2.2

Rscript my_sims.R
