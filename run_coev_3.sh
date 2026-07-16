#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --account=def-freder19
#SBATCH --cpus-per-task=19
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1

module load gcc
module load StdEnv


srun="srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task=19"

echo "First Run Started" 

parallel="parallel --delay 0.2 -j 20 --joblog runtask_coev_sim_3.log --resume"


$parallel ~/programs/SLiM_build/slim -d D_p={1} -d D_m={1} -d rep={2} coev_sim_niagara_2.slim :::  0.025 ::: 1 2 3 4 5 6 7 8 9 10

echo "Done"