#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --account=def-freder19
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192

module load python/3.11
source $HOME/python_env/my_env/bin/activate

python -u analysis_final_gs_2.py

echo "Done" 