#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --array=1-1000
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --mem-per-cpu=6144M
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END

sed "${SLURM_ARRAY_TASK_ID}q;d" Control_vs_Bipolar_permutations_density_threshold.cmdlist | bash