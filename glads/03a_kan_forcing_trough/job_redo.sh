#!/bin/bash
#SBATCH --time=00-01:30:00
#SBATCH --mem=8G
#SBATCH --account=def-gflowers
#SBATCH --ntasks=1

module load matlab/2022a

matlab -r "run_job(0.4583, 0.0296, 3); exit"
