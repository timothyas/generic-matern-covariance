#!/bin/bash

#SBATCH -J slurm-driver
#SBATCH -o slurm-driver.%j.out
#SBATCH -e slurm-driver.%j.err
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p normal
#SBATCH -t 72:00:00


eval "$(conda shell.bash hook)"
conda activate phd

python3 corrpacmap.py
