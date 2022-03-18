#!/bin/bash

#SBATCH -J sampledriver
#SBATCH -o sampledriver.%j.out
#SBATCH -e sampledriver.%j.err
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p normal
#SBATCH -t 72:00:00


eval "$(conda shell.bash hook)"
conda activate py38_tim

python3 sample_llc90.py
