#!/bin/sh
# SLURM directives
#
#SBATCH --job-name=slurm_job
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4

/home/dobreschkow/stingray/stingray make_all -parameterfile /home/dobreschkow/stingray/parameters.txt