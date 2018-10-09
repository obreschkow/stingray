#!/bin/bash
#
#SBATCH --job-name=medi_stingray

#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --ntasks-per-node=24

#SBATCH --mem=50GB


./stingray make_all -parameterfile parameters-gama.txt
