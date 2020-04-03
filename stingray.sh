#!/bin/bash
#
#SBATCH --job-name=medi_stingray

#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --ntasks-per-node=24

#SBATCH --mem=100GB


module load hdf5/1.10.2
module load gfortran/6.3.0
module load gcc/6.3.0

./stingray make_all -parameterfile parameters-test.txt
