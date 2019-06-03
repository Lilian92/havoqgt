#!/bin/bash
#SBATCH --job-name=gen_rmat
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --distribution=block
#SBATCH --output=./output/output-gen-rmat-%j.out

# example:
# sbatch --export=ARG=<arg>,P=<PROG> degree.sbatch

hostname; date

echo "Running $P $ARG"
date


# Now, run the program appropriately
srun $P $ARG

date
