#!/bin/bash
#SBATCH --job-name=pattern_matching
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --distribution=block
#SBATCH --output=./output/output-pattern-matching-%j.out

# example:
# sbatch --export=ARG=<arg>,P=<PROG> degree.sbatch

hostname; date

echo "Running: $P $ARG"
date


# Now, run the program appropriately
srun $P $ARG

date
