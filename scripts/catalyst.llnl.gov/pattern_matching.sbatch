#!/bin/bash
#SBATCH --job-name=pattern_matching
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:20:00
#SBATCH --distribution=block
#SBATCH --output=./output/output_pattern_matching_%j.out

# example:
# sbatch --export=ARG=<arg>,P=<PROG> degree.sbatch

hostname; date

echo "Running: $P $ARG"
echo "Running round: $round"
date


# Now, run the program appropriately
srun $P $ARG

date
