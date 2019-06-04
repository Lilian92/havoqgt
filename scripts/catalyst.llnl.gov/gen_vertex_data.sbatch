#!/bin/bash
#SBATCH --job-name=gen_vertex_data
#SBATCH --time=24:00:00
#SBATCH --distribution=block
#SBATCH --output=./output/output-gen-vertex-data-%j.out

hostname; date

echo "Running $P $ARG"
date


# Now, run the program appropriately
srun $P $ARG

date
