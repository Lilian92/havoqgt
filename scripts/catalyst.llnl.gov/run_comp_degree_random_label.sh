if [ "$#" -ne 2 ]; then echo "Usage: $0 program-dir $1 sbatch-dir"; exit 1; fi

program="$1"
scripts="$2"
patterns[0]="/p/lustre1/an4/patterns/tree/0"

if [ ! -d "$program" ]; then
    echo "Error: program-dir not found"
    exit 2
fi

if [ ! -d "$scripts" ]; then
    echo "Error: sbatch-dir not found"
    exit 2
fi


#example: (../src is the directory of program)
#./run_comp_degree_random_label ../src

#scales=(28 31 34 37)
#nodes=(2 16 128 1024)

scales=(18 20 22 24)
nodes=(2 2 2 2)

#when flags count is 0, it's running with degree based label
#when flags count is bigger than 0, it's running generating random label in range [0, flags_count)
for i in 0 1 2 3; do
    echo $i
    for pattern in ${patterns[*]}; do
        for flags_count in 0 2 4 8 16 32 64 128 256 512 1024; do
            echo $flags_count
            sbatch -N${nodes[i]}\
                "--export=ARG=-i /dev/shm/rmat${scales[i]} -b /p/lustre1/an4/graph/rmat${scales[i]} -u $flags_count -p $pattern -o ./usr/results,P=$program/run_pattern_matching_beta_1.1" $scripts/pattern_matching.sbatch
        done
    done
done
