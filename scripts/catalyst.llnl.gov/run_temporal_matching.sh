program="./src"
sbatch="./scripts/"
graph="/p/lustre1/an4/graph"
patternp="/p/lustre1/an4/patterns/tree"
if [ "$#" -ne 0 ]; then
    if [ "$#" -ne 4 ]; then
        echo "Usage: $0 program-dir $1 sbatch-dir $2 graph-dir $3 pattern-dir";
        exit 1;
    fi
    program="$1"
    sbatch="$2"
    graph="$3"
    patternp="$4"
fi

if [ ! -d "$program" ]; then
    echo "Error: program-dir not found"
    exit 2
fi

if [ ! -d "$sbatch" ]; then
    echo "Error: sbatch-dir not found"
    exit 2
fi

if [ ! -d "$graph" ]; then
    echo "Error: graph-dir not found"
    exit 2
fi

if [ ! -d "$patternp" ]; then
    echo "Error: pattern-dir not found"
    exit 2
fi

for file in  ../../examples/prunejuice/rmat_log2_tree_pattern/*
do 
    cp -r ${file} ${patternp}
done

patterns[0]=$patternp"/12"  #all local constraints, and global constraints are added to reduce token forwarding
patterns[1]=$patternp"/15"  #edge 0 < 1 < 2 < 3 < 4 < 5
patterns[1]=$patternp"/16"  #same as pattern 15, but with extra checking

#example: (../src is the directory of program)
#./run_comp_degree_random_label ../src

scales=(18 20 22 24)
nodes=(2 2 2 2)
labels=(1 2 3 6 20 100)

for round in 2; do
    for i in 0 1 2 3; do
        echo $i
        for pattern in ${patterns[*]}; do
            echo ${pattern}
            sbatch -N${nodes[i]}\
                "--export=ARG=-i /dev/shm/rmat${scales[i]} -d 1 -b $graph/rmatlabel${label}scale${scales[i]} -p $pattern -o ./usr/results,P=$program/run_pattern_matching_beta_1.1,round=$round" $sbatch/pattern_matching.sbatch
        done
    done
done
