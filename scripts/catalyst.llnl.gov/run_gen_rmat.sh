program="./src"
sbatch="./scripts/"
graph="/p/lustre1/an4/graph"

if [ "$#" -ne 0 ]; then
    if [ "$#" -ne 3 ]; then
        echo "Usage: $0 program-dir $1 sbatch-dir $2 graph-dir";
        exit 1;
    fi
    program="$1"
    sbatch="$2"
    graph="$3"
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

scales=(18 20 22 24)
nodes=(2 2 2 2)
labels=(1 2 3 6 20 100)

for i in 0 1 2 3; do
    for label in ${labels[*]}; do
        echo $i, ${scales[i]}, ${nodes[i]}, ${graph}, ${label}
        sbatch -N${nodes[i]}\
            "--export=ARG=-s ${scales[i]} -e ${label} -p 1 -f 1 -o /dev/shm/rmatlabel${label}scale${scales[i]} -b $graph/rmatlabel${label}scale${scales[i]},P=$program/generate_rmat" $sbatch/gen_rmat.sbatch
    done
done
