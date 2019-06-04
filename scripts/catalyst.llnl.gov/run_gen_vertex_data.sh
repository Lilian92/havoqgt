program="./src"
sbatch="./scripts"
graphmatedata="/p/lustre1/an4/graph"

if [ "$#" -ne 0 ]; then
    if [ "$#" -ne 3 ]; then
        echo "Usage: $0 program-dir $1 sbatch-dir $2 graph-matedata-dir";
        exit 1;
    fi
    program="$1"
    sbatch="$2"
    graphmatedata="$3"
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
    echo "Error: graph-matedata-dir not found"
    exit 2
fi

#scales=(28 31 34 37)
#nodes=(2 16 128 1024)

scales=(18 20 22 24 28)

#TODO: need to generate head file
for scale in ${scales[*]}; do
    echo $scale
    n=$(dc -e "2 $scale ^ p");
    echo $n
    sbatch -N1\
        "--export=ARG=-n $n -s 10 -f 10 -o $graphmatedata/rmat${scale},P=$program/generate_vertex_data" $sbatch/gen_vertex_data.sbatch
done
