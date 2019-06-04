if [ "$#" -ne 2 ]; then echo "Usage: $0 program-dir $1 sbatch-dir"; exit 1; fi

program="$1"
scripts="$2"

if [ ! -d "$program" ]; then
    echo "Error: program-dir not found"
    exit 2
fi

if [ ! -d "$scripts" ]; then
    echo "Error: sbatch-dir not found"
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
        "--export=ARG=-n $n -s 10 -f 10 -o /p/lustre1/an4/graph/matedata/rmat${scale},P=$program/generate_vertex_data" $scripts/gen_vertex_data.sbatch
done
