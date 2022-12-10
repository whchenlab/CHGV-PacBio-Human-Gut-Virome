R1=$1
R2=$2
output=$3
threads=$4
path=/mnt/raid9/sunchuqing/0000_Softwares/2023_Benchmark4_ViromeAssembler/
export PATH=$PATH:$path/metaSPAdes_v3.15.4/SPAdes-3.15.4-Linux/bin/
if [ ! -d $output ]
then
    mkdir -p $output
fi
metaspades.py -1 $R1 -2 $R2 -o $output -t $threads