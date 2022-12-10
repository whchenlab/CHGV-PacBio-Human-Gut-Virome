pacbio=$1
output=$2
threads=$3

path=/mnt/raid9/sunchuqing/0000_Softwares/2023_Benchmark4_ViromeAssembler/
sample=`basename $pacbio|cut -d "_" -f 1|cut -d "." -f 1`
if [ ! -d $output ]
then
    mkdir -p $output
fi
export PATH=$PATH:$path/Canu_v2.2/canu/build/bin
#CCSed hifi reads fa
canu  \
    -p ${sample} \
    -d $output genomeSize=20k corOutCoverage=1 \
    -corrected \
    -pacbio $pacbio  \
    useGrid=false