pacbio=$1
output=$2
threads=$3

path=/mnt/raid9/sunchuqing/0000_Softwares/2023_Benchmark4_ViromeAssembler/
sample=`basename $pacbio|cut -d "_" -f 1|cut -d "." -f 1`
if [ ! -d $output ]
then
    mkdir -p $output
fi
export PATH=$PATH:$path/metaFlye_v2.9.1/Flye/bin
#CCSed hifi reads fa
flye --meta --pacbio-corr $pacbio -t $threads -r $output