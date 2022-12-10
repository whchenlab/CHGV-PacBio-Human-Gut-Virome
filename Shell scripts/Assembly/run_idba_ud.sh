R1=$1
R2=$2
output=$3
threads=$4
sample=`basename $R1|cut -d "_" -f 1|cut -d "." -f 1`
path=/mnt/raid5/sunchuqing/Softwares/idba/bin/
export PATH=$PATH:$path
if [ ! -d $output ]
then
    mkdir -p $output
fi
fq2fa --merge --filter $R1 $R2 $output/$sample.fa
idba_ud -r $sample.fa -o $output --num_threads $threads \
    --maxk 120 --step 10 --min_contig 1500