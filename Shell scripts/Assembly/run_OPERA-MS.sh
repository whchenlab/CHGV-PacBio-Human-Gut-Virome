R1=$1
R2=$2
pacbio=$3
output=$4
threads=$5

path=<Software path>
sample=`basename $pacbio|cut -d "_" -f 1|cut -d "." -f 1`

export PATH=$PATH:$path/OPERA-MS/
#CCSed hifi reads fa
rawpath=`pwd`
if [ ! -d $output ]
then
    mkdir -p $output
fi
cd $output
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/perl $path/OPERA-MS/OPERA-MS.pl \
    --short-read1 $R1 \
    --short-read2 $R2 \
    --long-read $pacbio \
    --out-dir $output \
    --contig-len-thr 1500 \
    --num-processors $threads \
    --polishing --no-strain-clustering --no-ref-clustering

cd $rawpath
