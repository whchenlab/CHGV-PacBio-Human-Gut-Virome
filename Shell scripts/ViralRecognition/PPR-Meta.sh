asm=`realpath $1`
threads=$2

mkdir -p 03_pprmeta_result
pathh=`pwd`
conda activate /mnt/raid6/sunchuqing/Softwares/miniconda3/envs/tensorflow
cd /mnt/raid6/sunchuqing/Softwares/PPR-Meta
./PPR_Meta $asm $pathh/03_pprmeta_result/${infile}.csv
cd $pathh
awk 'BEGIN {FS=","} $3>0.7 {print $1}'  03_pprmeta_result/${infile}.csv |awk 'BEGIN {FS=" "} {print $1}' > 03_pprmeta_result/${infile}.pprmeta.genome
seqkit grep -f 03_pprmeta_result/${infile}.pprmeta.genome $asm > 03_pprmeta_result/${infile}.pprmeta.fa
pprmetafa=`realpath 03_pprmeta_result/${infile}.pprmeta.fa`
conda deactivate 
