asm=`realpath $1`
threads=$2
#2. VirFinder
mkdir -p 02_virfinder_result
/usr/bin/Rscript  /mnt/raid5/sunchuqing/Buffalo_gut/VirFinder.R $asm 02_virfinder_result/${sample}.csv
awk 'BEGIN {FS=","} $3>0.6 {print $1}' 02_virfinder_result/${sample}.csv|\
    awk 'BEGIN {FS=" "} {print $1}' > 02_virfinder_result/${sample}.virfiner.genome
seqkit grep -f 02_virfinder_result/${sample}.virfiner.genome $asm > 02_virfinder_result/${sample}.virfinder.fa
