asm=`realpath $1`
threads=$2

cd 05_ViralRecognition
# 1. VirSorter2
mkdir -p 01_virsorter_result
sample=`basename $asm|cut -d "." -f 1|cut -d "_" -f 1`
conda activate vs2
virsorter run \
    -w ./01_virsorter_result/${sample} \
    -i $asm \
    --db-dir /mnt/raid8/sunchuqing/Softwares/Virsorterdb\
    -j $threads --rm-tmpdir\
    --tmpdir ./01_virsorter_result/${sample}/temp --min-score 0.7
conda deactivate
virsorterfa=`ls ./01_virsorter_result/${sample}/Predicted_viral_sequences/*.fa|realpath`
