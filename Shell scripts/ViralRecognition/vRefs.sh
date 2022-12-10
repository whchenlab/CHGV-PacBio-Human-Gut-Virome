mkdir -p 03_Blast_m8/${infile}
blastall -p blastn \
    -i $asm \
    -d /mnt/raid6/sunchuqing/Database/Virus/virus \
    -o 03_Blast_m8/${infile}/${infile}.m8.tab \
    -m 8 -e 1e-10 -a 16
grep -v -wFf /mnt/raid6/sunchuqing/Database/Virus/not.list 03_Blast_m8/${infile}/${infile}.m8.tab >tmp && mv tmp 03_Blast_m8/${infile}/${infile}.m8.tab
#echo "chrom start   end" > 03_Blast_m8/${infile}/${infile}.bed
awk 'BEGIN {FS="\t"} $3>50 {print $1 "\t" $7 "\t" $8 }' 03_Blast_m8/${infile}/${infile}.m8.tab |sort -k 1 -t $'\t' |sed 's/[[:space:]]*$//'|tr -d " "|sort -k1,1 -k2,2n> 03_Blast_m8/${infile}/${infile}.bed.1
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/bedtools merge -i 03_Blast_m8/${infile}/${infile}.bed.1 >03_Blast_m8/${infile}/${infile}.bed
sed 's/,/\t/g' 03_circle/${infile}/${infile}.len > 03_Blast_m8/${infile}/${infile}.len
cut -f 1 03_Blast_m8/${infile}/${infile}.len > 03_Blast_m8/${infile}/${infile}.list
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/bedtools genomecov -i 03_Blast_m8/${infile}/${infile}.bed -g 03_Blast_m8/${infile}/${infile}.len |awk 'BEGIN {FS="\t"} $2>=1 {print $0}'|grep -v "genome" > 03_Blast_m8/${infile}/${infile}.bedtools
rm 03_Blast_m8/${infile}/${infile}.genomecov.csv                                                 
awk 'BEGIN {FS=","} $NF>0.9 {print $1}' 03_Blast_m8/${infile}/${infile}.genomecov.csv > 04_Positive/${infile}.ref.genome
