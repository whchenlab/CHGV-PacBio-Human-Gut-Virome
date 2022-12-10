mkdir -p ${infile}
blastall -p blastn \
    -i $seq\
    -d db/${infile}\
    -o ${infile}/${infile}.m8.tab.1 \
    -m 8 -e 1e-10 -a 40


awk 'BEGIN {FS="\t"} $3>90 && $10>$9 {print $1 "," $2 "\t" $9 "\t" $10 }' ${infile}/${infile}.m8.tab.1 > ${infile}/${infile}.bed.1
awk 'BEGIN {FS="\t"} $3>90 && $9>$10 {print $1 "," $2 "\t" $10 "\t" $9 }' ${infile}/${infile}.m8.tab.1 >> ${infile}/${infile}.bed.1
cat ${infile}/${infile}.bed.1|sort -k 1 -t $'\t' |sed 's/[[:space:]]*$//'|tr -d " "|sort -k1,1 -k2,2n > ${infile}/${infile}.bed.2
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/bedtools merge -i ${infile}/${infile}.bed.2 >${infile}/${infile}.bed
#sed 's/,/\t/g' 03_circle/${infile}/00${count}.gvd.len > ${infile}/${infile}.len
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' $seq|sed 's/\t/ /g' |awk '{print $1","length($NF)}'|sed 's/>//g'|sed 's/,/\t/g'|sort -k 1b,1 > ${infile}/${infile}.len
# cut -f 1 ${infile}/${infile}.bed > 
cut -f 1 ${infile}/${infile}.bed |sort|uniq |sed 's/,/\t/g' |sort -t $'\t' -k 2|sort -k 1b,1  > ${infile}/${infile}.list
join -1 2 -2 1 ${infile}/${infile}.list ${infile}/${infile}.len -o 1.1,1.2,2.2  | sed 's/ /\t/g'|sed 's/\t/,/'> ${infile}/${infile}.len.1
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/bedtools genomecov -i ${infile}/${infile}.bed -g ${infile}/${infile}.len.1 |awk 'BEGIN {FS="\t"} $2>=1 {print $0}'|grep -v "genome" > ${infile}/${infile}.bedtools

awk 'BEGIN {FS="\t"} $3>90 {print $1 "," $2 "\t" $7 "\t" $8 }' ${infile}/${infile}.m8.tab.1 |sort -k 1 -t $'\t' |sed 's/[[:space:]]*$//'|tr -d " "|sort -k1,1 -k2,2n> ${infile}/${infile}.bed.1
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/bedtools merge -i ${infile}/${infile}.bed.1 >${infile}/${infile}.bed
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' $seq|awk '{print $1","length($NF)}'|sed 's/>//g'|sed 's/,/\t/g'|sort -k 1b,1 > ${infile}/${infile}.len
cut -f  1 ${infile}/${infile}.bed |sort|uniq |sed 's/,/\t/g' |sort -t $'\t' -k 2 > ${infile}/${infile}.list
join -1 1 -2 1 ${infile}/${infile}.list ${infile}/${infile}.len -o 1.1,1.2,2.2  | sed 's/ /\t/g'|sed 's/\t/,/'> ${infile}/${infile}.len.1
/mnt/raid6/sunchuqing/Softwares/miniconda3/bin/bedtools genomecov -i ${infile}/${infile}.bed -g ${infile}/${infile}.len.1 |awk 'BEGIN {FS="\t"} $2>=1 {print $0}'|grep -v "genome" > ${infile}/${infile}.bedtools.1    
sort -k 1b,1 ${infile}/${infile}.bedtools  -o ${infile}/${infile}.bedtools 
sort -k 1b,1 ${infile}/${infile}.bedtools.1 -o ${infile}/${infile}.bedtools.1
join ${infile}/${infile}.bedtools  ${infile}/${infile}.bedtools.1 -o 1.1,1.5,2.5 > ${infile}/all.bedtools

cat ${infile}/all.bedtools ../all.withncbi.bedtools ../mgv/mgv/all.bedtools > HGV.all.bedtools
awk 'BEGIN {FS=" "} ($2>0.7 && $3>0.9) || ($2>0.9 && $3>0.7) {print $0}' HGV.all.bedtools  > ./all_cov.csv
sed 's/,/\t/g' ./all_cov.csv | sed 's/ /\t/g' | awk 'BEGIN {FS="\t"} $1!=$2 { if($3>=$4) {print $1 "\t" $2 "\t" $3 } else {print $1 "\t" $2 "\t" $4} }' > mcl.input

mcl  mcl.input.2 --abc -o mcl.out -I 6
/usr/bin/python /mnt/raid8/sunchuqing/2020_Human_phage/0000_scripts/get_VC.py -i mcl.out -o VC.csv
awk 'BEGIN {FS=","} {print $2}' VC.csv > VC_contig.list
grep -vw -Ff VC_contig.list all.contig.list > single.list
sed -i 's/\t/_/g' single.list
cat mcl.out single.list > mcl.sin.output
/usr/bin/python /mnt/raid8/sunchuqing/2020_Human_phage/0000_scripts/get_VC.py -i mcl.sin.output -o VC_sin.csv
grep -wFf crass.list VC_sin.csv
