mkdir -p 03_Blast_pVOG/${infile}

ls ../01_Glimmer/| grep -v "\." | grep "." |awk 'BEGIN {FS="."} {print $1}'|sort|uniq > ../${infile}.list
cd $pathh
blastall -p blastx -i ../01_Glimmer/${infile}/${infile}.faa -d /mnt/raid6/sunchuqing/Database/Virus/blastdb/POGseqs \
        -o 03_Blast_pVOG/${infile}/${infile}.m8.tab \
        -m 8 -e 1e-10 -a 16
while read -r line
do
    file="../01_Glimmer/${infile}/${line}"
    filename=`echo "$file" | awk 'BEGIN {FS="/"} {print $NF}'`
    CDSnum=`grep '>' ${file}_predict.fasta |wc -l `
    name=`grep '>' ../01_Glimmer/${infile}/${line} |head -n 1|sed 's/>//g' `
    length=`grep -w "$name" 03_Blast_m8/${infile}/${infile}.len`
    hitnum=`grep "${filename}_" 03_Blast_pVOG/${infile}/${infile}.m8.tab |awk 'BEGIN {FS="\t"} $3>50 {print $1} ' |sort |uniq  |wc -l`
    echo "$length,$CDSnum,$hitnum"|sed 's/\t/,/g'|awk 'BEGIN {FS=","} $3>3 && $2/5000<$3 && $2/5000<$4 {print $0}' >> 03_Blast_pVOG/${infile}/${infile}_cds_num.csv
done <"../01_Glimmer/${infile}.list"

