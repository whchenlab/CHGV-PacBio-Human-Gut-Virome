asm=$1
infile=$2
mkdir -p 03_circle/db 03_circle/${infile}
makeblastdb -in $asm -out 03_circle/db/${infile} -dbtype nucl 
blastall -p blastn \
    -i $asm  \
    -d 03_circle/db/${infile} \
    -o 03_circle/${infile}/${infile}.cir.tab \
    -m 8 -e 1e-5 -a 16
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' $asm |awk '{print $1","length($NF)}'|sed 's/>//g' > 03_circle/${infile}/${infile}.len
awk 'BEGIN {FS="\t"} $1==$2 && $7==1 && $3==100.00 {print $0}' 03_circle/${infile}/${infile}.cir.tab > 03_circle/${infile}/${infile}.cir711.tab
rm 03_circle/${infile}/${infile}.cir.len.tab
while read -r line
do
    contig=`echo $line |awk 'BEGIN {FS=","} {print $1}'`
    sed "s/^$contig\t/$line\t/g" 03_circle/${infile}/${infile}.cir711.tab |grep "$line" >> 03_circle/${infile}/${infile}.cir.len.tab
done <"03_circle/${infile}/${infile}.len"
awk 'BEGIN {FS="[,\t]"} ( $2==$10 || $2==$11 ) && $2!=$9 {print $0}' 03_circle/${infile}/${infile}.cir.len.tab > 03_circle/${infile}/${infile}.circle.tab
awk 'BEGIN {FS=","} {print $1}' 03_circle/${infile}/${infile}.circle.tab |sed 's/^/>/g'  > 03_circle/${infile}/${infile}.complete
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0 }' $asm > 03_circle/${infile}/${infile}.fna
