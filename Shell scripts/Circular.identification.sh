db=CHGV.filtered
seq=CHGV.filtered.fa
bowtie2-build $seq db/$db
while read -r infile
do
    path=<path to rawdata>
    R1=`ls $path/$infile*R1*`
    R2=`ls $path/$infile*R2*`
    bowtie2 -x db/$db \
        -1 $R1\
        -2 $R2 \
        -S 04_mapping/${infile}_NGS.sam -p 40

done < "sample.list"


#Circular----
while read -r id
do
    bam=$id.bam
    if [ -s "$id.circular" ];then
        continue
    fi
    if [ ! -s $bam ];then
        continue
    fi
    echo "$id"
    if [ ! -s $id.start ];then
      # bedtools.start : contig.id  1 2
        bedtools intersect -a $bam -b bedtools.start|samtools view|awk 'BEGIN {FS="\t"} {print $1 "," $3}'|sort > $id.start
    fi
    if [ ! -s $id.end ];then
      # bedtools.end : contig.id  contigLength-1 contigLength
        bedtools intersect -a $bam -b bedtools.end |samtools view|awk 'BEGIN {FS="\t"} {print $1 "," $3}'|sort> $id.end
    fi
    comm -12 $id.start $id.end > $id.circular
done < "sample.list"
while read -r id
do
    sam=${id}_NGS.sam
    bam=${id}_NGS.bam
    if [ -s $id.NGS.circular ];then
        continue
    fi
    if [ ! -s $sam ];then
        continue
    fi
    echo "$id"
    if [ ! -s $bam ];then
        samtools view -b $sam > $bam
    fi

    bedtools intersect -a $bam -b bedtools.start|samtools view|awk 'BEGIN {FS="\t"} {print $1 "," $3}'|sort > $id.NGS.start
    bedtools intersect -a $bam -b bedtools.end |samtools view|awk 'BEGIN {FS="\t"} {print $1 "," $3}'|sort > $id.NGS.end
    comm -12 $id.NGS.start $id.NGS.end > $id.NGS.circular
done < "sample.list" 
cat *.circular > ../Mapping.circular 
