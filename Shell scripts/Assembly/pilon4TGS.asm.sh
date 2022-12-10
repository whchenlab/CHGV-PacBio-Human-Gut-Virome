infile=$1
NGS_PATH=$2
if [ -s ${NGS_PATH}/*${infile}*R1* ];then
  if [ ! -s /02_trimmed/*${infile}*clean.1.fq ];then
    TRIMMO_JAR_FILE='/mnt/raid1/tools/ngs_tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
    TRIMMO_ADAPTOR_FILE_PE='/mnt/raid1/tools/ngs_tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa'
    R1=`ls ${NGS_PATH}/*${infile}*R1*`
    R2=`ls ${NGS_PATH}/*${infile}*R2*`
    mkdir -p ../NGS/02_trimmed
    java -jar $TRIMMO_JAR_FILE PE -threads 4 $R1 $R2 02_trimmed/${infile}_clean.1.fq 02_trimmed/${infile}_clean_unpaired.1.fq 02_trimmed/${infile}_clean.2.fq 02_trimmed/${infile}_clean_unpaired.2.fq ILLUMINACLIP:$TRIMMO_ADAPTOR_FILE_PE:2:15:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50
  fi
  mkdir -p 07_Pilon/${infile}
  cd 07_Pilon/${infile}
  mkdir -p index
  bwa index -p index/draft ../../06_CD-HIT/${infile}/${infile}.fa
  R1=`ls ../../../NGS/02_trimmed/*${infile}*clean.1.fq`
  R2=`ls ../../../NGS/02_trimmed/*${infile}*clean.2.fq`
  bwa mem -t 16 index/draft ${R1} ${R2} | /mnt/raid6/sunchuqing/Softwares/miniconda3/bin/samtools sort -@ 10 -O bam -o align.bam
  samtools index -@ 10 align.bam
  samtools sort -n align.bam > align.sort.bam
  samtools fixmate -m align.sort.bam fixmate.bam
  samtools sort -o align.som.bam fixmate.bam
  samtools markdup align.som.bam align_markdup.bam
  samtools view -@ 10 -q 30 -b align_markdup.bam > align_filter.bam
  samtools index -@ 10 align_filter.bam

  java -Xmx5G -jar ${Software}/pilon-1.23.jar --genome ../../06_CD-HIT/${infile}/${infile}.fa --frags align_filter.bam \
      --fix snps,indels \
      --output ${infile}.pilon 
  cd ../..
fi
