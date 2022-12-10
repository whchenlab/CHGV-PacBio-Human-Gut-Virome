infile=$1
NGS_PATH=$2

mkdir -p 02_trimmed 03_bac_cpn60 03_human_hg38 04_Stats
mkdir -p 05_Removed 06_Assembly 07_CD-HIT

TRIMMO_JAR_FILE='/mnt/raid1/tools/ngs_tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
TRIMMO_ADAPTOR_FILE_PE='/mnt/raid1/tools/ngs_tools/Trimmomatic-0.38/adapters/TruSeq3-PE.fa'
R1=`ls ${NGS_PATH}/*${infile}*R1*|head -n 1 `
R2=`ls ${NGS_PATH}/*${infile}*R2*|head -n 1 `

if [ ! -s 02_trimmed/${infile}_clean.1.fq ];then
    java -jar $TRIMMO_JAR_FILE PE -threads 16 ${R1} ${R2} 02_trimmed/${infile}_clean.1.fq 02_trimmed/${infile}_clean_unpaired.1.fq 02_trimmed/${infile}_clean.2.fq 02_trimmed/${infile}_clean_unpaired.2.fq ILLUMINACLIP:$TRIMMO_ADAPTOR_FILE_PE:2:15:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50
fi


bowtie2 -p 16 --un-conc 05_Removed/${infile}_%.fastq --no-unal -k 20 -x hg38_ref -1 02_trimmed/${infile}_clean.1.fq -2 02_trimmed/${infile}_clean.2.fq >log

G3_PATH=$3

mkdir -p 01_ccs 02_Removed/${infile}
#Run CCS corrction
CCS=`ls ${G3_PATH}/*${infile}*.subreads.bam `

if [ ! -s ${G3_PATH}/CCS/${infile}.subreads.bam ];then
  ${Software}/miniconda3/bin/ccs ${CCS} 01_ccs/${infile}.ccs.fastq -j 16
else 
  CCS=`ls ${G3_PATH}/CCS/${infile}.subreads.bam`
  bedtools bamtofastq -i ${CCS} -fq 01_ccs/${infile}.ccs.fastq
fi
#Remove human genome
bowtie2 -p 16 --un 02_Removed/${infile}/${infile}.fastq  -x /mnt/raid5/sunchuqing/Human_Gut_Phage/ref/hg38_ref -U 01_ccs/${infile}.ccs.fastq  > log
seqtk seq -a 02_Removed/${infile}/${infile}.fastq   |seqkit rmdup -s -o  02_Removed/${infile}/${infile}.fa

