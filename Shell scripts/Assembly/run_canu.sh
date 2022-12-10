pacbio=$1
output=$2
threads=$3

path=<Software path>
sample=`basename $pacbio|cut -d "_" -f 1|cut -d "." -f 1`
if [ ! -d $output ]
then
    mkdir -p $output
fi
export PATH=$PATH:$path/Canu_v2.2/canu/build/bin
#CCSed hifi reads fa
canu  \
    -p ${sample} \
    -d $output genomeSize=20k corOutCoverage=1 \
    -corrected \
    -pacbio $pacbio  \
    useGrid=false

if [ -s $output/${$sample}.unitigs.fasta ];then
  mkdir -p 04_MetaBAT_Assembly/${$sample}
  cd 04_MetaBAT_Assembly/${$sample}
  mkdir -p db
  bowtie2-build ../../$output/${$sample}.unitigs.fasta db/${$sample}
  bowtie2 -x db/${$sample} -U ../../01_ccs/${$sample}.ccs.fastq -S ${$sample}.sam
  samtools view -bS ${$sample}.sam -o ${$sample}.bam
  samtools sort ${$sample}.bam > ${$sample}.sort.bam
  ${Software}/berkeleylab-metabat*/bin/jgi_summarize_bam_contig_depths --outputDepth depth_var.txt ${$sample}.sort.bam
  ${Software}/berkeleylab-metabat*/bin/metabat -i ../../$output/${$sample}.unitigs.fasta -a depth_var.txt  -o metabat -v 

  for file in ./metabat.*.fa
    do
      num=${file//[!0-9]/}
      #echo $num
      sed -e "/^>/ s/$/ ${num}/" metabat.$num.fa  >> metabat_binned.concat.fasta 
  done
  grep '>' metabat_binned.concat.fasta | sed 's/>//g' > metabat_binned.info

  cd ../../05_Assembly/${$sample}
  cut -f1 ../../$output/${$sample}.unitigs.bed | sort|uniq -c > contig.list
  while read -r contig
  do
    num=`echo ${contig} | awk 'BEGIN {FS=" "} {print $1}'`
    tig=`echo ${contig} | awk 'BEGIN {FS="ctg"} {print $2}'`
    if [ $num -eq 1 ];then
      awk '/^>/ {printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' < ../../$output/${$sample}.contigs.fasta |egrep -v '^$'|tr "\t" "\n" >../../$output/${$sample}.n1.contigs.fasta
      grep "$tig" ../../$output/${$sample}.n1.contigs.fasta -A 1|sed 's/>tig/>ctg/g' >> ../${$sample}.fasta
      echo $tig >> infa.contig.list  
      grep "ctg${tig}" ../../$output/${$sample}.unitigs.bed |cut -f4 >>infa.unitig.list
    else
      grep "ctg${tig}" ../../$output/${$sample}.unitigs.bed |cut -f4 >unitig.list
      sed 's/utg/tig/g' unitig.list > uni.list
      binnum=`grep -Ff uni.list ../../04_MetaBAT_Assembly/${$sample}/metabat_binned.info | awk 'BEGIN {FS=" "} {print $2}' |sort |uniq |wc -l`
      inbin=`grep -Ff uni.list ../../04_MetaBAT_Assembly/${$sample}/metabat_binned.info|wc -l`
      uninum=`cat unitig.list | wc -l `
      #echo "$binnum,$uninum,${inbin}"
      if [[ $binnum == 1 && $uninum == $inbin ]];then
      
        grep "$tig" ../../$output/${$sample}.contigs.fasta -A 100| awk -v RS='>' 'NR>1{i++}i==1{print ">"$0}' >> ../${$sample}.fasta
        echo $tig >> infa.contig.list
        cat unitig.list >> infa.unitig.list
      fi
    fi
  
  done <"contig.list" 
  cut -f4 ../../$output/${$sample}.unitigs.bed >unitig.list 
  awk '/^>/ {printf("\n%s\t",$0);next;} {printf("%s",$0);} END {printf("\n");}' < ../../$output/${$sample}.unitigs.fasta |egrep -v '^$'|tr "\t" "\n" >../../$output/${$sample}.n1.unitigs.fasta
  for unitig in `grep -v -Ff infa.unitig.list unitig.list `
  do
    tig=`echo ${unitig} | awk 'BEGIN {FS="utg"} {print $2}'`
    #echo $tig
    grep "$tig" ../../$output/${$sample}.n1.unitigs.fasta -A 1 | sed 's/class=contig/class=unitig/g'|sed 's/>tig/>utg/g' >> ../${$sample}.fasta
  done
