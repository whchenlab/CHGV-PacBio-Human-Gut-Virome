fasta=$1
checkv contamination $fasta 05_checkV/$pj -t 16
checkv completeness $fasta 05_checkV/$pj -t 16
checkv complete_genomes $fasta 05_checkV/$pj
checkv quality_summary $fasta 05_checkV/$pj
