# Shell scripts
> Shell scripts used in genome assembly and viral recognition.

Raw next generation sequencing of viral reads (referred to as vNGS hereafter) were processed with Trimmomatic v0.3864 (with parameter LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50) to remove adaptors and trim low-quality bases; reads of 50 bp or less after trimming were discarded. The third generation sequencing of viral reads (referred to as vTGS) reads were corrected with CCS using pbccs (v4.0.0, https://github.com/nlhepler/pbccs) with the default parameters.
Putative human reads were identified from the trimmed/CCSed reads by aligning the latter to the human reference genome (hg38; GCA_000001405.15) using Bowtie2 (v2.4.2, --end-to-end) with default parameters and removed from further analysis.

Briefly, IDBA-UD (Release 1.1.3, parameters: --maxk 120 --step 10 –min_contig 1500) was used to assemble the filtered vNGS data. Canu (v2.0-, parameters: genomeSize=20k corOutCoverage=1 -corrected) and Flye (v2.8.2, parameters: --meta --genome-size 20k --min-overlap 1000) were used to assemble the filtered vTGS CCS reads. Because Canu does not have a meta-assembly mode and tends to extend contigs by merging DNA sequences from different viral species to generate erroneous contigs, unitigs were used for subsequent analysis; unitigs are basic blocks of contigs that are shorter but more reliable than contigs ('unitigs' are derived from contigs; wherever a contig end intersects the middle of another contig, the contig is split)69. To further extend the sequences, MetaBAT2 (version 2, default parameters) was used to group unitigs into bins. If all unitigs from one contig could be grouped into the same bin, contigs instead of unitigs were used for further analysis. OPERA-MS (v0.9.0, parameters: -contig-len-thr 1000 --polishing --no-strain-clustering --no-ref-clustering) and metaSpades (v3.13.1, default parameters) were employed for hybrid assemblies using both the vTGS and vNGS datasets from the same samples(Figure S7).
Contigs/unitigs obtained from all the above three strategies were merged; for samples that did not have vTGS data, contigs from the IDBA-UD assembler were used. 
The merged dataset was dereplicated using CD-HIT (v4.8.1, parameters: -c 0.95 -n 8) using a global identity threshold of 95%. 

To identify viral contigs, six independent state-of-the-art viral identification pipelines were used, including VirSorter v2.0 (--min-score 0.7), VirFinder v1.1 (default parameters), and PPR-Meta v1.1ART (default parameters). A BLAST search against the Viral RefSeq genomes was also performed using BLASTn v.2.7.1 with the default parameters and an E-value cutoff of <1e-10; Release 201 (Jul 06, 2020) of the Viral RefSeq database contained 13,148 viral genomes. In addition, the annotated protein sequences were used for BLAST searches against the NCBI POG (Phage Orthologous Groups) database 2013.

A contig was annotated as a virus if it was circular/met at least two out of the following criteria 1-5, adopted from the Gut Virome Database (GVD):
1. VirSorter score ≥ 0.7,
2. VirFinder score > 0.6,
3. PPR-Meta phage score > 0.7,
4. Hits to Viral RefSeq with > 50% identity & > 90% coverage,
5. Minimum of three ORFs, producing BLAST hits to the NCBI POG database 2013 with an E-value of ≤ 1e-5, with at least two per 10 kb of contig length.
6. Alternatively, contigs met one of the above criterium and were annotated as high-quality (≥ 90% completeness) by CheckV were also annotated as viruses.
As short contigs may only represent fragments of viral genomes, contigs that were longer than 5 kb or circular contigs longer than 1.5 kb were selected for further analyses; this dataset was referred to as the Chinese Human Gut Virome (CHGV) dataset, which consisted of a total of 23,513 viral populations 
 
