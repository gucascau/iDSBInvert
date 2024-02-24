#### 

id=4N8MH-JKM139-3_S3

softpath=/scratch/ch220812/Software

InvertSrc=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/RawDatasets_20221111Inverted/Update_V210924/InvertedInsertionScripts

srcDir=${softpath}/iDSBins/src
# default genome sequence
genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

# default genome annotation
genomeann=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_version6.bed


Chr3Masked=${softpath}/InsMicro/Database/Chr3.masked.fasta

echo ""
echo "iDSBdetectionInverted: The detection job of ssDNA is started..."

cd /project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/RawDatasets_20221111Inverted/Update_V210924/${id}/Detection

### we scan all the assembled reads that contain CATTGAACAA, CAACATGTTG, TTGCTGTAAG, ACTCATAGTA, GTGTATCTGG, TACCAGATAC, TACTATGAGT,  CTTACAGCAA, then measured the coverage of each sequences that are complete identical, and generate the read information. Here we choose the highest read quality as the representative reads

perl ${InvertSrc}/ExtactssDNAInsertionReads2.pl -i $PWD -m ${id} -o ${id} -q ../QualityControl/${id}_QC.Allquality.stat

echo ""

echo "iDSBdetectionInverted: The detection jobs of ssDNA is finished."
### Then we performed the deduplications of these events using the self blast
# 
 
echo ""
echo "iDSBdetectionInverted: The deduplication of these ssDNA insertions is beginning..."

echo "Define UMI for each ssDNA insertion"

perl ${InvertSrc}/DefineUMI.pl  -i ${id}.ssDNAInsertion.fasta -o ${id}.ssDNAInsertion

echo "Finished the UMI scanning"

echo ""
echo "Check the BLAST results of junction sequences"

# ### define the inserted site
#
# #### two rounds of the blast against the MATa region
#
#
#  ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${id}.ssDNAInsertion.fasta -db ${Chr3Masked} -num_threads 8  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${id}.ssDNAInsertion.MAT.blast
#
#   ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${id}.ssDNAInsertion.fasta -db ${Chr3Masked} -num_threads 8 -task blastn-short  -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${id}.ssDNAInsertion.MAT.short.blast
# #${softpath}/ncbi-blast-2.8.1+/bin/blastn -query   ${SampleID}_final.one.fasta -db ${Chr3Masked} -num_threads 8 -task blastn-short  -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${SampleID}_final.one.twoside.short.blast
#
# #cat ${SampleID}_final.one.twoside.noshortblast.blast ${SampleID}_final.one.twoside.short.blast >${SampleID}_final.one.twoside.finalblast.blast
#
# cat ${id}.ssDNAInsertion.MAT.blast ${id}.ssDNAInsertion.MAT.short.blast >${id}.ssDNAInsertion.MAT.final.blast
#
# ### then extract two edge at the junction regions



#### 1.2 Blast against itself
${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in ${id}.ssDNAInsertion.umi.fasta

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${id}.ssDNAInsertion.umi.fasta -db ${id}.ssDNAInsertion.umi.fasta -task blastn-short -num_threads 8  -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${id}.ssDNAInsertion.umi.blast 



#### Due to he high error rate of insertion we cut the junctions for blast, but the performance is not good enough. 

### 1.1  we combined the two junction sequences of all detected reads
#perl ${srcDir}/Cut_twoedge_sequences_forunassembled_v2.pl  -f ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq -r ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R2.fastq -i ${id}.ssDNAInsertion.umi.txt -o${id}.ssDNAInsertion.umi.cut60.fasta -u ${Cutsize} -t ${CutstartF} -e ${CutstartR}

#${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in ${id}.ssDNAInsertion.umi.cut60.fasta
#${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${id}.ssDNAInsertion.umi.cut60.fasta -db ${id}.ssDNAInsertion.umi.cut60.fasta -task blastn-short -num_threads 8  -word_size 5 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -out ${id}.ssDNAInsertion.umi.cut60.blast 

 # instead of using whole insertion blast, we used the junction sequences
 
#perl DuplciateWithUMI_BLAST.pl  -b ${id}.ssDNAInsertion.umi.cut60.blast  -i ${id}.ssDNAInsertion.umi.txt -f ${id}.ssDNAInsertion.umi.fasta -o ${id}.ssDNAInsertion.umi_removedduplicate.cut60.txt
 

#### we mainly used two information to remove the duplications: one the UMI and one the BLAST results


### If the UMI is same and the junction between the HO end and insertion is same, they are from the same event (same colony)
### If the UMI is same but the junction between the HO end and insertion is different (=>2bp difference. If only 1bp difference, they are likely resulted from sequencing error), they are from independent events
### If the UMI is different, they are independent events

perl ${InvertSrc}/DuplciateWithUMI_BLAST_Cls.pl  -b ${id}.ssDNAInsertion.umi.blast -i ${id}.ssDNAInsertion.umi.txt -f ${id}.ssDNAInsertion.umi.fasta -o ${id}.ssDNAInsertion.umi_removedduplicate_cls.txt  -g 8 -m 8 

#perl DuplciateWithUMI_BLAST.pl  -b ${id}.ssDNAInsertion.umi.blast -i ${id}.ssDNAInsertion.umi.txt -f ${id}.ssDNAInsertion.umi.fasta -o ${id}.ssDNAInsertion.umi_removedduplicate.txt
## WE added the forward and reverse insertion

perl -ne '{chomp; my $seq=(split/\t/,$_)[0]; my $inserttype; if ($seq=~/GGCGTGTATC/||$seq=~/TCTGGTATTG/){$inserttype.="F;"}elsif($seq=~/CAATACCAGA/||$seq=~/GATACACGCC/){$inserttype.="R;"} print "$_\t$inserttype\n"}' ${id}.ssDNAInsertion.umi_removedduplicate_cls.txt >${id}.ssDNAInsertion.umi_removedduplicate_withtype_cls.txt 

### here we removed the insertions that only supported by one read
perl -ne '{chomp; my $ReadCount=(split/\t/,$_)[2]; print "$_\n" if ($ReadCount>1)}' ${id}.ssDNAInsertion.umi_removedduplicate_withtype_cls.txt   >${id}.ssDNAInsertion.umi_removedduplicate_hql_cls.txt 

perl -ne '{chomp; my $ReadCount=(split/\t/,$_)[2]; print "$_\n" if ($ReadCount<=1)}' ${id}.ssDNAInsertion.umi_removedduplicate_withtype_cls.txt   >${id}.ssDNAInsertion.umi_removedduplicate_lql_cls.txt 
### we categorized them into confident clear reverted insertion



#1. clear with two reverted insertion
 perl -ne '{chomp; my ($a,$b)=(split/\t/,$_)[4,5]; if (length $a ==5  && length $b ==5) {print "$_\n"}}' ${id}.ssDNAInsertion.umi_removedduplicate_hql_cls.txt  >${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.txt


 ### we added the complementary sequence 
 
 perl -ne '{chomp; my ($a,$b)=(split/\t/,$_)[4,5]; $a=~s/;//; $b=~s/;//; my $rev=reverse $a; $rev=~tr/ATGCatgc/TACGtacg/; print "$_\t$rev\n"}' ${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.txt >${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.withrev.txt
 
 # 1.1 Clear with two reverted insertion that used model 2:
  perl -ne '{chomp; my ($a,$b)=(split/\t/,$_)[4,5]; $a=~s/;//; $b=~s/;//; my $rev=reverse $a; $rev=~tr/ATGCatgc/TACGtacg/; print "$_\t$rev\n" if ($rev eq $b)}' ${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.txt >${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.Mechanism2.txt
 
  # 1.2 Clear with two reverted insertion that follows model 1:
  
  perl -ne '{chomp; my ($a,$b)=(split/\t/,$_)[4,5]; $a=~s/;//; $b=~s/;//; my $rev=reverse $a; $rev=~tr/ATGCatgc/TACGtacg/; print "$_\t$rev\n" if ($rev ne $b)}' ${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.txt >${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.Mechanism1.txt
 

  mkdir Final_inverted
  
  cd Final_inverted
  cp ../${id}.ssDNAInsertion.umi_removedduplicate_hql_cls.txt ./
  cp ../${id}.ssDNAInsertion.umi_removedduplicate_lql_cls.txt  ./
  
  ##
  ## two or more 
  cp ../${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.withrev.txt ./
  
  cp ../${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.Mechanism2.txt ./
  
  cp ../${id}.ssDNAInsertion.umi_removedduplicate_cls.cleaninvereted.Mechanism1.txt ./
  

