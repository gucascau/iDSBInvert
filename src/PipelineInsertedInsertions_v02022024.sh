#!/bin/bash
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2022 Dr. Kaifu lab                                   
# PI: Kaifu Chen                                                   
# Description: 
#  The package is to characterize the inverted repeats from ssDNA transformation
# Method: 
# Experimental design: The 69-mer reverse-phase cartridge-purified oligos ordered from Sigma-Aldrich were dissolved to 100 uM in TE buffer. The oligoes CATTGAACAACATGTTGCTGTAAGNNNNACTCATAGTACAGACGGCGTGTATCTGGTATTGGTCTGTAC contain 4-mer random nucleotides in the middle, 8-nt inverted repeats (underlined) at 3’ terminal with 18-nt spacer. To analyze the insertions of transformed DNA, 10mL cells were collected when the density reached ~1.5 × 107 cells/mL in YEP-Raffinose, washed twice with water and mixed with 240 ul 50% PEG3350, 36 ul 1M lithium acetate, 40ul 100 uM DNA and 44 ul water to a total of 360 usl. After incubating at 30°C for 30 min followed by 42°C for 30 min, the cells were spun down, resuspended in 4.2 ml water, and spread 700 ul on 150x15 mm YEP-Galactose plates. A total of 6 plates were spread. The plates were incubated at 30°C for 5 days. All colonies on YP-galactose plates were collected, pooled, and mixed vigorously. Approximately 80 μl cells were spun down for isolating genomic DNA for sequencing as above.

# Bioinformatic analysis: Raw reads were pre-processed with quality control and merged with PEAR v0.9.111 as described above. Reads with any substring of transformed ssDNA (e.g., CATTGAACAA, CAACATGTTG, TTGCTGTAAG, ACTCATAGTA, GTGTATCTGG, TACCAGATAC, TACTATGAGT, or CTTACAGCAA) were recognized as insertions with ssDNA. UMI sequences were defined between 5bp upstream and 5bp downstream seed sequences (GTAAG (….) ACTCA, and TGAGT (….) CTTAC). Perfect matches were required for these seed sequences. Thereafter, UMI features were used to distinguish between a sequencing duplication of one unique insertion versus two independent insertions. In brief, we only considered the reads with identical UMIs and junction sequences (between HO and insertions) as duplicated reads. If with any difference of either the UMIs or the junction sequences, these inserted reads were as independent events. If the unique events cannot be solved by UMIs, we used BLASTN results for the deduplication as the method described above. For other unclassified reads, i.e., without any UMI or with erroneous UMIs, we performed the BLASTN to build a duplication index, in which two reads were clustered when had >95%, <= 2 mismatches and 2 gap size, and can span 95% of read length. Notably, among the duplications of ssDNA inserts, we selected one pair of reads with highest read quality to represent each unique event.
#
# To further characterize two possible mechanisms of the inverted repeats at DSBs, we specifically concentrated on the insertion events that have two clear 4bp UMIs at each side. These insertions with complementary UMIs were retrieved as fold back (inter ssDNA), whereas others were as two different ssDNA (intra ssDNA).


################################################################


### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use DSBInverted ***********************************"
	echo "Usage: sh $0 -a Sample ID -b Work Directory  -f Insertion events -r Output  -p Software installed Directory"
	echo "Note: The package is used to characterize the inverted repeats from ssDNA transformation "
	echo ""
	echo -e "Request Parameters:"
	echo -e "\t-a Sample Id (Example: yYY398-B_S10)"
	echo -e "\t-b The working directory, where you put the insertion events"
	echo -e "\t-f The insertion event (Note: We only considered the insertion defined as single donor as the insertions with multiple donors might not be able to have a clean Ty1 locus.)"
	echo -e "\t-p Software installed Path (Default:"", it is the path where you store your blast software)" 

	
	echo ""
	echo ""
	echo -e "Optional Parameters:"
	echo "" 
	echo -e "Optional Parameters -- Ty1 Blast setting up:"
	echo -e "\t-ms identity (Default: 60)"
	echo -e "\t-mc Gapsize  (Default t)"
	echo -e "\t-mb E value  (Default t)"
	echo "" 
	echo -e "Optional Parameters of Ty1 reference"
	echo -e "\t-c Ty1 blast index (Note: if you don't have blast index, please run the 'makeblastdb -in Ty1.fasta -dbtype nucl') (Default: ./LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.fasta)"
	echo -e "\t-d The Ty1 annotation file (Default: ./LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.annotation.bed)"
	echo -e "\t-r Script Stored Path (Default: ./LargeInsertionFeature/Ty1Nuceotide/src)" 
	echo ""

	echo -e "\t-h help"
	
	echo "For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu"
	echo "**************************************************"
   exit 1 # Exit script after printing help
}


# Get absolute path for scripts and check if required scripts exist
in=$PWD

### Set default options
# some general paramaters including mata locus thread, and software path

# default software path 
softwarepath=''

# default parameter parameter for selfblast
### 
DepIden=60
DepGapsize=5
DeEvalue=0.01


while getopts "a:b:c:d:p:f:r:ms:mc:mb" opt
do
   case "$opt" in
      a ) id="$OPTARG" ;;
      b ) WDir="$OPTARG" ;;
	  p ) softpath="$OPTARG" ;;
  	  f ) Insert="$OPTARG" ;;


      ms ) DepIden="$OPTARG" ;;
	  mc ) DepGapsize="$OPTARG" ;;
      mb ) DeEvalue="$OPTARG" ;;
	  c ) TYIndex="$OPTARG" ;;
      d ) TYAnn="$OPTARG" ;;
	  r ) Dscript="$OPTARG" ;;

      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


###single insertion:

## updated the name
#cp wt 3 days.One.txt ${SampleID}.One.txt

# default script path

srcDir=${softpath}/LargeInsertionFeature/Ty1Nuceotide/src

# default Ty1-1 reference sequence
Ty1=${softpath}/LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.fasta

# default Ty1-1 reference annotation
Ty1Ann=${softpath}/LargeInsertionFeature/Ty1Nuceotide/Database/YPLWTY1-1.annotation.bed

# Print helpFunction in case parameters are empty
if [ -z "${id}" ] 
then
   echo "*** error: input Sample ID must be provided ***";
   helpFunction
fi

if [ -z "${WDir}" ] 
then
   echo "*** error: input work path must be provided ***";
   helpFunction
fi


if [ -z "${Insert}" ] 
then
   echo "*** error: input Insertion events with Retrotransposons must be defined ***";
   helpFunction
fi


if [ -z "${softpath}" ] 
then
   echo "*** error: input software path must be provided ***";
   helpFunction
fi


# Begin script in case all parameters are correct
echo "${id}"
echo "${WDir}"
echo "${TYIndex}"
echo "${TYAnn}"
echo "${Insert}"
echo "${Dscript}"
echo "${softpath}"
echo "${Ty1}"
echo "${Ty1Ann}"
echo "${DepIden}"
echo "${DepGapsize}"
echo "${DeEvalue}"
# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"




#### 

id=4N8MH-JKM139-3_S3

softpath=/scratch/ch220812/Software

#InvertSrc=/project/RC_Cardio-Chen-e2/ch220812/Project/Insertion/RawDatasets_20221111Inverted/Update_V210924/InvertedInsertionScripts

srcDir=${softpath}/iDSBInverted/src
# default genome sequence
genomeseq=${softpath}/iDSBInverted/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

# default genome annotation
genomeann=${softpath}/iDSBInverted/Database/saccharomyces_cerevisiae_R64-2-1_20150113_version6.bed


Chr3Masked=${softpath}/iDSBInverted/Database/Chr3.masked.fasta




#### Here is the last version of large insertion events


echo "The job array is started ..."
date


### Set up path file:

echo "Change to the Working Path, where you store your raw fastq reads"

cd ${in}

echo "Create a Sample ID folder that store all the results"
mkdir ${SampleID}
cd ${SampleID}

###############################################################################################################
#### 1. Quality Control
###############################################################################################################
echo ""
echo "iDSBquality: Quality Control is beginning ..."
date

#### 1.1 we removed trimmed low quailty sequences (less than 10)
### also removed all reads that have a 31-mer match to PhiX, allowing two mismatches. 
####detemrined the indext information

### First round the quality control
mkdir QualityControl
cd ${in}/${SampleID}/QualityControl

#fastqc ${in}/../${SampleID}_L001_R1_001.fastq   ${in}/../${SampleID}_L001_R2_001.fastq

### measure the raw reads number :
wc -l ${in}/../${Fread} |awk '{ave=$1/4;print "RawReadnumber\t"ave}' >${in}/../${SampleID}/${SampleID}_FinalReadStat.txt

perl  ${srcDir}/Extract_fastq_withindex.pl -f ${in}/../${Fread} -r ${in}/../${Rread} -i ${Findex} -g ${Rindex} -o ${SampleID}_filter

wc -l ${in}/${SampleID}/QualityControl/${SampleID}_filter.R1.fastq |awk '{ave=$1/4;print "RawRIndexReadnumber\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt


###measure the corrected index reads number:

#### 1.2 filter phix sequences

${softpath}/bbmap/bbduk.sh -in1=${SampleID}_filter.R1.fastq in2=${SampleID}_filter.R2.fastq out1=${SampleID}_filter.unmaped.R1.fastq out2=${SampleID}_filter.unmaped.R2.fastq outm1=${SampleID}_filter.mappedphix.R1.fastq outm2=${SampleID}_filter.mappedphix.R2.fastq ref=${softpath}/bbmap/resources/phix_adapters.fa.gz k=31 hdist=2 stats=stats.txt

wc -l ${in}/${SampleID}/QualityControl/${SampleID}_filter.mappedphix.R1.fastq |awk '{ave=$1/4;print "PhixReads\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt
wc -l ${in}/${SampleID}/QualityControl/${SampleID}_filter.unmaped.R1.fastq |awk '{ave=$1/4;print "RawFilterdReadnumber\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt

#### 1.3 set the stringint quality of beginning of 25 and ending of 15, here we also generate the fasta file of each type of reads 

perl ${srcDir}/Divide_high_lowQuality_v2.pl -f ${SampleID}_filter.unmaped.R1.fastq -r ${SampleID}_filter.unmaped.R2.fastq -u ${Mqmin} -g ${Iqmin} -o ${SampleID}_QC -t ${Matasize}

wc -l ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq |awk '{ave=$1/4;print "FinalhighQualityReads\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt

### meausre the high quality reads number:
echo ""
echo "iDSBquality: Quality Control is finished ..."
date

###############################################################################################################
# 2 .Detection
###############################################################################################################
echo ""
echo "iDSBdetection: Detaction is beginning ..."
date

cd ${in}/${SampleID}/

mkdir Detection

cd ${in}/${SampleID}/Detection

#### 2.1  First round with the unmerged reads (In the aging software, we eliminated the unmerged strategy. We are going to update the method in the second verstion)

## ## We perform the blast for each forward and reverse read
# This step we firstly developed to compare the assembly and unassembly strategies.
# ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ../QualityControl/${SampleID}_QC.highquality.R1.fasta -out ${SampleID}.forward.blast.tbl  -db ${genomeseq} -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1
#
# ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ../QualityControl/${SampleID}_QC.highquality.R2.fasta -out ${SampleID}.reverse.blast.tbl -db ${genomeseq} -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1
#
#
# ### Then detemine the insertion events based on the forward and reverse reads, we could set up different parameter accordingly.
# ### We generated a folder that contains the reads with large insertion events and a statistical file with the information of read quality, identity, mapped size.
# ### This statistical file will be further used to detemine the unique large insertion event.
# ### Here we used all the reads: the reads information put to the fastq file
#
# perl ${srcDir}/Determination_Insertion_nonInsertion-july1_v3.pl -i ${SampleID}.forward.blast.tbl -g ${SampleID}.reverse.blast.tbl  -f ../QualityControl/${SampleID}_QC.highquality.R1.fastq -r ../QualityControl/${SampleID}_QC.highquality.R2.fastq -q ../QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_detected_first -j First

#date

### 2.1 To further confirm the insertion event of each read, second round with the merge methods were applied.
## Assemble All the reads, in this case we would eliminate lots of sequence errors generated at the end of reads

# This will generate some assembled files and unassembled files
${softpath}/pear-0.9.11-linux-x86_64/bin/pear -f ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq -r ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R2.fastq -o ${SampleID}_merged

grep ">" ${SampleID}_merged.assembled.fasta -c |awk '{print "MergedReads\t"$0}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt
wc -l ${SampleID}_merged.unassembled.forward.fastq |awk '{ave=$1/4;print "UnMergedReads\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt


### tranfer into fastq
perl ${srcDir}/fastq2fasta.pl -i ${SampleID}_merged.assembled.fastq -o ${SampleID}_merged.assembled.fasta
perl ${srcDir}/fastq2fasta.pl -i ${SampleID}_merged.unassembled.forward.fastq -o ${SampleID}_merged.unassembled.forward.fasta
perl ${srcDir}/fastq2fasta.pl -i ${SampleID}_merged.unassembled.reverse.fastq -o ${SampleID}_merged.unassembled.reverse.fasta

 

### perform the blast analyes
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_merged.assembled.fasta -out ${SampleID}_merged.assembled.blast.tbl  -db  ${genomeseq} -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_merged.unassembled.forward.fasta -out ${SampleID}_merged.unassembled.forward.blast.tbl  -db ${genomeseq} -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_merged.unassembled.reverse.fasta -out ${SampleID}_merged.unassembled.reverse.blast.tbl -db ${genomeseq} -num_threads ${nproc} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1


cut -f 1 ${SampleID}_merged.assembled.blast.tbl |sort|uniq|wc |awk '{print "MergedMappedReads\t"$1}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt


### Measure the aligned reads:


### Detect the insertion for the second methods

# for assmbled one
perl ${srcDir}/Determination_Insertion_nonInsertion_single_blast_after_merge_v3.pl -i ${SampleID}_merged.assembled.blast.tbl -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -f ${SampleID}_merged.assembled.fasta -o ${SampleID}_detected_second_assembled -j SecondAssembly -c ${Matachr} -t ${Matastart} -e ${Mataend} -m ${Ilength} -n ${Matasize}

# for unassmbled ones
### Here we used unassembled the reads: the reads information put to the fastq file
perl ${srcDir}/Determination_Insertion_nonInsertion-july1_v3.pl -i ${SampleID}_merged.unassembled.forward.blast.tbl -g ${SampleID}_merged.unassembled.reverse.blast.tbl  -f ${SampleID}_merged.unassembled.forward.fastq -r ${SampleID}_merged.unassembled.forward.fastq -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o  ${SampleID}_detected_second_unassembled -j SecondUnAssembly -c ${Matachr} -t ${Matastart} -e ${Mataend} -m ${Ilength} -n ${Matasize}

## 2.3 We further combined two round of large insertion events and confirm whether each read contain insertion event consistently. 

## If reads identified in both method, we confirm these reads as inserted reads. If reads only confirmed only in one method, we will check their read read and determined based on the identity (90) and quality ()

mkdir ${SampleID}_detected_final
cd ${SampleID}_detected_final

cp ../${SampleID}_detected_second_assembled/${SampleID}_detected_second_assembled.assmebled.Ainsertion.quality.txt ../${SampleID}_detected_second_unassembled/${SampleID}_detected_second_unassembled.Ainsertion.quality.txt  ./
cat ${SampleID}_detected_second_assembled.assmebled.Ainsertion.quality.txt ${SampleID}_detected_second_unassembled.Ainsertion.quality.txt >${SampleID}_detected_combined.highqual.txt
## check consistency, we will not check in the first method
#perl ${srcDir}/CheckConsistentInsertion.pl -i ${SampleID}_detected_first.Ainsertion.quality.txt -g ${SampleID}_detected_second.Ainsertion.quality.txt -o ${SampleID}_detected_combined
grep "ID" -v ${SampleID}_detected_combined.highqual.txt|wc |awk '{print "LargeInsertedReads\t"$1}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt

echo ""
echo "iDSBdetection: First detection job array is finished ..."


echo ""
echo "iDSBInverted: The detection job of ssDNA is started..."

### As the unassembled reads are not clear to figure out the ssDNA, we only concentrated on the assembled reads to address the inverted ssDNA events,K
cd ${WDir}

### we scan all the assembled reads that contain CATTGAACAA, CAACATGTTG, TTGCTGTAAG, ACTCATAGTA, GTGTATCTGG, TACCAGATAC, TACTATGAGT,  CTTACAGCAA, then measured the coverage of each sequences that are complete identical, and generate the read information. Here we choose the highest read quality as the representative reads

perl ${srcDir}/ExtactssDNAInsertionReads2.pl -i $PWD -m ${id} -o ${id} -q ../QualityControl/${id}_QC.Allquality.stat

echo ""

echo "iDSBInverted: The detection jobs of ssDNA is finished."
### Then we performed the deduplications of these events using the self blast
# 
 
echo ""
echo "iDSBInverted: The deduplication of these ssDNA insertions is beginning..."

echo "Define UMI for each ssDNA insertion"

perl ${srcDir}/DefineUMI.pl  -i ${id}.ssDNAInsertion.fasta -o ${id}.ssDNAInsertion

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

perl ${srcDir}/DuplciateWithUMI_BLAST_Cls.pl  -b ${id}.ssDNAInsertion.umi.blast -i ${id}.ssDNAInsertion.umi.txt -f ${id}.ssDNAInsertion.umi.fasta -o ${id}.ssDNAInsertion.umi_removedduplicate_cls.txt  -g 8 -m 8 

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
  

