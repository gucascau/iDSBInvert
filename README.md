# iDSBInvert
The package is to characterize the inverted repeats from ssDNA transformation


## Method
1. Technology design
   ![TechniqueDesign](https://github.com/gucascau/iDSBInvert/assets/23031126/6260030a-10ae-4778-824e-3cfa1ec348ee)

The 69-mer reverse-phase cartridge-purified oligos ordered from Sigma-Aldrich were dissolved to 100 uM in TE buffer. The oligoes CATTGAACAACATGTTGCTGTAAGNNNNACTCATAGTACAGACGGCGTGTATCTGGTATTGGTCTGTAC contain 4-mer random nucleotides in the middle, 8-nt inverted repeats (underlined) at 3’ terminal with 18-nt spacer. To analyze the insertions of transformed DNA, 10mL cells were collected when the density reached ~1.5 × 107 cells/mL in YEP-Raffinose, washed twice with water and mixed with 240 ul 50% PEG3350, 36 ul 1M lithium acetate, 40ul 100 uM DNA and 44 ul water to a total of 360 usl. After incubating at 30°C for 30 min followed by 42°C for 30 min, the cells were spun down, resuspended in 4.2 ml water, and spread 700 ul on 150x15 mm YEP-Galactose plates. A total of 6 plates were spread. The plates were incubated at 30°C for 5 days. All colonies on YP-galactose plates were collected, pooled, and mixed vigorously. Approximately 80 μl cells were spun down for isolating genomic DNA for sequencing as above.

3. Bioinformatic analyses
Raw reads were pre-processed with quality control and merged with PEAR v0.9.111 as described above. Reads with any substring of transformed ssDNA (e.g., CATTGAACAA, CAACATGTTG, TTGCTGTAAG, ACTCATAGTA, GTGTATCTGG, TACCAGATAC, TACTATGAGT, or CTTACAGCAA) were recognized as insertions with ssDNA. UMI sequences were defined between 5bp upstream and 5bp downstream seed sequences (GTAAG (….) ACTCA, and TGAGT (….) CTTAC). Perfect matches were required for these seed sequences. Thereafter, UMI features were used to distinguish between a sequencing duplication of one unique insertion versus two independent insertions. In brief, we only considered the reads with identical UMIs and junction sequences (between HO and insertions) as duplicated reads. If with any difference of either the UMIs or the junction sequences, these inserted reads were as independent events. If the unique events cannot be solved by UMIs, we used BLASTN results for the deduplication as the method described above. For other unclassified reads, i.e., without any UMI or with erroneous UMIs, we performed the BLASTN to build a duplication index, in which two reads were clustered when had >95%, <= 2 mismatches and 2 gap size, and can span 95% of read length. Notably, among the duplications of ssDNA inserts, we selected one pair of reads with highest read quality to represent each unique event.
 
To further characterize two possible mechanisms of the inverted repeats at DSBs, we specifically concentrated on the insertion events that have two clear 4bp UMIs at each side. These insertions with complementary UMIs were retrieved as fold back (inter ssDNA), whereas others were as two different ssDNA (intra ssDNA).


##  Availability 
1. Detect reads with transformed ssDNA
2. Define the unique event of transformed ssDNA
3. Characterize possible mechanisms


## Dependencies

Perl is used to run the scripts. The following softwares are also required:

. bwa
. bbduk
. pear
. blastn

## Usage

````
sh iDSBinvert -a Sample ID -b Work Directory  -f Insertion events -r Output  -p Software installed Directory
````

# Contact


For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu

This project is licensed under the terms of the MIT license.

Copyright (c) 2023 Dr. Kaifu Chen lab

Current version v1.0






