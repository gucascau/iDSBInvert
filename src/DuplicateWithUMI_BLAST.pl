#!/usr/bin/perl
use strict;
use warnings;
#use Data::Dump qw(dump);
#Author:Xin Wang
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen
### function: This script is use the UMI and SelfBLast to remove duplicates and identify the representative insertion.


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f:s","o:s","i:s","b:s","h:s","t:s","g:s","m:s","c:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{f}  ||!defined $opts{o} ||!defined $opts{i}||!defined $opts{b}|| defined $opts{h}) {
	die "************************************************************************
	Usage: $0.pl -f fasta with ssDNA -i UMI index -b SelfBlast -o Output of Quality Control Reads

	Request Parameters:
	-b Blast Results (Self Blast results))
	-i The UMI information of each inserted reads 
	-f fasta of each inserted reads 
	-o The final results that removed duplicated reads
	
	
	Optional Parameters:
	-t Identity of two reads (default 95)
	-g Gap size (default 2)
	-m mismatches (default 2)
	-c Coverage (Matched size/Full length, default 95%)
	-h Help
************************************************************************\n";
}



########################################################################################
#### The first step to build a potential duplication indexs, here will provide the detailed information whether two reads might be identical.

#### Criterion : Identity more than 95%, less than 2 mismatches and 2 gapsize, two reads coverage more than 95%
########################################################################################

my $blast=$opts{b};
my $umi=$opts{i};
my $fasta=$opts{f};
my $ident=(defined $opts{t})?$opts{t}:95;
my $cover=(defined $opts{c})?$opts{c}:0.95;
my $mismatches = (defined $opts{m})?$opts{m}:2;
my $gap=(defined $opts{g})?$opts{g}:2;
my $output=$opts{o};

#my %remove; my %con; my %contain; my %hash; my %name; my $n=0;

my %hash;
open BLAST,"$blast" or die $!;
while (<BLAST>){
	chomp;
	my ($id1,$id2,$qidentity,$qmatches,$qmismatches,$qgapsize,$qlenth)=(split/\t/,$_)[0,1,2,3,4,5,8];
	next if ($id1 eq $id2);
	my $av=$qmatches/$qlenth;
	### set up stringent paramter, can be adjustable
	next unless ($av>=$cover && $qidentity >=$ident && $qmismatches <=$mismatches && $qgapsize<=$gap);
	
	### put the two reads into the hash
	
	my $str1=join ";",($id1,$id2);
	$hash{$str1}++;
	my $str2=join ";",($id2,$id1);
	$hash{$str2}++;
}

### read the sequene
my %sequence; my $Rid;
open FASTA, "$fasta" or die $!;
while (<FASTA>){
	chomp;
	if(/>(\S+)/){
		$Rid=$1;
	}else{
		$sequence{$Rid}.=$_;
	}
	
}

close FASTA;


########################################################################################
#### The second step to use the UMI to define the unique events

#### Criterion : 1. blast results of two read required to be highly identical, 2. Share the same UMIs 

## 				Attention: Even if some might have one detected UMI
########################################################################################


my %inf; my %qual; my %coverage; my %umiin; my %umipos; my %Fumi;
open UMI ,"$umi" or die $!;
while (<UMI>){
	chomp;
	my ($reads,$counts,$quality,$Um1,$Um2,$Um1Pos,$Um2Pos)=split/\t/,$_;
	
	# put the information into inf index;
	my $strumi =join "\t", ($Um1,$Um2);
	my $strumps =join "\t", ($Um1Pos,$Um2Pos);
	
	$qual{$reads}=$quality;
	$coverage{$reads}=$counts;
	$inf{$reads}=$_;
	$umiin{$reads}=$strumi;
	$umipos{$reads}=$strumps;
	
	# $inf{$reads}->{umi}=$strumi;
# 	$inf{$reads}->{umipos}=$strumi;
	
	# generate the UMI inf;
	
	$Fumi{$strumi}.="$reads;";
	
}


open OUT, ">$output" or die $!;

foreach my $q (keys %Fumi){
	
	my @array=split/;/,$Fumi{$q};
	### we eliminate the blank vaue
	@array= grep { $_ ne '' } @array;
	

	
	### print out those UMI that are only identifed in one read.
	if($#array==1){
		print OUT "$sequence{$array[0]}\t$inf{$array[0]}\n";
		next;
	}
	
	
	my %line;
	
	foreach my $inarray (@array){
		$line{$inarray}++;
	}

	
	### the following is to check the same UMI and whether these reads shared the same insertion, we also select the respresentative reads
	
	my $clsN=@array;
	
	### here we rank the reads first with quality and then read counts
	
	my @keyF = sort { $coverage{$b} <=> $coverage{$a}|| $qual{$b} <=> $qual{$a}  } keys %line;
	
	### push the representative reads to hash, and put the number of reads from each cluster into hash value
	#$uniq{$keyF[0]}=$clsN;
	
	my $reRead=$keyF[0];
	
	my $Fcov=$coverage{$reRead};
	my $Fcovstr=$coverage{$reRead}.";";
	
	my $Fidstr=$reRead.";";
	
	### we integrated the BLAST results and see wether they can truely from the same insertion
	
	my @OtherReads;
	foreach my $f (@keyF){
		
		next if ($f eq $reRead);
		
		my $Fstring=join ";",($f,$reRead);
		
		# the criterion is that whether they can be well blasted.
		if (exists $hash{$Fstring}){		
			$Fcov += $coverage{$f};
			$Fcovstr .= $coverage{$f}.";";
			$Fidstr .=$f.";";
		}else{
			push @OtherReads,$f;
		}	
	}
	
	
	
	### print out the representative reads and read counts
	
	print OUT "$sequence{$reRead}\t$reRead\t$Fcov\t$qual{$reRead}\t$umiin{$reRead}\t$umipos{$reRead}\t$Fidstr\t$Fcovstr\n";
	
	foreach my $q (@OtherReads){
		
		print OUT "$sequence{$q}\t$inf{$q}\n"
	}
	
	
}


close OUT;

