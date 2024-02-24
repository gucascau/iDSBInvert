#!/usr/bin/perl
#author:wangxin
### function: The script is to extract the MAT sequences and calcuated the reads counts

use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","m:s","q:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{m} || !defined $opts{q}  || !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
			-i: input dir
			-m: index of sample
			-q: insertion quality of reads
			-o: Output string of ssDNA with coverage
************************************************************************\n";
}

my $dir=$opts{i};
my $input1=$opts{m};
#my $input2=$opts{n};
my $output=$opts{o};
my $quality=$opts{q};

my %stringA; my $id1;
print "$dir/$input1\_merged.assembled.fasta\n";
open IN1, "$dir/$input1\_merged.assembled.fasta" or die $!;
while (<IN1>){
	chomp;
	if (/>(\S+)/){
		$id1=$1;
	}else {
		$stringA{$id1}.=substr $_,3, -3;
	}
}

close IN1;

my %hashqual;
open QU, "$quality" or die $!;
while (<QU>){
	chomp;
	my ($Aread,$Rqual)=(split/\t/,$_)[0,5];
	$hashqual{$Aread}=$Rqual;
}

close QU;


### here is for reads with ssDNA insertion
my %hash; my %string;
foreach my $i (keys %stringA){
	next unless ($stringA{$i}=~/CATTGAACAA/|| $stringA{$i}=~/CAACATGTTG/ || $stringA{$i}=~/TTGCTGTAAG/ || $stringA{$i}=~/ACTCATAGTA/|| $stringA{$i}=~/GTGTATCTGG/ || $stringA{$i}=~/TACCAGATAC/|| $stringA{$i}=~/TACTATGAGT/ || $stringA{$i}=~/CTTACAGCAA/);
	my $inf1=$stringA{$i};
	$hash{$inf1}->{A}++;
	$string{$inf1}.="$i;";

}

open OUT, ">$output.ssDNAInsertion.txt" or die $!;
open FA, ">$output.ssDNAInsertion.fasta" or die $!;
print OUT "ssDNARead\tReadQuality\tRidentity\tRMatches\tRLength\tRtype\tRcounts\tssInsertions\tIdenticalReads\tIdentialReadsQuality\n";
foreach my $n (keys %hash){

	my $numA=(exists $hash{$n}->{A})?$hash{$n}->{A}:0;
	#my $numB=(exists $hash{$n}->{B})?$hash{$n}->{B}:0;
	
	my @AllReads=split/;/,$string{$n};
	
	my $Fquality=0; my $readid; my $FqualityString;
	foreach my $Fsort (@AllReads){
		
		next if ($Fsort eq "");
		### sorted the reads
		$readid=($hashqual{$Fsort}>$Fquality)?$Fsort:$readid;
		$FqualityString.= "$hashqual{$Fsort};";
		$Fquality=($hashqual{$Fsort}>$Fquality)?$hashqual{$Fsort}:$Fquality;
		
	}
	#my @AllReadsSorted= sort {$hashqual{$b}<=>$hashqual{$a} } @AllReads;
	my $Flength=length $n;
	print OUT "$readid\t$Fquality\tNA\tNA\t$Flength\tInverted\t$numA\t$n\t$string{$n}\t$FqualityString\n";
	#print OUT "$n\t$output\t$numA\t$string{$n}\t$Fquality\t$FqualityString\n";
	
	print FA ">$readid\t$numA\t$Fquality\n$n\n";
}

close OUT;
close FA;