#!/usr/bin/perl
#author:wangxin
### function: The script is to define the UMI in the ssDNA insertion events
use strict;
use warnings;


my $version="1.0 version";
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i:s","o:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i} || !defined $opts{o}) {
       	die "************************************************************************
       	Usage: $0.pl
			-i: input inserted sequences
			-o: Output string of ssDNA with UMI, one for inserted sequence with UMI and another for txt with UMI details
************************************************************************\n";
}

my $input=$opts{i};
#my $input2=$opts{n};
my $output=$opts{o};


my %stringA; my $id1;

open IN, "$input" or die $!;

open FA, ">$output.umi.fasta" or die $!;
open OUT, ">$output.umi.txt" or die $!;
while (<IN>){
	chomp;
	if (/>(\S+)/){
		$id1=$_;
		$id1=~s/>//;
	}else {
		#$stringA{$id1}.=substr $_,3, -3;
		
		#determine the UMI 
		
		my $UMI1=''; my $UMI2=''; my $UMI1num='';my $UMI2num='';
		while($_=~/GTAAG(\w+?)ACTCA/g){
			$UMI1.="$1;";
			$UMI1num.=pos($_).";";
			
		}
		
		while ($_=~/TGAGT(\w+?)CTTAC/g){
			$UMI2.="$1;";
			$UMI2num.=pos($_).";";
		}
		
		$UMI1=($UMI1 ne "")?$UMI1:"NO";
		$UMI2=($UMI2 ne "")?$UMI2:"NO";
		
		$UMI1num=($UMI1num ne "")?$UMI1num:"NO";
		$UMI2num=($UMI2num ne "")?$UMI2num:"NO";
		print OUT "$id1\t$UMI1\t$UMI2\t$UMI1num\t$UMI2num\n";
		print FA ">$id1\t$UMI1\t$UMI2\t$UMI1num\t$UMI2num\n$_\n";
	}
}

close IN;
close OUT;
close FA;