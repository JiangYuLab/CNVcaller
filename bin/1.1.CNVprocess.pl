#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $program=`basename $0`; chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program: count the reads depth on 800bp/given window based on BWA mapped BAM file, need to install samtools software

Usage: $program <refdb> <BAM>

	-k1	window_size(bp)
		default = 800, keep the same with fasta.db.window size

	-k2	read_size, used for correction as putting most part of the read onto the window
		default = 100
		
USAGE

my %opts;
GetOptions(\%opts, "k1:i","k2:i","help!");
die $usage if ( @ARGV != 2 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#

$opts{k1} = (defined $opts{k1}) ? $opts{k1} : 800;
$opts{k2} = (defined $opts{k2}) ? $opts{k2} : 100;
=pod
#aal,HWI-D00621:25:HAJ30ADXX:2:2216:19030:101341 99      1       237420928       60      80M     =       237421248       388     ATCTTTACTGGCATAAATACTAGTCTATATTTTTAAGATGTAAAAAATAATAGAAATTTCTGTTAAAATTATGTAGGTCA        +:1BDBB?DCDF?F<<AEEIHIIDBHGGBHCEHII9:E?DF@:?C9@@B?FC>>?8?B8B4BF@8/.)/.8)8=7=C=)=     RG:Z:aal        XT:A:U  NM:i:0  SM:i:37 AM:i:37 X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:80
#aal,HWI-D00621:25:HAJ30ADXX:2:2216:19030:101341 147     1       237421248       60      12S68M  =       237420928       -388    CATCCAGTCTGTTCTCTGCCTAGAGTGTCATTTTTTTTCCTTGTCTTGTCTGGCTAACTCATTAGTCTTTTAATTCTCAG        ############@@D?;@=7)=).).).);A76IIIEDB?9ED?EE@EFA?@FEBEC<,?<+<EIEEBDDC=:B;D????     RG:Z:aal        XC:i:68 XT:A:U  NM:i:0  SM:i:37 AM:i:37 X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:68
#HWI-ST807:180:D0HW0ACXX:2:2101:3488:47705       99      OAR13   62503587        60      101M    =       62503620

#HF>D:B<BF=CCEHIBFEHEH@4=7=3?@@?>D>AAB6;=?@:   X0:i:1  X1:i:0  MD:Z:101        RG:Z:0  XG:i:0  AM:i:37 NM:i:0  SM:i:3
#JJJJJIJIJIIJHHHHHHFFFFFFEEEDECDDDDDDDDDDDDA   X0:i:1  X1:i:2  XA:Z:OAR8,-64727766,101M,3;OAR1,-71312757,101M,3;
=cut

my $window_size = $opts{k1};
my $step_size = $opts{k1}/2;
my $read_size = $opts{k2}/2;
my %uniq;
my %mul;
print STDERR "stats read depth per window start!\n";
print STDERR scalar localtime,"\n";	
##########################BAM file loading and reads depth stat###############################################################
my $header;
open (HEAD, "-|", "samtools view -H $ARGV[1]") or die "$!\n";
while (<HEAD>){
	chomp;
	next unless /^\@RG/;
	($header) = $_ =~ /SM:(\S+)/;
}
close HEAD;
die "Unexpected bam format...\n" unless defined $header;
print STDERR "Parsing sample $header ...\n";
open (IN0, "-|", "samtools view -F 0x504 $ARGV[1]") or die "$!\n";
while (<IN0>){
	chomp;
	my @read_inf = split/\t/;
	$read_inf[3] += $read_size;
	my $start = int($read_inf[3]/$opts{k1})*$opts{k1}+1;
	my $end = int($read_inf[3]/$opts{k1}+0.5)*$opts{k1}-$step_size+1;
	if ($_ =~ /\sXA:Z/){
		$mul{$read_inf[2]}{$start}++;	
		$mul{$read_inf[2]}{$end}++;
	}
	else{
		$uniq{$read_inf[2]}{$start}++;
		$uniq{$read_inf[2]}{$end}++;
	}
}
close IN0;
print STDERR "$ARGV[1] file loading done!";
print STDERR scalar localtime,"\n";
##############################################################################################################################
open (REFDB, "$ARGV[0]") or die "Failed to open fasta.db file!\n";
open (OUT, ">$header\_raw") or die "raw file writing failed\n";
print OUT "#chr\tindex\tpos\tuniq\tmul\tGC\trepeat\n";
while (<REFDB>){
	chomp;
	my @genome_inf = split/\t/;
	if (exists $uniq{$genome_inf[0]}{$genome_inf[2]} and exists $mul{$genome_inf[0]}{$genome_inf[2]}){
		print OUT $genome_inf[0],"\t",$genome_inf[1],"\t",$genome_inf[2],"\t",$uniq{$genome_inf[0]}{$genome_inf[2]},"\t",$mul{$genome_inf[0]}{$genome_inf[2]},"\t",$genome_inf[3],"\t",$genome_inf[4],"\n";
	}
	elsif (exists $uniq{$genome_inf[0]}{$genome_inf[2]} and not exists $mul{$genome_inf[0]}{$genome_inf[2]}){
		print OUT $genome_inf[0],"\t",$genome_inf[1],"\t",$genome_inf[2],"\t",$uniq{$genome_inf[0]}{$genome_inf[2]},"\t",0,"\t",$genome_inf[3],"\t",$genome_inf[4],"\n";
	}
	elsif (not exists $uniq{$genome_inf[0]}{$genome_inf[2]} and exists $mul{$genome_inf[0]}{$genome_inf[2]}){
		print OUT $genome_inf[0],"\t",$genome_inf[1],"\t",$genome_inf[2],"\t",0,"\t",$mul{$genome_inf[0]}{$genome_inf[2]},"\t",$genome_inf[3],"\t",$genome_inf[4],"\n";
	}
	else{
		print OUT $genome_inf[0],"\t",$genome_inf[1],"\t",$genome_inf[2],"\t",0,"\t",0,"\t",$genome_inf[3],"\t",$genome_inf[4],"\n";
	}
}
close REFDB;
close OUT;
print STDERR "completed bam parsing!\n";
