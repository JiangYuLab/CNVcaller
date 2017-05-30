#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
my $program=basename($0);
#variable
my $infile = '';
my $mclust = '';
my $out = '';
my $number = '';
#usage
my $usage=<<USAGE; #******* Instruction of this program *********#
Program:
Usage: $program
  -gc  genotype file by constrained Gaussian mixture mode
  -gu  genotype file by unsupervised Gaussian mixture mode
  -n   sample number
  -o   output
USAGE
GetOptions(
	'gc|constrained=s' => \$infile,
	'gu|mclust=s' => \$mclust,
	'n|number=i' => \$number,
	'o|out=s'    => \$out
) or pod2usage(2);
pod2usage("Error: Cannot read your input genotype file by constrained Gaussian mixture mode '$infile'...\n$usage") unless( -r $infile );
pod2usage("Error: Cannot read your input genotype file by unsupervised Gaussian mixture mode '$mclust'...\n$usage") unless( -r $mclust );
pod2usage("Error: You must provide your sample number '$number'...\n$usage") unless( defined $number );
pod2usage("Error: Cannot write out file '$out'...\n$usage") unless( defined $out );
open (IN, "$infile");
open (OUT, ">$out");
print OUT "##fileformat=VCFv4.1\n##CNVcaller\n##ALT=<ID=DEL,Description=\"Deletion\">\n##ALT=<ID=DUP,Description=\"Duplication\">\n##ALT=<ID=mCNV,Description=\"multi-allelic copy number variations\">\n##FORMAT=<ID=CP,Number=1,Type=Float,Description=\"Normalized copy number\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
my $header = <IN>;
my @tmp_header = split/\s+/, $header;
my $sample_start = 8+$number+2;
my @sample_name = splice(@tmp_header, $sample_start, $number);
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach my $tmp (@sample_name){
	my @sample_inf = split/\./, $tmp;
	print OUT "\t",$sample_inf[0]
}
print OUT "\n";
my %gen;
while(<IN>){
	chomp;
	my @tmp = split/\t/;
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $end = $tmp[2];
	my $average = $tmp[8+$number];
	my $type;
	my %genotype;
	next unless $average <= 2.00;
	if ($average >= 0 and $average < 1.00){
		$type = "DEL";
		$genotype{dd} = "1/1";
		$genotype{Ad} = "0/1";
		$genotype{AA} = "0/0";
	}
	elsif ($average >= 1.00 and $average <= 2.00){
		$type = "DUP";
		$genotype{AA} = "0/0";
		$genotype{AB} = "0/1";
		$genotype{BB} = "1/1";
	}
	my @genotype = splice(@tmp, $sample_start, $number);
	my @tmpgenotype = ();
	foreach my $tmp_index (0..$#genotype){
		push @tmpgenotype, "$genotype{$genotype[$tmp_index]}:$tmp[8+$tmp_index]";
	}
	$gen{$chr}{$start} = ".\t<$type>\t.\tPASS\tEND=$end;SVTYPE=$type\tGT:CP\t".join("\t", @tmpgenotype);
}
open (MCLUST, "$mclust");
<MCLUST>;
while(<MCLUST>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $chr = $tmp[0];
	my $start = $tmp[1];
	my $end = $tmp[2];
	my $average = $tmp[8+$number];
	if ( $average > 2 ){
		my $class_number = 1;
		my $class1_average = $tmp[-9];
		my $class2_average = $tmp[-8];
		my $class3_average = $tmp[-7];
		if ($class1_average eq "NA" and $class2_average eq "NA" and $class3_average eq "NA"){
			die "error! Maybe this classfication is impossible!\n";
		}
		my $class1_num = 0;
		my $class2_num = 0;
		my $class3_num = 0;
		my @classficaton = splice(@tmp, $sample_start, $number);
		foreach my $tmp_classfication (@classficaton){
			$class1_num++ if $tmp_classfication == 1;
			$class2_num++ if $tmp_classfication == 2;
			$class3_num++ if $tmp_classfication == 3;
		}
		die "sample number error\n" if $class1_num+$class2_num+$class3_num != $number;
		my %genotype;
		$genotype{1} = "0/0";
		$genotype{2} = "0/1";
		$genotype{3} = "1/1";
		print OUT $chr,"\t",$start,"\t",$start,"\t",".","\t","<mCNV>","\t.\t","PASS","\t","END=$end;SVTYPE=mCNV","\t","GT:CP";
		foreach my $tmp_index (0..$#classficaton){
			print OUT "\t","$genotype{$classficaton[$tmp_index]}:$tmp[8+$tmp_index]";
		}
		print OUT "\n";
	}
	else{
		print OUT $chr,"\t",$start,"\t",$start,"\t",$gen{$chr}{$start},"\n";
	}
}
close MCLUST;
close OUT;
