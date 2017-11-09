#!/usr/bin/perl -w
use strict;
my $program=`basename $0`; chomp $program;
die "Program: Sum the raw reads count from the high similarity windows record in the duplicated window record file\n\nUsage: perl $program <raw> <dup> <output>\n" unless @ARGV == 3;

open (IN, "$ARGV[0]") or die "Error: Cannot read your input raw reads count file $ARGV[0].\n";
open (LINK, "$ARGV[1]") or die "Error: Cannot read your input duplicated window record file $ARGV[1].\n";
open (OUT, ">$ARGV[2]") or die "Permission denied!\n";

my %raw;
while(<IN>){
	next if /^#/;
	chomp;
	my @raw_inf = split/\t/;
	$raw{$raw_inf[0]}{$raw_inf[2]} = $raw_inf[3]+$raw_inf[4];
}
print STDERR "raw reads count file loading finished!\n";
my %link; #duplicated windows
while(<LINK>){
	chomp;
	my @link_inf = split/\s+/;
	my ($chr, $pos) = split/:/, shift @link_inf;
	$link{$chr}{$pos} = \@link_inf;
}
close LINK;
print STDERR "duplicated window record file loading finished!\n";
seek IN, 0, 0;
my $correct = 0;
while(<IN>){
	next if /^#/;
	chomp;
	my @raw_inf = split/\t/;
	my $absolute_cp = 0;
	my $genome_copy = 0;
	if (exists $link{$raw_inf[0]}{$raw_inf[2]}){
		$correct++;
		foreach my $pos (@{$link{$raw_inf[0]}{$raw_inf[2]}}){
			$genome_copy += 1;
			my ($tmp_chr, $tmp_pos) = split/:/, $pos;
			if (exists $raw{$tmp_chr}{$tmp_pos}){
				$absolute_cp += $raw{$tmp_chr}{$tmp_pos};
			}
			else{
				$absolute_cp += 0;
			}
		}
		print OUT "$raw_inf[0]\t$raw_inf[1]\t",$absolute_cp,"\t$raw_inf[5]\t$raw_inf[6]\t",$genome_copy,"\n";
	}
	else{
		print OUT "$raw_inf[0]\t$raw_inf[1]\t",$raw_inf[3]+$raw_inf[4],"\t$raw_inf[5]\t$raw_inf[6]\t",1,"\n";
	}
}
close IN;
close OUT;
print STDERR "corrected window number: $correct\n";
print STDERR "absolute reads counts correction finished!\n";
