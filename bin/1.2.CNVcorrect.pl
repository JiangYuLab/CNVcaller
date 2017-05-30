#!/usr/bin/perl -w
use strict;
use Data::Dumper;
die "Program: Sum the raw reads count from the high similarity windows record in the duplicated window record file\nUsage : perl $0 <raw> <dup> <output>\n" unless @ARGV == 3;
open (IN, "$ARGV[0]") or die "raw file required\n";
open (LINK, "$ARGV[1]") or die "duplicated file required!\n";
open (OUT, ">$ARGV[2]\_link") or die "permission denied!\n";
my %raw;
while(<IN>){
	next if /^#/;
	chomp;
	my @raw_inf = split/\t/;
	$raw{$raw_inf[0]}{$raw_inf[2]} = $raw_inf[3]+$raw_inf[4];
}
print STDERR "raw file loading finished!\n";
my %link; #duplicated windows
my $candidate_correct = 0;
while(<LINK>){
	chomp;
	my @link_inf = split/\s+/;
	my ($chr,$pos) = split/-/, shift @link_inf;
	next unless ($#link_inf+1)%2 eq 0;
	$candidate_correct++;
	$link{$chr}{$pos} = \@link_inf;
}
close LINK;
print STDERR "candidate correct window $candidate_correct\n";
print STDERR "dup file loading finished!\n";
seek IN, 0, 0;
my $correct = 0;
while(<IN>){
	next if /^#/;
	chomp;
	my @raw_inf = split/\t/;
	my $absolute_cp = 0;
	my $genome_copy = 0;
	if (exists $link{$raw_inf[0]}{$raw_inf[2]}){
#		print STDERR $raw_inf[0],"\t",$raw_inf[2],"\n";
		$correct++;
		foreach my $pos (@{$link{$raw_inf[0]}{$raw_inf[2]}}){
			$genome_copy += 0.5;
			my ($tmp_chr,$tmp_pos) = split/-/, $pos;
			if (exists $raw{$tmp_chr}{$tmp_pos}){
				$absolute_cp += $raw{$tmp_chr}{$tmp_pos}/2;
			}
			else{
				$absolute_cp += 0;
			}
		}
		print OUT "$raw_inf[0]\t$raw_inf[1]\t",$absolute_cp,"\t$raw_inf[5]\t$raw_inf[6]\t",$genome_copy,"\n";
	}
	else{
		print OUT "$raw_inf[0]\t$raw_inf[1]\t",$raw_inf[3]+$raw_inf[4],"\t$raw_inf[5]\t$raw_inf[6]\t",1+$genome_copy,"\n";
	}
}
close IN;
close OUT;
print STDERR "correct window $correct\n";
print STDERR "read depth refine finished!\n";
