#!/usr/bin/perl -w
use strict;
die "Usage : perl $0 <PSL> <window_size> <output>\n" unless @ARGV == 3;
open (IN, "$ARGV[0]") or die "blat PSL file required!\n";
open (OUT, ">$ARGV[2]") or die "permission denied!\n";
my $step_size = $ARGV[1]/2;
my $kmer_size = $ARGV[1]/2;
my %hash;
while(<IN>){
	next unless /^\d+/;
	chomp;
	my @readinf = split/\t/, $_;
	next if $readinf[7] >= 5;
	my @read = split/_/, $readinf[9];
	my $pos = pop @read;
	my $chr = join("_", @read);
	my $from = int(($pos+1)/$step_size)*$step_size+1;
	$readinf[15] += $kmer_size/2;
	my $mapped1 = int($readinf[15]/$ARGV[1])*$ARGV[1]+1;
	my $mapped2 = int($readinf[15]/$ARGV[1]+0.5)*$ARGV[1]-$step_size+1;
	next if ($mapped1 < 0 or $mapped2 < 0);
	$hash{$chr."-".$from}{$readinf[13]}{$mapped1} = 1;
	$hash{$chr."-".$from}{$readinf[13]}{$mapped2} = 1;
}
close IN;
my %link;
foreach my $from_chr (sort keys %hash){
	foreach my $map_chr (sort keys %{$hash{$from_chr}}){
		my @map_pos = sort {$a <=> $b} keys %{$hash{$from_chr}{$map_chr}};
		for (my $index = 0; $index < $#map_pos; $index++){
			if ($map_pos[$index+1]-$map_pos[$index] == $step_size){
				$link{$from_chr}{$map_chr."-".$map_pos[$index]} = 1;
				$link{$from_chr}{$map_chr."-".$map_pos[$index+1]} = 1;
			}
		}
	}
}
foreach my $fromchr(sort keys %link){
	next if scalar(keys %{$link{$fromchr}}) <= 2 or scalar(keys %{$link{$fromchr}}) > 20;
	next unless scalar(keys %{$link{$fromchr}})%2 eq 0;
	print OUT $fromchr,"\t",join("\t", keys %{$link{$fromchr}}),"\n";
}
close OUT;
