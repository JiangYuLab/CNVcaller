#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <ref> <window_size> <output>\n" unless @ARGV == 3;
my $storefile = &FastaReader($ARGV[0]) or die "fasta file is required!\n";
print STDERR "ref loading finished!\n";
my $length = $ARGV[1]/2;
my $step = $ARGV[1]/2;
open (OUT, ">$ARGV[2]") or die "permission denied!\n";
foreach my $chr (sort {$a cmp $b} keys %{$storefile}){
	if (length($storefile->{$chr}) <= $length){
		print OUT ">$chr","_","0","\n",$storefile->{$chr},"\n";
	}
	else{
		for (my $pos = 0; $pos <= int(length($storefile->{$chr})/$step); $pos++){
			my $start = $pos*$step;
			my $seq = substr($storefile->{$chr}, $start, $length);
			if (length($seq) eq $length){
				print OUT ">$chr","_","$start","\n",$seq,"\n";
			}
			else{
				next unless length($seq) > 0;
				print OUT ">$chr","_","$start","\n",$seq,"\n";
			}
		}
	}
	print STDERR "$chr completed!\n";
}
close OUT;
sub FastaReader {
	my ($file) = @_;
	open IN, "<", $file or die "Fail to open file: $file!\n";
	local $/ = '>';
	<IN>;
	my ($head, $seq, %hash);
	while (<IN>){
		s/\r?\n>?$//;
		( $head, $seq ) = split /\r?\n/, $_, 2;
		my $tmp = (split/\s+/,$head)[0];
		$seq =~ s/\s+//g;
		$hash{$tmp} = $seq;
	}
	close IN;
	$/ = "\n";
	return \%hash;
}
