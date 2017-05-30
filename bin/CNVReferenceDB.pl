#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $program=`basename $0`; chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program: format reference genome into state table of GC and repeat content per window, as removing gapped regions.

Usage: $program <ref>
		
  -w  window_size
      default = 800

  -l  GC content lower limit
      default = 0.2

  -u  GC content upper limit
      default = 0.7

  -g  gap content
      default = 0.5			
USAGE
my $help;
my %opts = (w=>800, l=>0.2, u=>0.7, g=>0.5);
GetOptions(\%opts, "w:i","l:f","u:f","g:f","help!");
die $usage if ( @ARGV!=1 || defined($opts{"help"}));
#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#
my $window_size=$opts{w};
my $step_size=$opts{w}/2;
my $GC_count_max=$window_size*$opts{u};
my $GC_count_min=$window_size*$opts{l};
print STDERR "Window_size: $opts{w}\n";
print STDERR "GC lower limit :$GC_count_min\nGC upper limit :$GC_count_max\n";
print STDERR "genome format start!\n";
print STDERR scalar localtime,"\n";
#######################calculate GC content, repeat content for genome############
my $reffile = &FastaReader($ARGV[0]) or die "ref fasta file required!\n";
open (OUT, ">referenceDB.$opts{w}") or die "$!\n";
my $number;
my $window_number;
foreach my $chr (sort keys %{$reffile}){
	$window_number = int(length($reffile->{$chr})/$step_size);
	my $i = 0;
	for $number (0..$window_number){
		my $position=$number*$step_size;
		my $sub_seq = substr($reffile->{$chr}, $position, $window_size);
		my $gap_count = () = $sub_seq =~ /N/gi;
		my $gap_content = $window_size > 0 ? sprintf "%.2f",$gap_count/$window_size : 0;
		next if $gap_content > $opts{g};
		my $GC_count = () = $sub_seq =~ /[GC]/gi;
		next if ($GC_count > $GC_count_max || $GC_count < $GC_count_min);
		my $repeat_count = () = $sub_seq =~ /[atgc]/g ;
		$i++;############used window NO. per chromosome.
		print OUT "$chr\t$i\t",$position+1,"\t$GC_count\t$repeat_count\t$gap_content\n";
	}
}
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
close OUT;
local $/ = "\n";
print STDERR "genome format done!\n";
print STDERR scalar localtime,"\n";
