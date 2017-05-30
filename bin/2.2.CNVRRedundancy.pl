#!/usr/bin/perl -w
use Data::Dumper;
use strict;
die "Usage: $0 <CNVR> <correlation> <output>\n" unless @ARGV == 3;
open (IN, "$ARGV[0]") or die "population CNVR file required!\n";
open (OUT,">$ARGV[2]") or die "permission denied!\n";
my $percent = 0.1; #the percent for the distance of two adjacent CNV in their CNV length when merge
my $correlation = $ARGV[1]; # the pearson correlation of copy numbers between two adjacent CNV when merge
#############################################input format####################################
=pod
chr	start	end	number	duplication	repeat	gc	abnormal ID1 ID2 ... AVERAGE SD
=cut
my %hash;   #chr -> $start -> @array
my $header = <IN>;
print OUT $header;
my $firstline = <IN>;
chomp $firstline;
my @tmp_start_array = split/\s+/,$firstline;
my $tmp_chr = shift  @tmp_start_array;
my $tmp_start = shift @tmp_start_array;
my $tmp_end = shift @tmp_start_array;
my $tmp_effective_windows = shift @tmp_start_array;
my $tmp_length = $tmp_end - $tmp_start + 1;
my @tmp_weight_array1 = ($tmp_effective_windows, @tmp_start_array);
my $tmp_gap = shift @tmp_start_array;
my $tmp_repeat = shift @tmp_start_array;
my $tmp_gc = shift @tmp_start_array;
my $tmp_kmer = shift @tmp_start_array;
my $tmp_sd = pop @tmp_start_array;
my $tmp_average = pop @tmp_start_array;
$hash{$tmp_chr}{$tmp_start} = [$tmp_end,@tmp_weight_array1];
##############################################################################################
while (<IN>) {
	chomp;
	my @tmp = split/\s+/,$_;
	my $chr = shift @tmp;
	my $start = shift @tmp;
	my $end = shift @tmp;
	my $effective_windows = shift @tmp;
	my $length = $end - $start + 1;
	my @tmp_weight_array2 = ($effective_windows, @tmp);
	my $gap = shift @tmp;
	my $repeat = shift @tmp;
	my $gc = shift @tmp;
	my $kmer = shift @tmp;
	my $sd = pop @tmp;
	my $average = pop @tmp;
	die "unexpected sample number!\n" unless length(@tmp_start_array) == length(@tmp);
	my $tmp_cor = &pearson(\@tmp,\@tmp_start_array);
#	print STDERR "$tmp_cor\n";
	if ($tmp_chr eq $chr && (($start-$tmp_end+1) <= ($tmp_length + $length)*$percent or $start-$tmp_end+1 <= 300) && $tmp_cor >= $correlation){
		my @tmp_weight = &weight_average(\@tmp_weight_array1, \@tmp_weight_array2);
		$hash{$tmp_chr}{$tmp_start} = [$end,@tmp_weight];
		@tmp_weight_array1 = @tmp_weight;
		shift @tmp_weight;
		shift @tmp_weight;
		shift @tmp_weight;
		shift @tmp_weight;
		shift @tmp_weight;
		pop @tmp_weight;
		pop @tmp_weight;
		@tmp_start_array = @tmp_weight;
		$tmp_chr = $chr;
		$tmp_end = $end;
		$tmp_length = $end - $tmp_start + 1;
	}else{
		undef @tmp_weight_array1;
		$hash{$chr}{$start} = [$end,@tmp_weight_array2];
		@tmp_weight_array1 = @tmp_weight_array2;
		@tmp_start_array = @tmp;
		$tmp_chr = $chr;
		$tmp_end = $end;
		$tmp_start = $start;
		$tmp_length = $start - $end + 1;
	}
}
foreach my $chr(sort keys %hash){
	foreach my $pos_start (sort {$a <=> $b} keys %{$hash{$chr}}){
		print OUT $chr,"\t",$pos_start,"\t",join("\t",@{$hash{$chr}{$pos_start}}),"\n";
	}
}
############################################subroutine######################
sub pearson {
	my ($ref_a, $ref_b) = @_;
	my @x = @{$ref_a};
	my @y = @{$ref_b};
	if($#x == $#y){
		my $N = $#x;
		my $sum_sq_x = 0;
		my $sum_sq_y = 0;
		my $sum_coproduct = 0;
		my $mean_x = &get_average(@x);
		my $mean_y = &get_average(@y);
		for(my $i=0;$i<=$N;$i++){
			my $delta_x = $x[$i] - $mean_x;
			my $delta_y = $y[$i] - $mean_y;
			$sum_sq_x += $delta_x * $delta_x;
			$sum_sq_y += $delta_y * $delta_y;
			$sum_coproduct += $delta_x * $delta_y;
		}
		my $pop_sd_x = sqrt($sum_sq_x);
		my $pop_sd_y = sqrt($sum_sq_y);
		my $cov_x_y = $sum_coproduct;
		my $correlation = $pop_sd_x * $pop_sd_y > 0 ? $cov_x_y / ($pop_sd_x * $pop_sd_y) : 0;
		return $correlation;
	}
}
sub weight_average {
	my ($array1, $array2) = @_;
	my @weight;
	my $effective_windows_array1 = shift @{$array1};
	my $effective_windows_array2 = shift @{$array2};
	my $sum_effective_windows = $effective_windows_array1+$effective_windows_array2;
	push @weight, $sum_effective_windows;
	for my $i (0..$#{$array1}){
		push @weight, $sum_effective_windows ? sprintf "%.2f",($array1->[$i]*$effective_windows_array1+$array2->[$i]*$effective_windows_array2)/$sum_effective_windows : 0;
	}
	return @weight;
}
sub get_average {
	my @array = @_;
	my $sum;
	for my $value (@array) {
		$sum += $value;
	}
	my $average = @array > 0 ? $sum/@array : 0;
	return $average;
}
