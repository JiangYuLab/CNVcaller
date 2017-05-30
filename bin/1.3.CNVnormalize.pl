#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
my $program=`basename $0`; chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program: GC correction and RD normalization

Usage: $program <corrected raw>

	-k1	window_size
		default = 800, keep the same with fasta.db.window size

	-k2	sex_chr_name[X,chrX,Z,chrZ...]; known half coverage chromosome/region, as chrX in mammals and chrZ in birds
		default = X

	-k3	percent of gain[5-20]%, pre-estimation ratio of copy gain region and errs in genome
		default = 5

	-k4	percent of loss[5-20]%, pre-estimation ratio of copy loss region and errs in genome
		default = 5

USAGE

my %opts;
my @sex = ("X");
GetOptions(\%opts, "k1:i","k2:s@","k3:i","k4:i","k5:f","help!");
die $usage if ( @ARGV<1 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#

$opts{k1} = (defined $opts{k1}) ? $opts{k1} : 800;
$opts{k2} = (defined $opts{k2}) ? $opts{k2} : \@sex;
$opts{k3} = (defined $opts{k3}) ? $opts{k3} : 5;
$opts{k4} = (defined $opts{k4}) ? $opts{k4} : 5;
##############################################################################################################################
my $repeat_max=$opts{k1}*0.3;
my %GC_region_average;
my @global_array=();my @X_array=();
open (IN, "$ARGV[0]") or die "Failed to open raw file!\n";
my $header = (split/_link/, (split/\//,$ARGV[0])[-1])[0];
print STDERR "parsing $header ......\n";
my %genome_copy;
while(<IN>){
	chomp;
	my @RAW = split/\t/;
	$genome_copy{$RAW[0]}{$RAW[1]} = $RAW[-1];
	if (grep (/^$RAW[0]$/ ,@{$opts{k2}})){
		if ($RAW[4] <= $repeat_max and $RAW[-1] eq '1') {
			push @X_array, $RAW[2];
		}
	}
	else{
		if ($RAW[4] <= $repeat_max and $RAW[-1] eq '1' ) {
			push @{$GC_region_average{$RAW[3]}}, $RAW[2];
			push @global_array, $RAW[2];
		}
	}
}
close IN;
print STDERR "$ARGV[0] file loading finished ...\n";
##############################################################################################################################
@global_array = sort {$a <=> $b} @global_array;
my $median90 = $global_array[int(0.9*@global_array)];
my $median50 = $global_array[int(0.5*@global_array)];
my $median10 = $global_array[int(0.1*@global_array)];
my $median8 = $global_array[int(0.08*@global_array)];
my $median6 = $global_array[int(0.06*@global_array)];
my $median4 = $global_array[int(0.04*@global_array)];
my $median2 = $global_array[int(0.02*@global_array)];
print STDERR "tail 2 percent window is lower than $median2\n";
print STDERR "tail 4 percent window is lower than $median4\n";
print STDERR "tail 6 percent window is lower than $median6\n";
print STDERR "tail 8 percent window is lower than $median8\n";
print STDERR "tail 10 percent window is lower than $median10\n";
print STDERR "tail 50 percent window is lower than $median50\n";
print STDERR "tail 90 percent window is lower than $median90\n";
my $Xmedian50 = 0;
if (@X_array >= 1){
	@X_array = sort {$a <=> $b} @X_array;
	my $Xmedian90= $X_array[int(0.9*@X_array)];
	$Xmedian50= $X_array[int(0.5*@X_array)];
	my $Xmedian10= $X_array[int(0.1*@X_array)];
	print STDERR "tail 10 percent ",join(" ", @{$opts{k2}}), " window is lower than $Xmedian10\n";
	print STDERR "tail 50 percent ",join(" ", @{$opts{k2}}), " window is lower than $Xmedian50\n";
	print STDERR "tail 90 percent ",join(" ", @{$opts{k2}}), " window is lower than $Xmedian90\n";
}

my $sex_correct_fold;
if ($Xmedian50>$median50*0.75 && $Xmedian50<$median50*1.5 ){
	print STDERR "sex chromosome ",join(" ", @{$opts{k2}}), " show similar coverage of autosome!!!\n\n";
	$sex_correct_fold=1;
}
elsif($Xmedian50<$median50*0.75 && $Xmedian50>$median50*0.25){
	print STDERR "sex chromosome ",join(" ", @{$opts{k2}}), " show ~ half coverage of autosome!!!\n\n";
	$sex_correct_fold=2;
}
else {
	$sex_correct_fold=1;
	print STDERR "sex chromosome ",join(" ", @{$opts{k2}}), " show UNKOWN relationship with autosome???\n\n";
}

print STDERR "raw RD file printing done!";
print STDERR scalar localtime,"\n";
@global_array=();@X_array=();
####################################calculate average value for each GC content region #######################################
my %region_average_sd;
my @sorted_tmp;
my $sum;
my $sum_window;
foreach my $tmp_gc (sort {$a <=> $b} keys %GC_region_average){
	$sum_window = 0;
	$sum = 0;
	@sorted_tmp = ();
	@sorted_tmp = sort {$a <=> $b} @{$GC_region_average{$tmp_gc}};
	for my $i (int(($opts{k4}/100)*@sorted_tmp)..int((1-$opts{k3}/100)*@sorted_tmp)) {
		$sum_window += 1;
		$sum += $sorted_tmp[$i];
	}
	$region_average_sd{$tmp_gc} = $sum_window > 0 ? sprintf "%.2f", ($sum/$sum_window) : 0 ;
	print STDERR "$tmp_gc\t$region_average_sd{$tmp_gc}\t$sum_window\t$#sorted_tmp\n";######################################GC group output test
}
%GC_region_average=();

print STDERR "calculate average value for each GC content region done!\n";
print STDERR scalar localtime,"\n";
print STDERR "correct read depth ...\n";
#######################################correct read depth according to GC content, into 40% average GC########################
my $standard_average = $region_average_sd{$opts{k1}*0.4};
foreach my $tmp_gc (sort {$a <=> $b} keys %region_average_sd){
	$region_average_sd{$tmp_gc} = $region_average_sd{$tmp_gc} > 0 ? sprintf "%.2f",($standard_average/$region_average_sd{$tmp_gc}) : 0;
#	print STDERR "$tmp_gc\t$region_average_sd{$tmp_gc}\n";
}
my %clean_record;
open (RAW, "$ARGV[0]") or die "Failed to open raw file!\n";
while (<RAW>){
	chomp;
	my @read_inf = split/\t/;
	$region_average_sd{$read_inf[3]} = 0 if not exists $region_average_sd{$read_inf[3]};
	if (grep (/^$read_inf[0]$/ ,@{$opts{k2}})) {
		$read_inf[2]=sprintf "%.2f",($sex_correct_fold*$read_inf[2]*$region_average_sd{$read_inf[3]});
	}
	else{
		$read_inf[2]=sprintf "%.2f",($read_inf[2]*$region_average_sd{$read_inf[3]});
		push @global_array, $read_inf[2];
	}
	$clean_record{$read_inf[0]}{$read_inf[1]}=$read_inf[2];
}
close RAW;
print STDERR "correct read count according to GC content done!\n";
print STDERR scalar localtime,"\n";
@global_array = sort {$a <=> $b} @global_array;
my $correct_median50 = sprintf "%.2f", $global_array[int(0.5*@global_array)];
print STDERR "global median50 $correct_median50\n";
##########################################calculate the global average and SD#################################################
$opts{k3}=100-$opts{k3};
my $global_max = &top_percent($opts{k3}, @global_array);
my $global_min = &top_percent($opts{k4}, @global_array);
my ($global_average, $global_sd);
@global_array=();
foreach my $tmp_chr (keys %clean_record){
	foreach my $tmp_pos(keys %{$clean_record{$tmp_chr}}){
		next if ($clean_record{$tmp_chr}{$tmp_pos} > $global_max || $clean_record{$tmp_chr}{$tmp_pos} < $global_min);
		push @global_array, $clean_record{$tmp_chr}{$tmp_pos};
	}
}
	$sum=0;
	for my $value (@global_array) {
		$sum += $value;
	}
	$global_average =$sum/@global_array;
	$sum=0;
	for my $value (@global_array) {
		$sum += ($value-$global_average)**2;
	}
	$global_sd =sqrt($sum/(@global_array-1));
	
	$global_average = sprintf "%.2f",$global_average;
	$global_sd = sprintf "%.2f",$global_sd;
print STDERR "mapped reads stats across the genome-wide windows:\n";
print STDERR "global average: $global_average\tglobal SD: $global_sd\n";
print STDERR "$opts{k3} percent window is more than $global_max\n";
print STDERR "$opts{k4} percent window is lower than $global_min\n";
print STDERR scalar localtime,"\n";
####################################################correct window and candiate CNVR##########################################
open (OUT2, ">$header\_mean_$correct_median50\_SD_$global_sd\_sex_$sex_correct_fold") or die "permission denied!\n";
foreach my $tmp_chr (sort keys %clean_record){
	foreach my $i (sort {$a <=> $b} keys %{$clean_record{$tmp_chr}}){
		my $copynumber = sprintf "%.2f",$clean_record{$tmp_chr}{$i}/$correct_median50;
		print OUT2 $tmp_chr,"\t",$i,"\t",$copynumber,"\t",$genome_copy{$tmp_chr}{$i},"\n";
	}
}
close OUT2;
print STDERR "output correct window copy done!\n";
print STDERR scalar localtime,"\n";
##########################################subroutine##########################################################################
sub top_percent {
	my $cutoff;
	my ($percent,@tmp) = @_;
	my @sorted_tmp = sort {$a <=> $b} @tmp;
	my $int_cutoff = int($percent/100*(@sorted_tmp-1));
	$cutoff = $sorted_tmp[$int_cutoff];
	return $cutoff;
}
