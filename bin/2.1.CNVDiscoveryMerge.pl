#!/usr/bin/perl -w
use strict;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Getopt::Long;
my $program=`basename $0`; chomp $program;
my $usage=<<USAGE; #******* Instruction of this program *********#

Program: merge output RD_corr files into one file based on refgenome.fasta.db file, then output CNVR.

Usage: $program <refdb> <normalizedFileList> <excludedFileList> <primaryCNVR>
	-k1	window size
		default = 800

	-k2	percent of gain/loss samples in all samples per window samples
		default = 0.1
			
	-k3	number of homozygous gain/loss samples in all samples per window [2 - 5]
		default = 3
			
	-k4	homozygous gain/loss copy number cutoff [>1.75 or <0.25]
		default = 0.75

	-k5	CNV copy number hard cutoff
		default = 0.35

	-k6	pearson correlation R value, Cnvcaller was set according to the following criterion:
		0.5  for sample size (0, 30]
		0.4  for sample size (30, 50]
		0.3  for sample size (50, 100]
		0.2  for sample size (100, 200]
		0.15 for sample size (200, 500]
		0.1  for sample size (500,+∞)
		You can also set this value by yourself.
USAGE
my %opts;
GetOptions(\%opts, "k1:i","k2:f","k3:i","k4:f","k5:f","k6:f","help!");
die $usage if ( @ARGV!=4 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#

$opts{k1} = (defined $opts{k1}) ? $opts{k1} : 800;
$opts{k2} = (defined $opts{k2}) ? $opts{k2} : 0.05;
$opts{k3} = (defined $opts{k3}) ? $opts{k3} : 3;
$opts{k4} = (defined $opts{k4}) ? $opts{k4} : 0.75;
$opts{k5} = (defined $opts{k5}) ? $opts{k5} : 0.35;
######################################################################################################################################
my %filter_sample;
open (Filter, "$ARGV[2]") or die "corrected/raw filter sample removed!\n";
while(<Filter>){
	chomp;
	$filter_sample{$_} = 1;
}
close Filter;
##########################read all of the corrected coverage per window and +/- results from multiple files###########################
my %hash_res_gain;
my %hash_res_loss;
my %hash_res_homo_gain;
my %hash_res_homo_loss;
my %hash_cnv_num;
my $sample_num = 0;
my $effective_sample = 0;
######################################################################################################################################
open (IN, "$ARGV[1]") or die "corr file list required\n";
my @sv_array;
while(<IN>){
	chomp;
	next if (/^\s+$/);
	next if (exists $filter_sample{$_});
	my ($average, $sd) = $_ =~ /(\d+\.?\d+)_SD_(\d+\.?\d+)/;
	my $sv = $sd > 0 ? $opts{k5}*$average/$sd : 0;
	push @sv_array, $sv;
}
@sv_array = sort {$b <=> $a}@sv_array;
my $median75 = $sv_array[int(0.75*@sv_array)];
print STDERR $median75,"\n";
seek IN, 0, 0;
my %genome_copy;
while(<IN>){
	chomp;
	my $remove = 0;
	next if (/^\s+$/);
	$sample_num++;
	if (exists $filter_sample{$_}){$remove = 1}
	$effective_sample++ if $remove == 0;
	my ($average, $sd) = $_ =~ /(\d+\.?\d+)_SD_(\d+\.?\d+)/;
	my $sv = $sd > 0 ? $opts{k5}*$average/$sd : 0;
	my $hetcutoff = ($sv <= $median75) ? $median75*$sd/$average : $opts{k5};
	print STDERR $_,"\t",$hetcutoff,"\n";
	open (TMP, "$_") or die "$_ file not found\n";
	while(<TMP>){
		chomp;
		my @tmp = split/\t/;
		$genome_copy{$tmp[0]}{$tmp[1]} = $tmp[3] if not defined $genome_copy{$tmp[0]}{$tmp[1]}; ##genome copies according to duplicated window file
		next unless $remove == 0;
		if ($tmp[2] >= 1+$hetcutoff ) {
			$hash_res_gain{$tmp[0]}{$tmp[1]}++;
			if ($tmp[2] >= 1+$opts{k4}) {
				$hash_res_homo_gain{$tmp[0]}{$tmp[1]}++;
			}
		}elsif ($tmp[2] <= 1-$hetcutoff ){
			$hash_res_loss{$tmp[0]}{$tmp[1]}++;
			if ($tmp[2] <= 1-$opts{k4}){
				$hash_res_homo_loss{$tmp[0]}{$tmp[1]}++;
			}
		}
	}
	close TMP;
}
die "Error, none effective sample was detected!\n" if ($effective_sample == 0); 
print STDERR "sample number: $sample_num\neffective sample: $effective_sample\n";
######################################################################################################################################
if (defined $opts{k6}){
	$opts{k6} = $opts{k6};
}
elsif ($sample_num <= 30){
	$opts{k6} = 0.5;
}
elsif ($sample_num > 30 and $sample_num <= 50){
	$opts{k6} = 0.4;
}
elsif ($sample_num > 50 and $sample_num <= 100){
	$opts{k6} = 0.3;
}
elsif ($sample_num > 100 and $sample_num <= 200){
	$opts{k6} = 0.2;
}
elsif ($sample_num > 200 and $sample_num <= 500){
	$opts{k6} = 0.15;
}
elsif ($sample_num > 500){
	$opts{k6} = 0.1;
}
print STDERR "Pearson correlation: $opts{k6}\n";
########################################write cnv windows into one merged file #######################################################
my %refdb;
my %gain_effective;
my %loss_effective;
open (REFDB, "$ARGV[0]") or die "$ARGV[0] file required\n";
while(<REFDB>){
	chomp;
	next if(/^\s+$/);
	my @tmp = split/\t/;
	$refdb{$tmp[0]}{$tmp[1]} = $_;
	if ((defined $hash_res_gain{$tmp[0]}{$tmp[1]} and ($hash_res_gain{$tmp[0]}{$tmp[1]}/$effective_sample) >= $opts{k2}) || (defined $hash_res_homo_gain{$tmp[0]}{$tmp[1]} and $hash_res_homo_gain{$tmp[0]}{$tmp[1]} >= $opts{k3}) ) {
		$hash_res_gain{$tmp[0]}{$tmp[1]} = 1;
		$gain_effective{$tmp[0]}{$tmp[1]} = 1;
	}
	else{
		$hash_res_gain{$tmp[0]}{$tmp[1]} = 0;
	}
	if ((defined $hash_res_loss{$tmp[0]}{$tmp[1]} and ($hash_res_loss{$tmp[0]}{$tmp[1]}/$effective_sample) >= $opts{k2}) || (defined $hash_res_homo_loss{$tmp[0]}{$tmp[1]} and $hash_res_homo_loss{$tmp[0]}{$tmp[1]} >= $opts{k3}) ) {
		$hash_res_loss{$tmp[0]}{$tmp[1]} = 1;
		$loss_effective{$tmp[0]}{$tmp[1]} = 1;
	}
	else{
		$hash_res_loss{$tmp[0]}{$tmp[1]} = 0;
	}
}
close REFDB;
print STDERR "reference db loading finished!\n";
###########################merged cnv windows (as >=3 Yes windows in 5 windows) into Regions #########################################
my %CNV_region_gain;
my %CNV_region_loss;
foreach my $tmp_chr (sort keys %refdb){
	foreach my $i (sort {$a <=> $b} keys %{$refdb{$tmp_chr}}){
		next unless ($i > 2 and defined $hash_res_gain{$tmp_chr}{$i+1} and defined $hash_res_gain{$tmp_chr}{$i+2} and defined $hash_res_loss{$tmp_chr}{$i+1} and defined $hash_res_loss{$tmp_chr}{$i+2});
		if (($hash_res_gain{$tmp_chr}{$i-2}+$hash_res_gain{$tmp_chr}{$i-1}+$hash_res_gain{$tmp_chr}{$i}+$hash_res_gain{$tmp_chr}{$i+1}+$hash_res_gain{$tmp_chr}{$i+2}) >= 3){ 
			$CNV_region_gain{$tmp_chr}{$i} = "Y"; 
		}
		elsif (($hash_res_loss{$tmp_chr}{$i-2}+$hash_res_loss{$tmp_chr}{$i-1}+$hash_res_loss{$tmp_chr}{$i}+$hash_res_loss{$tmp_chr}{$i+1}+$hash_res_loss{$tmp_chr}{$i+2}) >= 3){
			$CNV_region_loss{$tmp_chr}{$i} = "Y";
		}
	}
}
#######################################################################################################################################
my %hash_num;
seek IN, 0, 0;
while(<IN>){
	chomp;
	next if (/^\s+$/);
	open (TMP, "$_") or die "$_ file not found\n";
	while(<TMP>){
		chomp;
		my @tmp = split/\t/;
		$hash_num{$tmp[0]}{$tmp[1]} .= "$tmp[2]\t" if (exists $CNV_region_gain{$tmp[0]}{$tmp[1]} or exists $CNV_region_loss{$tmp[0]}{$tmp[1]});
	}
}
close IN;
#######################################################################################################################################
my ($tmp_gain_start, $tmp_gain_end, $tmp_loss_start, $tmp_loss_end, $gain_index, $gain_index2, $loss_index, $loss_index2);
my %cnv_tmp_gain;
my %cnv_tmp_loss;
foreach my $tmp_chr (sort keys %refdb){
	$gain_index = $gain_index2 = $loss_index = $loss_index2 = 0;
	foreach my $i (sort {$a <=> $b} keys %{$refdb{$tmp_chr}}){
		if (exists $CNV_region_gain{$tmp_chr}{$i}){
			$gain_index++;
			if ($gain_index==1) {
				$tmp_gain_start = $i;
			}
			if (exists $gain_effective{$tmp_chr}{$i}){
				$gain_index2++;
			}
			$tmp_gain_end=$i;
		}else{
			if ($gain_index2 >= 3){
				$cnv_tmp_gain{$tmp_chr}{$tmp_gain_start} = $tmp_gain_end; ##{chr -> start window index ->  end window index}
			}
			$gain_index = $gain_index2 = 0;
		}
		if (exists $CNV_region_loss{$tmp_chr}{$i}){
			$loss_index++;
			if ($loss_index==1) {
				$tmp_loss_start = $i;
			}
			if (exists $loss_effective{$tmp_chr}{$i}){
				$loss_index2++;
			}
			$tmp_loss_end=$i;
		}else{
			if ($loss_index2 >= 3){
				$cnv_tmp_loss{$tmp_chr}{$tmp_loss_start} = $tmp_loss_end;
			}
			$loss_index = $loss_index2 = 0;
                }
	}
	if ($gain_index2 >= 3){
		$cnv_tmp_gain{$tmp_chr}{$tmp_gain_start} = $tmp_gain_end;
	}
	if ($loss_index2 >= 3){
		$cnv_tmp_loss{$tmp_chr}{$tmp_loss_start} = $tmp_loss_end;
	}
}
print STDERR "define effective window finished!\n";
###################################################CNV region refine###################################################################
my %cnv_gain;
foreach my $tmp_chr (sort keys %cnv_tmp_gain){
	my @tmp_effcetive = keys %{$gain_effective{$tmp_chr}};
	foreach my $tmp_index (sort {$a <=> $b} keys %{$cnv_tmp_gain{$tmp_chr}}){
		my @tmp_selected = sort {$a <=> $b} grep {$_ >= $tmp_index and $_ <= $cnv_tmp_gain{$tmp_chr}{$tmp_index}} @tmp_effcetive;
		foreach my $tmp (0..($#tmp_selected-2)){
			my @tmp_cn_a = split/\t/, $hash_num{$tmp_chr}{$tmp_selected[$tmp]};
			my @tmp_cn_b = split/\t/, $hash_num{$tmp_chr}{$tmp_selected[$tmp+2]};
			my $tmp_correlation = &pearson(\@tmp_cn_a, \@tmp_cn_b);
			if ($tmp_correlation >= $opts{k6}){
				$cnv_gain{$tmp_chr}{$tmp_selected[$tmp]} = 1;
				$cnv_gain{$tmp_chr}{$tmp_selected[$tmp+2]} = 1;
			}
			else{
				delete $cnv_gain{$tmp_chr}{$tmp_selected[$tmp]};
			}
		}
	}
}
my %cnv_loss;
foreach my $tmp_chr (sort keys %cnv_tmp_loss){
	my @tmp_effcetive = keys %{$loss_effective{$tmp_chr}};
	foreach my $tmp_index (sort {$a <=> $b} keys %{$cnv_tmp_loss{$tmp_chr}}){
		my @tmp_selected = sort {$a <=> $b} grep {$_ >= $tmp_index and $_ <= $cnv_tmp_loss{$tmp_chr}{$tmp_index}} @tmp_effcetive;
		foreach my $tmp (0..($#tmp_selected-2)){
			my @tmp_cn_a = split/\t/, $hash_num{$tmp_chr}{$tmp_selected[$tmp]};
			my @tmp_cn_b = split/\t/, $hash_num{$tmp_chr}{$tmp_selected[$tmp+2]};
			my $tmp_correlation = &pearson(\@tmp_cn_a, \@tmp_cn_b);
			if ($tmp_correlation >= $opts{k6}){
				$cnv_loss{$tmp_chr}{$tmp_selected[$tmp]} = 1;
				$cnv_loss{$tmp_chr}{$tmp_selected[$tmp+2]} = 1;
			}
			else{
				delete $cnv_loss{$tmp_chr}{$tmp_selected[$tmp]};
			}
		}
	}
}
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
#######################################################################################################################################
my %cnv;
foreach my $tmp_chr (sort keys %refdb){
	$gain_index = $gain_index2 = $loss_index = $loss_index2 = 0;
	foreach my $i (sort {$a <=> $b} keys %{$refdb{$tmp_chr}}){
		if (exists $cnv_gain{$tmp_chr}{$i}){
			$gain_index++;
			if ($gain_index==1) {
				$tmp_gain_start = $i;
			}
			if (exists $gain_effective{$tmp_chr}{$i}){
				$gain_index2++;
			}
			$tmp_gain_end=$i;
		}else{
			if ($gain_index2 >= 3){
				$cnv{$tmp_chr}{$tmp_gain_start} = $tmp_gain_end;
			}
			$gain_index = $gain_index2 = 0;
		}
		if (exists $cnv_loss{$tmp_chr}{$i}){
			$loss_index++;
			if ($loss_index==1) {
				$tmp_loss_start = $i;
			}
			if (exists $loss_effective{$tmp_chr}{$i}){
				$loss_index2++;
                        }
			$tmp_loss_end=$i;
		}else{
			if ($loss_index2 >= 3){
				$cnv{$tmp_chr}{$tmp_loss_start} = $tmp_loss_end;
			}
			$loss_index = $loss_index2 = 0;
		}
	}
	if ($gain_index2 >= 3){
		$cnv{$tmp_chr}{$tmp_gain_start} = $tmp_gain_end;
	}
	if ($loss_index2 >= 3){
		$cnv{$tmp_chr}{$tmp_loss_start} = $tmp_loss_end;
	}
}
print STDERR "refine effective window finished!\n";
##############################compute the GC ratio, repeat ratio, Duplication per CNV region###########################################
my %cnv_GC;
my %cnv_repeat;
my %cnv_gap;
my %cnv_all;
my %cnv_mul;
my %cnv_real_start;
my %cnv_real_end;
my %kmer_median;

foreach my $tmp_chr (keys %cnv){
	foreach my $start_index (keys %{$cnv{$tmp_chr}}){
		$cnv_real_start{$tmp_chr}{$start_index} = (split/\t/, $refdb{$tmp_chr}{$start_index})[2];
		my @tmp_median_kmer;
		foreach my $tmp_index ($start_index..$cnv{$tmp_chr}{$start_index}){
			if (exists $gain_effective{$tmp_chr}{$tmp_index} or exists $loss_effective{$tmp_chr}{$tmp_index}) {
				$cnv_GC{$tmp_chr}{$start_index} += (split/\t/, $refdb{$tmp_chr}{$tmp_index})[3];
				$cnv_repeat{$tmp_chr}{$start_index} += (split/\t/, $refdb{$tmp_chr}{$tmp_index})[4];
				$cnv_gap{$tmp_chr}{$start_index} += (split/\t/, $refdb{$tmp_chr}{$tmp_index})[5];
				push @tmp_median_kmer, $genome_copy{$tmp_chr}{$tmp_index};
			}
		}
		$cnv_real_end{$tmp_chr}{$start_index} = (split/\t/, $refdb{$tmp_chr}{$cnv{$tmp_chr}{$start_index}})[2] + $opts{k1}-1;
		@tmp_median_kmer = sort{$a <=> $b} @tmp_median_kmer;
		$kmer_median{$tmp_chr}{$start_index} = $tmp_median_kmer[int(0.5*@tmp_median_kmer)]; #the median of kmer
	}
}
###########################compute corrected coverage per sample per cnv region based on OUT1 file,using $median50###################
open (OUT, ">$ARGV[3]") or die "Merge CNVR output file name required!\n";
my %tmp_rd_sample; #store read depth in a CNVR for each sample
my $num;
print OUT "chr\tstart\tend\tnumber\tgap\trepeat\tgc\tkmer";
open (CORR, "$ARGV[1]") or die "corr list is required!\n";
my @tmp_file = <CORR>;
my $sample = @tmp_file;
for $num(0..$#tmp_file) {
	my $tmp_file_name = shift @tmp_file;
	my $suffix = (split/\//, $tmp_file_name)[-1];
	my @tmp = split/_mean_\d+\.?\d+?_SD/,$suffix;
	print OUT "\t$tmp[0]";
}
print OUT "\taverage\tsd\n";
foreach my $tmp_chr (sort{$a cmp $b} keys %cnv){
	foreach my $start_index (sort {$a <=> $b} keys %{$cnv{$tmp_chr}}){
		my $tmp_effective_wind_num = 0;
		foreach my $tmp_index ($start_index..$cnv{$tmp_chr}{$start_index}){
			if (exists $gain_effective{$tmp_chr}{$tmp_index} or exists $loss_effective{$tmp_chr}{$tmp_index}){
				$tmp_effective_wind_num++;
				for $num (0..($sample-1)){
					push @{$tmp_rd_sample{$num}}, (split/\t/,$hash_num{$tmp_chr}{$tmp_index})[$num];
				}
			}
		}
		my $repeat_ration = $tmp_effective_wind_num > 0 ? sprintf "%.2f",$cnv_repeat{$tmp_chr}{$start_index}/($tmp_effective_wind_num*$opts{k1}) : 0;
		my $gc_ration = $tmp_effective_wind_num > 0 ? sprintf "%.2f",$cnv_GC{$tmp_chr}{$start_index}/($tmp_effective_wind_num*$opts{k1}) : 0;
		my $gap_ration = $tmp_effective_wind_num > 0 ? sprintf "%.2f",$cnv_gap{$tmp_chr}{$start_index}/($tmp_effective_wind_num) : 0;
		print OUT "$tmp_chr\t$cnv_real_start{$tmp_chr}{$start_index}\t$cnv_real_end{$tmp_chr}{$start_index}\t$tmp_effective_wind_num\t$gap_ration\t$repeat_ration\t$gc_ration\t$kmer_median{$tmp_chr}{$start_index}";
		my @cp_median_CNVR;
		for $num (0..($sample-1)){
			my @sample_CNVR_cp = sort {$a <=> $b} @{$tmp_rd_sample{$num}};
			my $cp_median = $sample_CNVR_cp[int(0.5*@sample_CNVR_cp)];
			push @cp_median_CNVR, $cp_median;
			print OUT "\t$cp_median";
		}
		my $average_cp = &get_average(@cp_median_CNVR);
		my $sd_cp = &get_sd(@cp_median_CNVR);
		print OUT "\t$average_cp\t$sd_cp\n";
		%tmp_rd_sample = ();
	}
}
close OUT;
############################################################subroutine#################################################################
sub get_average {
	my @array = @_;
	my $sum;
	for my $value (@array) {
		$sum += $value;
	}
	my $average = @array > 0 ? $sum/@array : 0;
	return sprintf "%.2f",$average;
}
sub get_sd {
	my @array = @_;
	my $sum;
	for my $value (@array) {
		$sum += $value;
	}
	my $average = @array > 0 ? $sum/@array : 0;
	$sum = 0;
	for my $value (@array) {
		$sum += ($value-$average)**2;
	}
	my $sd = @array > 1 ? sqrt($sum/(@array-1)) : 0;
	return sprintf "%.2f",$sd;
}
