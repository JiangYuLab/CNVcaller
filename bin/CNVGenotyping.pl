#!/usr/bin/perl
use strict;
use FindBin '$Bin';
use Data::Dumper;
die "Usage : perl $0 <CNVR> <sample_number>\n" unless @ARGV == 2;
open (IN1, "$Bin/database") or die "database missing in database directory!\n";
open (IN2, "$ARGV[0]") or die "Input CNVR file!\n";
my $sample_num = $ARGV[1];
my $i;
my $n_01;my $n_02;my $n_11;my $n_12;my $n_21;my $n_22;my $n_31;my $n_32;
my %database;my %database2;my %database3;my %database4;my %database5;my %database6;my %database7;my %database8;
my $m1; my $m2; my $a; my $b;
my $x01;my $x02;my $x11;my $x12;my $x21;my $x22;my $x31;my $x32;
my $pl;my $average;
my $select; my $rd;
my $genotyping_aa;my $genotyping_Aa;my $genotyping_AA; my $genotyping_AB; my $genotyping_BB;my $genotyping_BC;my $genotyping_M;
my @biaotou;my @zhongjian;
my $n1;my $n2;my $n3;

while(<IN1>){
	next if /^#/;
	chomp;
	my @data = split/\s+/;
	$database{"$data[0]-$data[1]"} = $data[2];
	$database2{"$data[0]-$data[1]"} = $data[3];
	$database3{"$data[0]-$data[1]"} = $data[2]+1;
	$database4{"$data[0]-$data[1]"} = $data[3]+1;
	$database5{"$data[0]-$data[1]"} = $data[2]+1.5;
	$database6{"$data[0]-$data[1]"} = $data[3]+1.5;
	$database7{"$data[0]-$data[1]"} = $data[2]+2;
	$database8{"$data[0]-$data[1]"} = $data[3]+2;
	
}
close IN1;

my $header = <IN2>;
chomp $header;
my @biaotou = split/\s+/, $header;
my @zhongjian = splice(@biaotou,8, $sample_num);
print $header,"\t",join("\t",@zhongjian),"\t","aa","\t","Aa","\t","AA","\t","AB","\t","BB","\t","BC","\t","M","\t","average1\taverage2\taverage3\tsd1\tsd2\tsd3","\n";
while(<IN2>){
	chomp;
	print $_;
	my @tmpkaishi = split/\s+/;
	my @average = splice(@tmpkaishi, $sample_num+8, 1);
	my @tmp = splice(@tmpkaishi, 8, $sample_num);
	if( $average[0]>=0 && $average[0]< 1.00 ){
		$pl = 0;
	}
	elsif( $average[0]>=1.00 && $average[0]< 2.00){
		$pl = 1.00;
	}
	elsif( $average[0]>=2.00 && $average[0]< 2.5 ){
		$pl = 1.50;
	}
	elsif($average[0]>= 2.5){
		$pl = 2;
	}
	$n_01=$n_02=0;
	for $i(0..@tmp-1){
		if(@tmp[$i] <= 0.25+$pl){
			$n_01++;
		} 
		elsif(@tmp[$i] > 0.25+$pl && @tmp[$i ]<= 0.75+$pl){
			$n_02++;
		}	
	}
	$x01= sprintf"%.2f",$n_01/@tmp;
	$x02= sprintf"%.2f",$n_02/@tmp;
	if ($x01+$x02>1){$x02 = sprintf"%.2f",(1-$x01)};
	if( $average[0]>=0 && $average[0]< 1.00 ){
		$a=$database{"$x01-$x02"};                        
		$b=$database2{"$x01-$x02"};
	}
	elsif( $average[0]>=1.00 && $average[0]< 2.00){
		$a=$database3{"$x01-$x02"};                        
		$b=$database4{"$x01-$x02"};   
	}
	elsif( $average[0]>=2.00 && $average[0]< 2.5 ){
		$a=$database5{"$x01-$x02"};                        
		$b=$database6{"$x01-$x02"};
	}
	elsif($average[0]>= 2.5){
		$a=$database7{"$x01-$x02"};                        
		$b=$database8{"$x01-$x02"} ;
	}
	$n_11=$n_12=0;
	for $i(0..@tmp-1){
		if(@tmp[$i]<=$a){
			$n_11++;
		}
		elsif(@tmp[$i]>$a && @tmp[$i]<=$b){
			$n_12++;
		}
	}
	$x11= sprintf"%.2f",$n_11/@tmp;
	$x12= sprintf"%.2f",$n_12/@tmp;
	if ($x11+$x12>1){$x02 = sprintf"%.2f",(1-$x11)};
	if( $average[0]>=0 && $average[0]< 1.00 ){
		$a=$database{"$x11-$x12"};
		$b=$database2{"$x11-$x12"};
	}
	elsif( $average[0]>=1.00 && $average[0]< 2.00){
		$a=$database3{"$x11-$x12"};
		$b=$database4{"$x11-$x12"};   
	}
	elsif( $average[0]>=2.00 && $average[0]< 2.5 ){
		$a=$database5{"$x11-$x12"};
		$b=$database6{"$x11-$x12"};
	}
	elsif($average[0]>= 2.5){
		$a=$database7{"$x11-$x12"};
		$b=$database8{"$x11-$x12"};
	}
	$n_21=$n_22=0;
	for $i(0..@tmp-1){
		if(@tmp[$i]<=$a){
			$n_21+=1;
		}
		elsif(@tmp[$i]>$a && @tmp[$i]<=$b){
			$n_22++;
		}
	}
	$x21= sprintf"%.2f",$n_21/@tmp;
	$x22= sprintf"%.2f",$n_22/@tmp;
	if ($x21+$x22>1){$x22 = sprintf"%.2f",(1-$x21)};
	if( $average[0]>=0 && $average[0]< 1.00 ){
		$a=$database{"$x21-$x22"};
		$b=$database2{"$x21-$x22"};
	}
	elsif( $average[0]>=1.00 && $average[0]< 2.00){
		$a=$database3{"$x21-$x22"};                        
		$b=$database4{"$x21-$x22"};
	}
	elsif( $average[0]>=2.00 && $average[0]< 2.5 ){
		$a=$database5{"$x21-$x22"};                        
		$b=$database6{"$x21-$x22"};
	}
	elsif( $average[0]>= 2.5){
		$a=$database7{"$x21-$x22"};
		$b=$database8{"$x21-$x22"};
	}
	$m1=$x21;
	$m2=$x22;
	$a=$a;
	$b=$b;
	print  "\t";
	
#################################################################################################################genotyping
	my @select = @tmp;
	my (@mode1, @mode2, @mode3);
	$genotyping_aa=$genotyping_Aa=$genotyping_AA=$genotyping_AB=$genotyping_BB=$genotyping_BC=$genotyping_M=0;
	foreach my $rd (@select){
		if ($a < 0.45 && $a > 0.1  and $b > 0.49 && $b < 0.96 ){
			if ($rd <= $a){
				print "dd\t";
				$genotyping_aa++;
				push @mode1, $rd;
			}
			elsif( $rd > $a and $rd <= $b){
				print "Ad\t";
				$genotyping_Aa++;
				push @mode2, $rd;
			}
			elsif($rd > $b){
				print "AA\t";
				$genotyping_AA++;
				push @mode3, $rd;
			}
		}
		elsif ($a < 1.45 && $a > 1.1 and $b > 1.49 && $b < 1.96){            
			if($rd <= $a){
				print "AA\t";
				$genotyping_AA++;
				push @mode1, $rd;
			}				   
			elsif($rd >$a and $rd <= $b){
				print "AB\t";
				$genotyping_AB++;
				push @mode2, $rd;
			}
			elsif($rd > $b){
				print "BB\t";
				$genotyping_BB++;
				push @mode3, $rd;
			}
		}
		elsif ($a < 1.95 && $a > 1.6 and $b > 1.99 && $b <2.46 ){
			if($rd <= $a){
				print "AB\t";
				$genotyping_AB++;
				push @mode1, $rd;
			}
			elsif($rd >$a and $rd <= $b){
				print "BB\t";
				$genotyping_BB++;
				push @mode2, $rd;
			}
			elsif($rd > $b){
				print "BC\t";
				$genotyping_BC++;
				push @mode3, $rd;
			}
		}
		elsif ( $a < 2.45 && $a > 2.1 and $b < 2.99){
			if($rd <= $a){
				print "BB\t";
				$genotyping_BB++;
				push @mode1, $rd;
			}
			elsif($rd > $a and $rd <= $b){
				print "BC\t";
				$genotyping_BC++;
				push @mode2, $rd;
			}elsif($rd >$b ){
				print "M\t";
				$genotyping_M++;
				push @mode3, $rd;
			}
		}
	}
	my $averagemode1 = &get_average(@mode1);my $sdmode1 = &get_sd(@mode1);
	my $averagemode2 = &get_average(@mode2);my $sdmode2 = &get_sd(@mode2); 
	my $averagemode3 = &get_average(@mode3);my $sdmode3 = &get_sd(@mode3);
	print "$genotyping_aa\t$genotyping_Aa\t$genotyping_AA\t$genotyping_AB\t$genotyping_BB\t$genotyping_BC\t$genotyping_M\t$averagemode1\t$averagemode2\t$averagemode3\t$sdmode1\t$sdmode2\t$sdmode3\n"; 
}
close IN2;
sub get_average {
	my @array = @_;
	my $sum;
	for my $value (@array) {
		$sum += $value;
	}
	my $average = @array > 0 ? sprintf "%.2f",$sum/@array : "NA";
	return $average;
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
	my $sd = @array > 1 ? sprintf "%.2f",sqrt($sum/(@array-1)) : "NA";
	return $sd;
}
