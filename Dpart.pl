#!/usr/bin/perl


#Sept 2014. Updated to make denominator and numerator sum for each loci and then average
use warnings;
use strict;
use lib '/home/owens/bin'; #Location of countbadcolumns.pl
my %t;
$t{"N"} = "NN";
$t{"A"} = "AA";
$t{"T"} = "TT";
$t{"G"} = "GG";
$t{"C"} = "CC";
$t{"W"} = "TA";
$t{"R"} = "AG";
$t{"M"} = "AC";
$t{"S"} = "CG";
$t{"K"} = "TG";
$t{"Y"} = "CT";

my $in = $ARGV[0]; #In SNP table in iupac or not. Will detect the number of columns before data
my $pop = $ARGV[1]; #Population file, (samplename\tpopname\n)
my $poporder = $ARGV[2]; #The file specify which populations represent which samples in the ABBA BABA scheme (popname\t[1-4]\n)

require "countbadcolumns.pl";
my ($iupac_coding, $badcolumns) = count_bad_columns($in);
$. = 0;

my %pop;
my %poplist;
my %group;
my %samplegroup;
#Load population information
open POP, $pop;
while (<POP>){
	chomp;
	my @a = split (/\t/,$_);	
	$pop{$a[0]}=$a[1];
	$poplist{$a[1]}++;
}
close POP;
#Load group information
open GROUP, $poporder;
while (<GROUP>){
	chomp;
	my @a = split (/\t/,$_);
	$group{$a[0]} = $a[1];
}
close GROUP;

open IN, $in;
while (<IN>){
	chomp;
	next if /^\s*$/;
	my @a = split (/\t/,$_);
  	if ($. == 1){
  		foreach my $i ($badcolumns..$#a){ #Get sample names for each column
        		if ($pop{$a[$i]}){
				if ($group{$pop{$a[$i]}}){
	        			$samplegroup{$i} = $group{$pop{$a[$i]}};
				}
        		}
        	}
		print "chr\tpos\tD1num\tD1den\tD2num\tD2den\tD12num\tD12den";
	}else{
   		my $pos = "$a[0]\t$a[1]";
    		my %BC;
		my %total_alleles;
		#Load in data for each individual and sum up allele frequencies by group 
    		foreach my $i ($badcolumns..$#a){
			if ($samplegroup{$i}){
				$BC{"total"}{"total"}++;
				if ($iupac_coding eq "TRUE"){
					$a[$i] = $t{$a[$i]};
				}
				unless (($a[$i] eq "NN")or($a[$i] eq "XX")){
					my @bases = split(//, $a[$i]);
					$total_alleles{$bases[0]}++;
					$total_alleles{$bases[1]}++;
					
					$BC{"total"}{$bases[0]}++;
		        		$BC{"total"}{$bases[1]}++;
					$BC{$samplegroup{$i}}{$bases[0]}++;
		 			$BC{$samplegroup{$i}}{$bases[1]}++;	

					$BC{"total"}{"Calls"}++;
					$BC{$samplegroup{$i}}{"Calls"}++;
				}
			}
		}

		my %B;
		my %A;
		my $missing;
		
		if (keys %total_alleles == 2){ #Only look at biallelic sites
			my @bases = sort { $total_alleles{$b} <=> $total_alleles{$a} } keys %total_alleles ;
			my @outgroupbases = sort { $BC{5}{$b} <=> $BC{5}{$a} } keys %{$BC{5}}; #Sort bases in outgroup
			if ($outgroupbases[0]){
				my $A = $outgroupbases[0]; #Ancestral allele (A)
				my $B;
				if ($bases[0] eq $outgroupbases[0]){
					$B = $bases[1]; #Derived allele (B)
				}else{
					$B = $bases[0]; #Derived allele (B)
				}
				foreach my $i(1..5){
					if ($BC{$i}{"Calls"}){
						if ($BC{$i}{$B}){
							$B{$i} = $BC{$i}{$B}/($BC{$i}{"Calls"} * 2);
							$A{$i} = (1 - $B{$i});
						}else{
							$B{$i} = 0;
							$A{$i} = 1;
						}
					}else{
						$B{$i} = "NA";
						$A{$i} = "NA";
						$missing++;
					}
				}
                                if ($B{5} ne "0"){
                                	$missing++;
                                }
				unless ($missing){
					my $ABBBA = ($A{"1"} * $B{"2"} * $B{"3"} * $B{"4"} * $A{"5"});
					my $BABBA = ($B{"1"} * $A{"2"} * $B{"3"} * $B{"4"} * $A{"5"});
					my $ABBAA = ($A{"1"} * $B{"2"} * $B{"3"} * $A{"4"} * $A{"5"});
					my $BABAA = ($B{"1"} * $A{"2"} * $B{"3"} * $A{"4"} * $A{"5"});
					my $ABABA = ($A{"1"} * $B{"2"} * $A{"3"} * $B{"4"} * $A{"5"});
					my $BAABA = ($B{"1"} * $A{"2"} * $A{"3"} * $B{"4"} * $A{"5"});
					my $D1Num = ($ABBAA - $BABAA);
					my $D1Denom = ($ABBAA + $BABAA); 
					my $D2Num = ($ABABA - $BAABA);
					my $D2Denom = ($ABABA + $BAABA);
					my $D12Num = ($ABBBA - $BABBA);
					my $D12Denom = ($ABBBA + $BABBA);
					print "\n$pos";
					print "\t$D1Num\t$D1Denom\t$D2Num\t$D2Denom\t$D12Num\t$D12Denom";
					
				}
			}
		}
	}
}
close IN; 

