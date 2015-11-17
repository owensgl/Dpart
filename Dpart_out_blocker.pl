#!/usr/bin/perl

use warnings;
use strict;
#This takes in the ABBA BABA output, sums all applicable values within a scaffold


my $in = $ARGV[0]; #The output of ABBA_BABA.pl with 3


my $D1num = 0;
my $D1den = 0;
my $D2num = 0;
my $D2den = 0;
my $D12num = 0;
my $D12den = 0;
my $current_chrom;
my $window_size = 10000000;
my $current_end = $window_size;
my %colname;
open IN, $in;
while (<IN>){
	chomp;
	my @a = split (/\t/, $_);
	if ($. == 1){
		foreach my $i (0..$#a){
			$colname{$i}= $a[$i];
		}
		print "chr\tend\tD1num\tD1den\tD2num\tD2den\tD12num\tD12den";
	}else{
		unless($#a eq 0){
			my $chrom = $a[0];
			my $pos = $a[1];
			unless ($current_chrom){
				$current_chrom = $chrom;
				until($pos < $current_end){
					$current_end += $window_size;
				}
			}
			if (($current_chrom eq $chrom) and ($pos < $current_end)){
				foreach my $i (0..$#a){
					if ($colname{$i} eq "D1num"){
						$D1num += $a[$i];
					}elsif ($colname{$i} eq "D1den"){
						$D1den += $a[$i];
					}elsif ($colname{$i} eq "D2num"){
	                               	        $D2num += $a[$i];
					}elsif ($colname{$i} eq "D2den"){
	                                        $D2den += $a[$i];
					}elsif ($colname{$i} eq "D12num"){
	                                        $D12num += $a[$i];
					}elsif ($colname{$i} eq "D12den"){
						$D12den += $a[$i];
					}
				}
			}else{ #End of window
				
				print "\n$current_chrom\t$current_end\t$D1num\t$D1den\t$D2num\t$D2den\t$D12num\t$D12den";
				if ($current_chrom ne $chrom){
					$current_chrom = $chrom;
					$current_end = $window_size;
	                                until($pos < $current_end){
        	                                $current_end += $window_size;
                	                }
				}else{
					until($pos < $current_end){
                                                $current_end += $window_size;
                                        }
				}
                                foreach my $i (0..$#a){
                                        if ($colname{$i} eq "D1num"){
                                                $D1num = $a[$i];
                                        }elsif ($colname{$i} eq "D1den"){
                                                $D1den = $a[$i];
                                        }elsif ($colname{$i} eq "D2num"){
                                                $D2num = $a[$i];
                                        }elsif ($colname{$i} eq "D2den"){
                                                $D2den = $a[$i];
                                        }elsif ($colname{$i} eq "D12num"){
                                                $D12num = $a[$i];
                                        }elsif ($colname{$i} eq "D12den"){
                                                $D12den = $a[$i];
                                       }
                                }
                        }
		}
	}
}
			
	
