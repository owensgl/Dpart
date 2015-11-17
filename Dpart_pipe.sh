#!/bin/bash

Name=$1
binpath=/home/owens/bin/Dpart
#Calculate the numerator and denominator for the D-statistic for each loci 
perl $binpath/Dpart.pl /home/owens/exil/2015/Exil.GATK.2015.tab /home/owens/exil/2015/Exil.GATK.2015.poplist.txt ${Name}_groups.txt > ${Name}_out.txt
#Sum up numerator and denominator for each scaffold and merge together to form blocks
perl $binpath/Dpart_out_blocker.pl ${Name}_out.txt > ${Name}_out_block.txt
#Jackknife bootstrap through the blocks to get standard error
Rscript $binpath/Dpart_jackknife.R ${Name}_out_block.txt ${Name}_out_final.txt

