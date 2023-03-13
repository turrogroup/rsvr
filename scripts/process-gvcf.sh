#! /bin/bash

# process single sample gVCF to get (CHROM POS REF ALT genotype passing) quintuplets (with passing determined by a perl regex pattern to match on the FILTER field and quantity/threshold pair (e.g. GQ >= 30).

gvcf=$1
ref=$2
pass=$3
quantity=$4
threshold=$5
bcftools norm -c w -m - -f ${ref} ${gvcf} 2> /dev/null | 
bcftools query -i 'GT="alt"' -f "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%FILTER]\t[%${quantity}]\n" | 
sed 's/^chr//' |
perl -le 'while (<STDIN>) { chomp;@F=split /\t/; $c=$F[0]; $p=int($F[5]=~/$ARGV[0]/ and ($F[6] >= $ARGV[1]));print join("\t",(@F[0..4],$p)) if (($c=~/^\d+$/ and $c > 0 and $c <= 25) or ($c eq "X" or $c eq "Y" or $c eq "M" or $c eq "MT")) and $F[3]=~/^[ACGT]+$/ }' ${pass} ${threshold}
