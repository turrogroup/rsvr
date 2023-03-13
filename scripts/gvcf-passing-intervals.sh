#! /usr/bin/bash

GVCF=$1
INC=$2
if [ -z ${GVCF+x} ]; then { 
	echo "Usage: gvcf-passing-intervals.sh <GVCF file path> <bcftools include '-i' parameter>"
	exit 1
}; fi
bcftools query -i "${INC}" -f "%CHROM\t%POS\t[%INFO/END]\n" ${GVCF} | 
perl -le '$re=-1;while (<STDIN>) { /^(chr)?(\w+)\t(\d+)\t(\d+|\.)/; if ($2 eq "X") { $c=23 } elsif ($2 eq "Y") { $c=24 } elsif ($2 eq "MT") { $c=25 } else { $c=$2 };$n=$c << 28;$s=$n + $3; if ($4 eq ".") { $e = $s + 1 } else { $e = $n + $4 + 1 } if ($s > $re) { if ($re > -1) { print "$rs\t$re" } $rs = $s } if ($e > $re) { $re = $e } } print "$rs\t$re";'
