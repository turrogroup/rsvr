#! /bin/bash

if [ $# -eq 0 ] 
then
	echo "Usage: build-GENOTYPE.sh <AC threshold> <path to merged VCF> <target DB file>";
	exit;
fi

ACLIM=$1
VCF=$2
DB=$3

INC="ALT!~\"N\"&AC<${ACLIM}&GT=\"alt\"&FILTER=\"PASS\""

sqlite3 ${DB} "CREATE TABLE \`GENOTYPE\`(\`RSVR_ID\` BIGINT NOT NULL, \`SAMPLE_ID\` INTEGER NOT NULL, \`GT\` TINYINT NOT NULL)"

bcftools view -i ${INC} ${VCF} |
perl -lne 'if (/^#CHROM/) { @F=split /\t/; print join("\t", (@F[0..8], (1..((scalar @F)-9)))) } else { print }' |
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT]\n' |
sed 's/^chr//' |
rsvr enc -s -p -t 100 |
perl -lne 'while (<STDIN>) { chomp; @F=split /\t/; @g=map {/^(\d+):(\d|\.)(?:(\/|\|)(\d|\.))?/;$1."\t".($2+$4)} @F[5..$#F]; $vid=$F[0]; for $i (@g) { print $vid."\t".$i; } }' |
rsvr gtnorm |
inserts.pl BEGIN GENOTYPE 1 |
sqlite3 ${DB}

sqlite3 ${DB} "CREATE INDEX \`GENOTYPE-RSVR_ID\` ON \`GENOTYPE\`(\`RSVR_ID\`)"
sqlite3 ${DB} "CREATE INDEX \`GENOTYPE-SAMPLE_ID\` ON \`GENOTYPE\`(\`SAMPLE_ID\`)"

