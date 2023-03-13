#! /bin/bash

if [ $# -eq 0 ] 
then
	echo "Usage: build-VARIANT.sh <path to reference fasta> <path to CADD tsv> <path to gnomAD VCF> <genome build> <DB containing GENOTYPE table>";
	exit;
fi

REF=$1
CADD=$2
GNOMAD=$3
GENOME_BUILD=$4
DB=$5

INSERTS=`mktemp -u`
CADD_ANNO=`mktemp -u`
GNOM_ANNO=`mktemp -u`

REGS=`mktemp` 
perl -le 'for (1..22,"X","Y","MT") { print }' > ${REGS}

mkfifo ${CADD_ANNO}
tabix -R ${REGS} ${CADD} | rsvr enc -p -s -x | cut -d$'\t' -f 1,7 > ${CADD_ANNO} &

mkfifo ${GNOM_ANNO}
bcftools view -H ${GNOMAD} | sed 's/^chr//' | rsvr enc -s -p -i 1,2,4,5 | rsvrgnomfreq.pl ${GENOME_BUILD} | rsvr sort | rsvr pmaf > ${GNOM_ANNO} &

sqlite3 ${DB} "SELECT DISTINCT RSVR_ID FROM GENOTYPE ORDER BY RSVR_ID" | rsvr dec -p -f ${REF} | rsvr ann -1 -m null -f ${CADD_ANNO} | rsvr ann -1 -m 3 -f ${GNOM_ANNO} | perl -lne '@F=split /\t/;$F[3]="\x27$F[3]\x27";print join(",",@F[(0,3,5,6)])' | inserts.pl BEGIN VARIANT 1 | gzip -c > ${INSERTS}

sqlite3 ${DB} "CREATE TABLE \`VARIANT\`( \`RSVR_ID\` BIGINT NOT NULL PRIMARY KEY, \`REF\` VARCHAR(63) NOT NULL, \`CADD_PHRED\` FLOAT NULL, \`PMAF\` TINYINT NOT NULL)"

zcat ${INSERTS} | sqlite3 ${DB}

sqlite3 ${DB} "CREATE INDEX \`VARIANT-RSVR_ID\` ON \`VARIANT\`(\`RSVR_ID\`)"

