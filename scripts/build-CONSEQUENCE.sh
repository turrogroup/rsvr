#! /bin/bash

if [ $# -eq 0 ] 
then
	echo "Usage: build-CONSEQUENCE.sh <path to reference fasta> <path to GTF file> <target DB file containing VARIANT table>";
	exit;
fi

REF=$1
GTF=$2
DB=$3

FEATURES_DB=`mktemp`
gtf2featuresDB.R ${GTF} ${REF} ${FEATURES_DB}

FEATURES_TXT=`mktemp`
sqlite3 -separator $'\t' ${FEATURES_DB} "SELECT * FROM COMPLETED_FEATURE" > ${FEATURES_TXT}

IDS=`mktemp`
sqlite3 ${DB} "SELECT RSVR_ID FROM VARIANT ORDER BY RSVR_ID" > ${IDS}

sqlite3 ${DB} "CREATE TABLE CONSEQUENCE(RSVR_ID BIGINT NOT NULL, GENE_ID INTEGER NOT NULL, TX_ID INTEGER NOT NULL, CDNA INTEGER NULL, CDS INTEGER NULL, WORST_CSQ_ID BIGINT NOT NULL, ALL_CSQ_ID BIGINT NOT NULL, CANONICAL_CSQ_ID BIGINT NOT NULL)"

cat ${IDS} | rsvr seqfx -f ${REF} -g ${FEATURES_TXT} | perl -lne '@F=split /\t/; print join ",", map { $_ eq "" ? "null" : $_ } @F' | inserts.pl BEGIN CONSEQUENCE 1 | sqlite3 ${DB}

sqlite3 ${DB} "CREATE INDEX \`CONSEQUENCE-RSVR_ID\` ON \`CONSEQUENCE\`(\`RSVR_ID\`);"
sqlite3 ${DB} "CREATE INDEX \`CONSEQUENCE-GENE_ID\` ON \`CONSEQUENCE\`(\`GENE_ID\`);"
sqlite3 ${DB} "CREATE INDEX \`CONSEQUENCE-TX_ID\` ON \`CONSEQUENCE\`(\`TX_ID\`);"
sqlite3 ${FEATURES_DB} .dump | sqlite3 ${DB}
add-consequences-list.R ${DB}
sqlite3 ${DB} "DROP VIEW COMPLETED_FEATURE"

