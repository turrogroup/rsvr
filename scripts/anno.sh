#! /bin/bash

# annotate input RSVR IDs with gnomAD PMAF scores and CADD phred scores for given chromosomal regions.

REF=$1
CADD=$2
GNOMAD=$3
BUILD=$4
REGIONS=$5

CADD_ANNO=`mktemp -u`
GNOM_ANNO=`mktemp -u`
mkfifo ${CADD_ANNO}
tabix ${CADD} ${REGIONS} | rsvr enc -p -s -x | cut -d$'\t' -f 1,7 > ${CADD_ANNO} &

mkfifo ${GNOM_ANNO}
bcftools view -r ${REGIONS} -H ${GNOMAD} | sed 's/^chr//' | rsvr enc -s -p -i 1,2,4,5 |	rsvrgnomfreq.pl ${BUILD} | rsvr sort | rsvr pmaf > ${GNOM_ANNO} &

cat - | rsvr dec -p -f ${REF} | rsvr ann -1 -m null -f ${CADD_ANNO} | rsvr ann -1 -m 3 -f ${GNOM_ANNO} | perl -lne '@F=split /\t/;$F[3]="\x27$F[3]\x27";print join(",",@F[(0,3,5,6)])'

rm ${GNOM_ANNO}
rm ${CADD_ANNO}

