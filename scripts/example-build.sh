#! /usr/bin/bash 

if [ $# -eq 0 ]; then {
	echo -e "Please supply a config file specifying:\nREF=<path to reference fasta file>" 
}; fi

CONF=`mktemp`
cat $1 | perl -le 'while (<STDIN>) { chomp; print "export $_" }' > ${CONF}

source ${CONF}

if [ -z ${REF+x} ]; then { 
	echo "REF must be set to location of reference genome fasta file"
	exit 1
}; fi
if [ -z ${VCF+x} ]; then { 
	echo "VCF must be set to location of merged VCF file"
	exit 1
}; fi
if [ -z ${ACLIM+x} ]; then { 
	echo "ACLIM must be set to maximum internal allele frequency (typically 0.2% of total allele number)"
	exit 1
}; fi
if [ -z ${GNOMAD+x} ]; then { 
	echo "GNOMAD must be set to location of the gnomAD VCF file"
	exit 1
}; fi
if [ -z ${BUILD+x} ]; then { 
	echo "BUILD must be set to the genome build version, either v37 or v38"
	exit 1
}; fi
if [ -z ${TARGET+x} ]; then { 
	echo "TARGET must be set to the path for the target SQLite database to be built"
	exit 1
}; fi
if [ -z ${CADD+x} ]; then { 
	echo "CADD must be set to the location of the genome wide CADD score file for SNVs (see https://cadd.gs.washington.edu/)"
	exit 1
}; fi
if [ -z ${GTF+x} ]; then { 
	echo "GTF must be set to the location of a GTF file containing the transcript features for which consequences should be predicted for each variant"
	exit 1
}; fi

INC="ALT!~\"N\"&AC<${ACLIM}&GT=\"alt\""

if test -f "${TARGET}"; then {
	echo "TARGET file '${TARGET}' already exists!"
	exit 1
}; fi

if ! command -v bcftools &> /dev/null
then
    echo "bcftools could not be found"
    exit
fi
if ! command -v sqlite3 &> /dev/null
then
    echo "sqlite3 could not be found"
    exit
fi

build-GENOTYPE.sh ${ACLIM} ${VCF} ${TARGET}
build-VARIANT.sh ${REF} ${CADD} ${GNOMAD} ${BUILD} ${TARGET}
build-CONSEQUENCE.sh ${REF} ${GTF} ${TARGET}
