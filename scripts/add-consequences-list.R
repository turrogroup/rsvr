#! /usr/bin/Rscript

suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1)
	stop("Usage: Rscript static-data.R <target .db to add static data to>")

add_csqs <- function(db) {
	csq_list <- read.table(
		col.names=c("BIT", "NAME"),
		colClasses=c("integer", "character"),
		sep=",",
		text=c(
			'0,transcript_variant',
			'1,intron_variant',
			'2,exon_variant',
			'3,five_prime_UTR_variant',
			'4,three_prime_UTR_variant',
			'5,coding_sequence_variant',
			'6,synonymous_variant',
			'7,splice_region_variant',
			'8,extended_intronic_splice_region_variant',
			'9,exonic_splice_region_variant',
			'10,stop_retained_variant',
			'11,missense_variant',
			'12,conservative_missense_variant',
			'13,non_conservative_missense_variant',
			'14,inframe_deletion',
			'15,conservative_inframe_deletion',
			'16,disruptive_inframe_deletion',
			'17,inframe_insertion',
			'18,conservative_inframe_insertion',
			'19,disruptive_inframe_insertion',
			'20,start_lost',
			'21,stop_lost',
			'22,frameshift_variant',
			'23,stop_gained',
			'24,splice_donor_variant',
			'25,splice_acceptor_variant'
		)
	)
	dbExecute(db, "CREATE TABLE `CONSEQUENCES_LIST` ( `BIT` INTEGER NOT NULL PRIMARY KEY, `NAME` VARCHAR(255) NOT NULL)")
	dbWriteTable(db, name="CONSEQUENCES_LIST", value=csq_list, append=TRUE, overwrite=FALSE)
}

db <- dbConnect(SQLite(), args[1])
has_tabs <- dbListTables(db)
if (!("CONSEQUENCES_LIST" %in% has_tabs)) add_csqs(db)
dbDisconnect(db)

