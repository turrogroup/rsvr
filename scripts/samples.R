#! /usr/bin/R

library(DBI)
library(RSQLite)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
	stop("Usage: samples.R <tab delimited samples table> <.db file to create the table in>\n\tNote that the samples table should have the named columns in the first row, of which one should be NAME (typically the sample names given in a merged VCF), and samples must be in same order as they appear in columns in merged VCF (an ID column will be generated automatically). Given column names will be used in the database (for example, additional useful columns may include sex, family identifiers and ancestry).\n")
}

sample_file <- args[1]
db_file <- args[2]
if (!all(file.exists(sample_file))) {
	stop("file not found!")
}

sample_table <- read.table(comment.char="", header=TRUE, stringsAsFactors=FALSE, file=sample_file)
db <- dbConnect(SQLite(), db_file)
dbWriteTable(db, name="SAMPLE", value=cbind(ID=seq_len(nrow(sample_table)), sample_table), append=FALSE, overwrite=TRUE)
dbDisconnect(db)
