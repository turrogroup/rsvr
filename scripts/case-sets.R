#! /usr/bin/R

library(DBI)
library(RSQLite)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
	stop("Usage: case-sets.R <tab delimited text file sample name-case set identifier pairs> <DB file to create table in>\n\tNote that DB file must contain the SAMPLE table with a NAME column, as this is used to assign database IDs from the SAMPLE table to the case sets. Case set identifiers given in the text file can be strings which identify the diagnosis (e.g. 'Thrombocytopenia'). The text file should be formatted as lines containing <sample name><tab><case set identifier>.\n")
}

case_sets_file <- args[1]
db_file <- args[2]
if (!all(file.exists(args))) {
	stop("file not found!")
}

case_sets_table <- read.table(comment.char="", header=FALSE, col.names=c("sample", "case_set"), colClasses=rep("character", 2), stringsAsFactors=FALSE, file=case_sets_file)
db <- dbConnect(SQLite(), db_file)

if (!("SAMPLE" %in% dbListTables(db)) | !("NAME" %in% dbListFields(db, "SAMPLE"))) {
	stop("DB must contain a 'SAMPLE' table with a 'NAME' column")
}

sample_names <- dbGetQuery(db, "SELECT ID, NAME FROM SAMPLE")
dbWriteTable(db, name="CASE_SET", value=data.frame(SAMPLE_ID=sample_names$ID[match(case_sets_table$sample, sample_names$NAME)], SET_ID=case_sets_table$case_set), append=FALSE, overwrite=TRUE)
dbDisconnect(db)
