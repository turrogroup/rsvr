#! /usr/bin/Rscript

refseq_types <- read.table(stringsAsFactors=FALSE, col.names=c("TYPE","RSVR_TYPE"), sep=",", text=c(
	"CDS,CDS",
	"exon,exon",
	"five_prime_utr,five_prime_utr",
	"gene,gene",
	"Selenocysteine,",
	"start_codon,",
	"stop_codon,",
	"three_prime_utr,three_prime_utr",
	"transcript,transcript",
	"intron,transcript",
	"miRNA,transcript"
))

ensembl_types <- read.table(stringsAsFactors=FALSE, col.names=c("TYPE","RSVR_TYPE"), sep=",", text=c(
	"aberrant_processed_transcript,transcript",
	"CDS,CDS",
	"C_gene_segment,gene",
	"chromosome,",
	"exon,exon",
	"five_prime_UTR,five_prime_utr",
	"gene,gene",
	"J_gene_segment,gene",
	"lincRNA,transcript",
	"lincRNA_gene,transcript",
	"miRNA,transcript",
	"miRNA_gene,transcript",
	"mRNA,transcript",
	"mt_gene,gene",
	"nc_primary_transcript,transcript",
	"NMD_transcript_variant,transcript",
	"processed_pseudogene,transcript",
	"processed_transcript,transcript",
	"pseudogene,gene",
	"pseudogenic_transcript,transcript",
	"RNA,transcript",
	"rRNA,transcript",
	"rRNA_gene,gene",
	"snoRNA,transcript",
	"snoRNA_gene,gene",
	"snRNA,transcript",
	"snRNA_gene,gene",
	"three_prime_UTR,three_prime_utr",
	"transcript,transcript",
	"VD_gene_segment,gene",
	"V_gene_segment,transcript"
))

ensembl_gtf_types <- read.table(stringsAsFactors=FALSE, col.names=c("TYPE","RSVR_TYPE"), sep=",", text=c(
	"CDS,CDS",
	"exon,exon",
	"five_prime_utr,exon",
	"gene,",
	"Selenocysteine,transcript",
	"start_codon,CDS",
	"stop_codon,CDS",
	"three_prime_utr,exon",
	"transcript,transcript"
))

RSVR_TRANSCRIPT_FEATURE_TYPES <- c("CDS", "exon", "five_prime_utr", "three_prime_utr", "transcript", "intron")
chroms <- c(1:22, "X", "Y", "MT")

id_pattern <- "[a-zA-Z0-9.:*$@!+_?|,-]+"

get_prop <- function(strings, property_name, fmt="%s=(%s)") vapply(
	FUN.VALUE=character(1), 
	FUN=function(x) if (length(x) == 2L) x[2L] else NA_character_, 
	X=regmatches(x=strings, m=regexec(pattern=sprintf(fmt, property_name, id_pattern), text=strings))
)

gtf_colnames <- c("chr","source","type","start","end","score","strand","phase","desc")

gtf2rtf_base <- function(x) {
	tab <- if (is.character(x)) {
		read.table(
			file=x,
			comment.char="#", 
			quote="", 
			header=FALSE, 
			col.names=gtf_colnames,
			stringsAsFactors=FALSE,
			sep="\t"
		) 
	} else if (is.data.frame(x)) { 
		stopifnot(identical(gtf_colnames, colnames(x)) )
		x
	} else {
		stop("'x' must be the file name of a GTF file or a data.frame containing GTF table")
	}

	id_columns <- lapply(
		setNames(nm=c("gene_id","transcript_id","transcript_biotype")), 
		function(property_name) get_prop(tab$desc, property_name, fmt="%s \"(%s)\"")
	)
	stopifnot(max(sapply(id_columns, function(x) max(lengths(strsplit(x, split=",")))))==1L)
	subset(do.call(
		what=function(...) cbind(
			stringsAsFactors=FALSE,
			tab[-match(c("source","type","score","desc"), names(tab))],
			otype=tab$type,
			type=with(ensembl_gtf_types, RSVR_TYPE[match(tab$type, TYPE)]),
			...
		),
		id_columns
	), chr %in% chroms)
}

get_biotype_tables <- function(desc_col, just_tx=FALSE, fmt="%s \"(%s)\"") {
	prop_types <- c("id","biotype")
	tabs <- lapply(
		setNames(nm=if (just_tx) "transcript" else c("gene", "transcript")),
		function(el) {
			subset(
				do.call(
					what=function(...) data.frame(
						stringsAsFactors=FALSE,
						...
					),
					lapply(
						setNames(nm=prop_types, sprintf("%s_%s", el, prop_types)),
						function(property_name) get_prop(desc_col, property_name, fmt=fmt)
					)
				),
				!is.na(id) & !duplicated(id)
			)
		}
	)
	if (just_tx) tabs[[1]] else tabs
}

get_first_cds <- function(type, start, end, strand, transcript_id) {
	stopifnot(is.factor(strand) && identical(levels(strand), c("+","-")))
	stopifnot(is.factor(transcript_id))
	unsplit(f=transcript_id, Map(
		split(
			data.frame(stringsAsFactors=FALSE, 
				is_cds=type=="CDS", 
				start, 
				end
			), 
			f=transcript_id
		),
		tapply(X=strand=="+", INDEX=transcript_id, FUN=function(x) { stopifnot(all(x)==x[1L]); x[1L] }),
		f=function(df, positive) {
			cds_f <- with(df, {
				cdsf <- factor(ifelse(is_cds,"CDS", "other"), levels=c("CDS", "other"))
				unsplit(
					f=cdsf,
					list(
						CDS=local({
							cds_ints <- cbind(start, end)[cdsf=="CDS",,drop=FALSE]
							if (max(c(0L, depth(cds_ints))) > 1L) {
								stop("depth greater than 1 of CDS regions!")
							}
							seq_len(nrow(cds_ints)) == which.min(cds_ints[,1] * (if (positive) 1L else -1L))
						}),
						other=rep(NA, table(cdsf)[2L])
					)
				)
			})
		}
	))
}

get_end_codon_starts <- function(is_cds, start, end, phase, positive) {
	cdsf <- factor(ifelse(is_cds,"CDS","other"), levels=c("CDS", "other"))
	l <- end-start
	ncodons <- ifelse(phase == 0L, 0L, 1L) + ifelse((l-phase) %% 3L == 0L, 0L, 1L) + (l-phase) %/% 3L
	stopifnot(all(ncodons[is_cds] > 0L))
	tb <- table(cdsf)
	x <- unsplit(
		f=cdsf,
		list(
			CDS=local({
				df <- if (positive) { 
					first <- start-ifelse(phase==0L, 0L, 3L-phase)
					data.frame(row.names=NULL, first=first, last=first+(ncodons-1L)*3L)
				} else { 
					e <- ifelse(phase == 0L, end-3L, end-phase)
					data.frame(first=e-3L*(ncodons-1L), last=e)
				}
				dfs <- with(df, data.frame(fs=pmax(start, first), fe=first+3L, ls=last, le=pmin(last+3L, end)))[is_cds,]
				row.names(dfs) <- paste0(rep("CDS", tb["CDS"]), seq_len(tb["CDS"]))
				dfs
			}),
			other=data.frame(row.names=paste0(rep("other",tb["other"]), seq_len(tb["other"])), first=rep(NA_integer_, tb["other"]), last=rep(NA_integer_, tb["other"]))
		)
	)
	row.names(x) <- NULL
	x
}

get_first_and_last_codon <- function(type, chr, start, end, phase, strand, transcript_id, ref, mc.cores=1L) {
	if (system("rsvr -v", ignore.stdout=TRUE) != 0) {
		stop("cannot find 'rsvr' program")
	}
	chrf <- factor(chr, levels=chroms)
	stopifnot(all(!is.na(chrf)))
	txf <- factor(transcript_id)
	codpos_parts <- (if (mc.cores==1L) Map else function(...) mcMap(mc.cores=mc.cores, ...))(
		split(
			data.frame(stringsAsFactors=FALSE, is_cds=type=="CDS", start=start, end=end, phase=phase),
			f=txf
		),
		tapply(X=strand=="+", INDEX=txf, FUN=function(x) { stopifnot(all(x)==x[1L]); x[1L] }),
		levels(txf),
		f=function(df, pos, txnm) { 
			o <- with(df, get_end_codon_starts(is_cds, start, end, phase, pos)) 
			row.names(o) <- paste0(rep(txnm, nrow(o)), seq_len(nrow(o)))
			o
		}
	)
	codpos_df <- do.call(what=data.frame, lapply(setNames(nm=c("fs", "fe", "ls", "le"), 1:4), function(x) unsplit(f=txf, value=lapply(codpos_parts, "[[", x))))
	stopifnot(identical(type!="CDS", is.na(codpos_df$fs)))
	codpos <- with(lapply(codpos_df, function(x) as.integer64(2)^28L * as.integer(chrf) + x), data.frame(from=c(fs, ls), to=c(fe, le)))
	codf <- factor(ifelse(is.na(codpos[[1]]), "other", "CDS"), levels=c("CDS","other"))
	matrix(
		dimnames=list(NULL, c("first", "last")),
		byrow=FALSE,
		ncol=2L,
		nrow=nrow(codpos_df),
		data=unsplit(
			with(split(codpos, f=codf), {
				ord <- order(as.integer(CDS[[1]] %/% (2L^28)), as.integer(CDS[[1]] %% (2L^28)))
				bases <- system2(
					command="rsvr",
					args=c("gref", "-f", ref),
					stdout=TRUE,
					input=with(CDS[ord,], sprintf("%s %s", from, to))
				)
				list(CDS=bases[order(ord)], other=rep(NA_character_, length(other)))
			}),
			f=codf
		)
	)
}

transcript_cumsum <- function(strand, start, end, include, txf, mc.cores=1L) {
	unsplit(
		f=txf,
		value=(if (mc.cores==1L) Map else function(...) mcMap(mc.cores=mc.cores, ...))(
			tapply(X=strand=="+", INDEX=txf, FUN="[", 1L),
			split(start, f=txf),
			split(end, f=txf),
			split(include, f=txf),
			f=function(str, s, e, inc) {
				m0 <- if (str) cbind(s, e) else -cbind(e, s)
				m <- m0-min(m0)
				if (min(m[,1]) != 0L) browser()
				stopifnot(all(m[,2]>m[,1]))
				stopifnot(min(m[,1]) == 0L)
				fl <- flatten(m[inc,,drop=FALSE])
				j <- join(fl, cbind(0, m[,2]))
				ifelse(inc, tapply(X=j[,4]-j[,3], INDEX=factor(j[,2], levels=seq_len(nrow(m))), FUN=sum), NA_integer_)
			}
		)
	)
}


remove_redundant_features <- function(type, start, end, fail_suspicious=FALSE) {
	stop("not implemented!")
}

rtf_base2rtf <- function(rtf_base, mc.cores) {
	stopifnot(setequal(names(rtf_base), c("chr","type","start","end","strand","phase","gene_id","transcript_id")))
	stopifnot(all(chr %in% chroms))
	mp <- if (mc.cores==1L) Map else function(...) mcMap(mc.cores=mc.cores, ...)
	rtf_base %>%
	mutate(
		strand=factor(strand, levels=c("+","-")),
		phase=as.integer(ifelse(phase==".", NA_integer_, phase))
	) %>%
	subset(type %in% RSVR_TRANSCRIPT_FEATURE_TYPES & !is.na(transcript_id) & !is.na(chr) & !is.na(type) & !is.na(strand)) %>%
	arrange(start) %>%
	mutate(end=end+1L) %>%
	subset(end > start) %>%
	mutate(first_cds=get_first_cds(type, start, end, strand, factor(transcript_id))) %>%
	(function(x) {
		data.frame(stringsAsFactors=FALSE, x, with(x, get_first_and_last_codon(type, chr, start, end, phase, strand, transcript_id, ref, mc.cores=mc.cores)))
	}) %>%
	subset({
		txf <- factor(transcript_id)
		is_exon <- type=="exon"
		unsplit(
			f=txf,
			value=simplify2array(mp(
				split(data.frame(s=start, e=end)[is_exon,], f=txf[is_exon]),
				f=function(tx_ex) {
					if (nrow(tx_ex) == 0L) TRUE
					else detached_sorted_nonempty(as.matrix(arrange(tx_ex, s)))
				}
			))
		)
	}) %>%
	subset({
		txf <- factor(transcript_id)
		is_cds <- type=="CDS"
		unsplit(
			f=txf,
			value=simplify2array(mp(
				split(data.frame(s=start, e=end)[is_cds,], f=txf[is_cds]),
				f=function(tx_ex) {
					if (nrow(tx_ex) == 0L) TRUE
					else {
						detached_sorted_nonempty(as.matrix(arrange(tx_ex, s))) & (with(tx_ex, sum(e-s)) %% 3L == 0)
					}
				}
			))
		)
	}) %>%
	mutate(
		cdna_sum=transcript_cumsum(
			mc.cores=mc.cores,
			strand,
			start, 
			end, 
			type=="transcript"|type=="exon"|type=="three_prime_utr"|type=="five_prime_utr"|type=="CDS",
			factor(transcript_id) 
		),
		cds_sum=transcript_cumsum(
			mc.cores=mc.cores,
			strand,
			start, 
			end, 
			type=="CDS",
			factor(transcript_id)
		)
	)
}

partition_tx_features <- function(tx_df) {
	stopifnot(is.data.frame(tx_df))
	copy_features <- c("strand", "gene_id", "transcript_id", "chr")
	spec_features <- c("start", "end", "type")
	blank_features <- setdiff(names(tx_df), c(copy_features, spec_features))
	stopifnot(all(c(spec_features, copy_features) %in% names(tx_df)))
	empty <- do.call(what=function(...) data.frame(stringsAsFactors=FALSE, check.names=FALSE, row.names=NULL, ...), lapply(tx_df, function(x) vector(mode=mode(x), length=0L)))
	if (nrow(tx_df) == 0L) {
		empty
	} else {
		stopifnot((function(x) all(x==x[1L]))(tx_df$strand=="+"))
		make_feats <- function(m, type) {
			s <- m[,1]
			e <- m[,2]
			if (length(s) == 0L) {
				empty
			} else {
				stopifnot((length(type)==1L) || (length(type)==length(s)))
				do.call(
					what=function(...) data.frame(
						stringsAsFactors=FALSE, 
						check.names=FALSE,
						start=s, 
						end=e, 
						type=type,
						...
					),
					c(lapply(tx_df[copy_features], "[", 1L), Map(tx_df[blank_features], f=function(x) { v <- NA; mode(v) <- mode(x); v }))
				)
			}
		}
		tec_df <- subset(tx_df, type %in% c("transcript", "exon", "CDS"))
		if (all(tec_df$type == "transcript")) {
			if (nrow(tec_df) != 1L) stop("non-contiguous transcript space!")
		} else if (any(tec_df$type == "exon") && !any(tec_df$type == "CDS")) {
			with(tec_df, {
				mex <- flatten(cbind(start, end)[type=="exon",,drop=FALSE])
				ex_region <- t(range(mex))
				mall <- t(range(c(start, end)))
				rbind(
					make_feats(setdiffs(mall, ex_region), "transcript"),
					make_feats(mex, "exon"),
					make_feats(setdiffs(ex_region, mex), "intron")
				)	
			})
		} else {
			positive <- (tx_df$strand=="+")[1]
			with(tec_df, {
				mall <- flatten(cbind(start, end))
				mex <- flatten(cbind(start, end)[type=="exon",,drop=FALSE])
				mcds <- flatten(cbind(start, end)[type=="CDS",,drop=FALSE])
				stopifnot(nrow(mcds) > 0L)
				mutrs <- setdiffs(mex, mcds)
				rbind(
					make_feats(setdiffs(mall, t(range(mex))), "transcript"),
					make_feats(setdiffs(t(range(mex)), mex), "intron"),
					subset(tec_df, type=="CDS"),
					make_feats(mutrs, ifelse((mutrs[,1]-mcds[1,1])*(if (positive) 1L else -1L) < 0, "five_prime_utr", "three_prime_utr"))
				)
			})
		}
	}
}

complete_features <- function(type, start, end, phase, positive, check_phase=TRUE) {
	stopifnot(all(type %in% c("exon", "transcript", "CDS")))
	if (length(type) == 0L) {
		data.frame(stringsAsFactors=FALSE, type=character(0), start=integer(0), end=integer(0), phase=integer(0))
	} else {
		if (all(type == "transcript")) {
			if (length(type) != 1L) stop("non-contiguous transcript space!")
			data.frame(
				stringsAsFactors=FALSE, 
				type=type,
				start=start,
				end=end,
				phase=phase
			)
		} else if (any(type == "exon") && !any(type == "CDS")) {
			mex <- flatten(cbind(start, end)[type=="exon",,drop=FALSE])
			ex_region <- t(range(mex))
			mall <- t(range(c(start, end)))
			ints <- list(setdiffs(mall, ex_region), mex, setdiffs(ex_region, mex))
			m <- do.call(what=rbind, ints)
			data.frame(
				stringsAsFactors=FALSE, 
				type=rep(c("transcript", "exon", "intron"), times=vapply(FUN.VALUE=integer(1), FUN=nrow, X=ints)),
				start=m[,1],
				end=m[,2],
				phase=rep(NA_integer_, nrow(m))
			)
		} else {
			mall <- flatten(cbind(start, end))
			mex <- flatten(cbind(start, end)[type=="CDS"|type=="exon",,drop=FALSE])
			exonic <- t(range(mex))
			mcds <- flatten(cbind(start, end)[type=="CDS",,drop=FALSE])	
			cds_lens <- mcds[,2]-mcds[,1]
			cum_len_at_start <- (function(x) { 
				f <- if (positive) identity else rev
				f(cumsum(f(x))-f(x))
			})(cds_lens)
			cds_df <- data.frame(s=start, e=end, p=phase)[type=="CDS",]
			mcd_phase <- with(cds_df, (3L - (((3L-p[if (positive) match(min(mcds[,1]), s) else match(max(mcds[,2]), e)]) + cum_len_at_start) %% 3L)) %% 3L)
			stopifnot(!any(is.na(mcd_phase)))
			if (check_phase) {
				if (!all(with(if (positive) list(orig=start, coll=mcds[,1]) else list(orig=end, coll=mcds[,2]), ifelse(orig %in% coll & type=="CDS", phase==mcd_phase[match(orig, coll)], TRUE))))
					stop("mismatch phase")
			}
			mutrs <- setdiffs(mex, mcds)
			ints <- list(tx=setdiffs(mall, exonic), intron=setdiffs(exonic, mex), cds=mcds, utrs=mutrs)
			do.call(what=rbind, Map(
				ints,
				with(ints, list(
					rep("transcript", nrow(tx)),
					rep("intron", nrow(intron)),
					rep("CDS", nrow(cds)),
					ifelse((utrs[,1]-cds[1,1])*(if (positive) 1L else -1L) < 0, "five_prime_utr", "three_prime_utr")
				)),
				list(NULL,NULL,mcd_phase,NULL),
				f=function(m, type, phase) {
					data.frame(stringsAsFactors=FALSE,
						type=type, 
						start=m[,1], 
						end=m[,2], 
						phase=if (is.null(phase)) rep(NA_integer_, nrow(m)) else phase
					)
				}
			))
		}
	}
}

complete_features_all_tx <- function(x, check_phase=TRUE, mc.cores=1L) {
	feature_specific <- c("start", "end", "type", "phase")
	txf <- factor(x$transcript_id)
	pos <- as.logical(tapply(X=x$strand=="+", INDEX=txf, FUN="[", 1L))
	properties_by_tx <- lapply(
		x[setdiff(names(x), c(feature_specific, "strand", "transcript_id"))],
		function(col) unname(tapply(X=col, INDEX=txf, FUN="[", 1L))
	)
	tabs <- mcMap(
		mc.cores=mc.cores,
		pos,
		split(x[feature_specific], f=txf),
		levels(txf),
		f=function(positive, df, lev) {
			with(df, complete_features(
				type=type,
				start=start,
				end=end,
				phase=phase,
				positive=positive,
				check_phase=check_phase
			))
		}
	)
	features_by_tx <- vapply(FUN.VALUE=integer(1), FUN=nrow, X=tabs)
	stopifnot(min(features_by_tx) > 0L)
	do.call(
		what=function(...) data.frame(
			stringsAsFactors=FALSE,
			check.names=FALSE,
			transcript_id=rep(levels(txf), times=features_by_tx),
			strand=factor(
				ifelse(rep(pos, times=features_by_tx),"+","-"),
				levels=c("+", "-")
			),
			...
		),
		c(
			lapply(setNames(nm=feature_specific), function(prop) do.call(what=c, lapply(tabs, "[[", prop))),
			lapply(properties_by_tx, function(prop) rep(prop, times=features_by_tx))
		)
	)
}

process_all_transcripts <- function(rtf_df, mc.cores=1L) {
	txf <- factor(rtf_df$transcript_id)
	rtf_by_tx <- split(rtf_df, f=txf)
	blocks <- split(rtf_by_tx, gl(n=mc.cores, k=1, length=nlevels(txf)))

	res <- do.call(what=rbind, 
		(if (mc.cores==1L) Map else function(...) mcMap(mc.cores=mc.cores, ...))(
			blocks, 
			f=function(block) do.call(what=rbind, lapply(block, partition_tx_features))))
	row.names(res) <- NULL
	res
}

gff2rsvr <- function(file, format=c("ensembl","refseq")) {
	stop("not implemented!")
}

reqd_feature_table_columns <- c("start","end","type","strand","gene_id","transcript_id","chr","phase","cdna_sum","cds_sum","first","last","transcript_biotype")

area <- function(keep, df, grp=factor(df$transcript_id)) { with(df, { use <- type %in% keep; mapply(FUN=function(s, e) { m <- flatten(cbind(s, e)); sum(m[,2]-m[,1]) }, split(start[use], f=grp[use]), split(end[use], f=grp[use])) } ) }

view_txt <- "CREATE VIEW COMPLETED_FEATURE AS SELECT F1.TYPE, CASE WHEN FEATURE_TYPE.NAME = 'intron' THEN F1.START-3 ELSE F1.START END AS ASTART, CASE WHEN FEATURE_TYPE.NAME = 'intron' THEN F1.END+3 ELSE F1.END END AS AEND, TX.STRAND, F1.PHASE, TX.GENE_ID, F1.TX_ID, F1.ID, CASE WHEN NOT F1.TYPE = 0 THEN NULL WHEN NOT LENGTH(F1.FIRST_CODON)=3 THEN F3.LAST_CODON || F1.FIRST_CODON ELSE F1.FIRST_CODON END AS COMPLETE_FIRST, CASE WHEN NOT F1.TYPE = 0 THEN NULL WHEN NOT LENGTH(F1.LAST_CODON)=3 THEN F1.LAST_CODON || F2.FIRST_CODON ELSE F1.LAST_CODON END AS COMPLETE_LAST, ((ROW_NUMBER() OVER (PARTITION BY F1.TX_ID, F1.TYPE ORDER BY CASE WHEN STRAND = 0 THEN -F1.START ELSE F1.START END)) = 1) AND (FEATURE_TYPE.NAME = 'CDS') AS FIRST_CDS, F1.CDNA_END, F1.CDS_END, CDNA_LENGTH, CDS_LENGTH, ((ROW_NUMBER() OVER (PARTITION BY F1.TX_ID, F1.TYPE ORDER BY CASE WHEN STRAND = 0 THEN F1.START ELSE -F1.START END)) = 1) AND (FEATURE_TYPE.NAME = 'CDS') AS LAST_CDS, CASE WHEN CAN.TX_ID IS NULL THEN 0 ELSE 1 END AS PRIORITY_TX FROM FEATURE F1 LEFT JOIN FEATURE F2 ON F1.TX_ID = F2.TX_ID AND F1.CDS_NUM = (F2.CDS_NUM-1) LEFT JOIN FEATURE F3 ON F1.TX_ID = F3.TX_ID AND F1.CDS_NUM = (F3.CDS_NUM+1) LEFT JOIN CANONICAL CAN ON F1.TX_ID = CAN.TX_ID, FEATURE_TYPE, TX, BIOTYPE WHERE F1.TYPE = FEATURE_TYPE.ID AND F1.TX_ID = TX.ID AND BIOTYPE.ID = TX.BIOTYPE_ID AND BIOTYPE.NAME IN ('IG_C_gene','IG_D_gene','IG_J_gene','IG_V_gene','lincRNA','miRNA','misc_RNA','Mt_rRNA','Mt_tRNA','protein_coding','rRNA','snoRNA','snRNA','TR_C_gene','TR_D_gene','TR_J_gene','TR_V_gene') ORDER BY ASTART"

create_tables <- c(
	"CREATE TABLE `GENE` ( `ID` INTEGER NOT NULL PRIMARY KEY, `NAME` VARCHAR(255) NOT NULL)",
	"CREATE TABLE `FEATURE_TYPE` ( `ID` INTEGER, `NAME` VARCHAR(255))",
	"CREATE TABLE `FEATURE` ( `ID` INTEGER NOT NULL PRIMARY KEY, `TYPE` INTEGER, `START` BIGINT, `END` BIGINT, `PHASE` INTEGER, `TX_ID` INTEGER, `FIRST_CODON` VARCHAR(3), `LAST_CODON` VARCHAR(3), `CDNA_END` INTEGER, `CDS_END` INTEGER, `CDS_NUM` INTEGER)",
	"CREATE TABLE `BIOTYPE` ( `ID` INTEGER NOT NULL PRIMARY KEY, `NAME` VARCHAR(255))",
	"CREATE TABLE `TX` ( `ID` INTEGER NOT NULL PRIMARY KEY, `NAME` VARCHAR(255), `BIOTYPE_ID` INTEGER, `STRAND` INTEGER, `GENE_ID` INTEGER, `CDNA_LENGTH` INTEGER, `CDS_LENGTH` INTEGER)",
	"CREATE TABLE `CANONICAL` ( `TX_ID` INTEGER NOT NULL PRIMARY KEY)"
)

create_indices <- c(
	"CREATE INDEX `FEATURE_TYPE-ID` ON FEATURE_TYPE(ID)",
	"CREATE INDEX `FEATURE-START` ON FEATURE(START)",
	"CREATE INDEX `FEATURE-TYPE` ON FEATURE(TYPE)",
	"CREATE INDEX `FEATURE-END` ON FEATURE(END)",
	"CREATE INDEX `FEATURE-TX_ID` ON FEATURE(TX_ID)",
	"CREATE INDEX `GENE-NAME` ON GENE(NAME)",
	"CREATE INDEX `TX-GENE_ID` ON TX(GENE_ID)",
	"CREATE INDEX `TX-BIOTYPE` ON TX(BIOTYPE_ID)"
)

create_database_from_rtf <- function(feature_table, target_db, canonical_tx=NULL) {
	if (file.exists(target_db)) unlink(target_db)
	stopifnot(all(reqd_feature_table_columns %in% names(feature_table)))
	gf <- factor(feature_table$gene_id)
	txf <- factor(feature_table$transcript_id)
	biotypes <- levels(factor(feature_table$transcript_biotype))
	db <- dbConnect(RSQLite::SQLite(), target_db)
	on.exit(dbDisconnect(db))
	for (ct in create_tables) dbExecute(db, ct)
	dbWriteTable(db, append=TRUE, overwrite=FALSE, name="GENE", value=data.frame(
		stringsAsFactors=FALSE, 
		ID=seq_len(nlevels(gf)), 
		NAME=levels(gf)
	))
	dbWriteTable(db, append=TRUE, overwrite=FALSE, name="BIOTYPE", value=data.frame(
		stringsAsFactors=FALSE, 
		ID=seq_along(biotypes), 
		NAME=biotypes
	))
	dbWriteTable(db, append=TRUE, overwrite=FALSE, name="TX", value=data.frame(stringsAsFactors=FALSE, 
		ID=seq_len(nlevels(txf)), 
		NAME=levels(txf), 
		BIOTYPE_ID=as.integer(factor(levels=biotypes, tapply(X=feature_table$transcript_biotype, INDEX=txf, FUN="[", 1L))),
		STRAND=as.integer(tapply(X=feature_table$strand=="+", INDEX=txf, FUN="[", 1L)),
		GENE_ID=as.integer(factor(levels=levels(gf), tapply(X=feature_table$gene_id, INDEX=txf, FUN="[", 1L))),
		CDNA_LENGTH=area(c("transcript","exon","three_prime_utr","five_prime_utr","CDS"), feature_table, txf),
		CDS_LENGTH=area("CDS", feature_table, txf)
	))
	dbWriteTable(db, append=TRUE, overwrite=FALSE, name="FEATURE_TYPE", value=data.frame(stringsAsFactors=FALSE, 
		ID=seq_along(RSVR_TRANSCRIPT_FEATURE_TYPES)-1L, 
		NAME=RSVR_TRANSCRIPT_FEATURE_TYPES
	))
	dbWriteTable(db, append=TRUE, overwrite=FALSE, name="FEATURE", value=with(feature_table, {
		chrval <- (as.integer64(2)^28)*match(chr, chroms)
		data.frame(stringsAsFactors=FALSE, 
			ID=seq_along(type), 
			TYPE=match(type, RSVR_TRANSCRIPT_FEATURE_TYPES)-1L, 
			START=start+chrval, 
			END=end+chrval, 
			PHASE=ifelse(is.na(as.integer(phase)), 0L, as.integer(phase)), 
			TX_ID=as.integer(txf),
			FIRST_CODON=first,
			LAST_CODON=last,
			CDNA_END=as.integer(cdna_sum),
			CDS_END=as.integer(cds_sum),
			CDS_NUM=unsplit(f=txf, value=Map(
				split(feature_table$start, f=txf),
				split(feature_table$type=="CDS", f=txf),
				f=function(s, cds) {
					cdsf <- factor(ifelse(cds, "cds", "not"), levels=c("cds","not"))
					unsplit(f=cdsf, list(cds=order(s[cds]), not=rep(NA_integer_, sum(!cds))))
				}
			))
		)
	}))
	for (ci in create_indices) dbExecute(db, ci)
	if (!is.null(canonical_tx)) {
		dbWriteTable(db, name="CANONICAL", append=TRUE, overwrite=FALSE, value=data.frame(TX_ID=match(canonical_tx, levels(txf))))
	} else {
		dbExecute(db, "INSERT INTO CANONICAL SELECT ID FROM (SELECT ID, ROW_NUMBER() OVER (PARTITION BY GENE_ID ORDER BY CDS_LENGTH DESC, CDNA_LENGTH DESC) AS RN FROM TX) X WHERE RN = 1")
	}
	dbExecute(db, view_txt)
}

suppressPackageStartupMessages(library(IntervalSurgeon))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(bit64))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly=TRUE)
if (!(length(args)) %in% 2:3)
	stop("Usage: Rscript gtf2featuresDB.R <gtf> <ref fasta> <target file> <canonical transcripts: optional>")
gtf <- args[1]
ref <- args[2]
target_db <- args[3]
canonical_file <- args[4]

base_df <- gtf2rtf_base(gtf)

completed_features_df <- 
	mutate(
		base_df,
		end=end+1L,
		strand=factor(strand, levels=c("+","-")),
		phase=as.integer(ifelse(phase==".", NA_integer_, phase))
	) %>%
	subset(!is.na(type) & type %in% c("CDS", "exon", "transcript")) %>%
	subset(!is.na(strand) & !is.na(chr) & !is.na(transcript_id) & !is.na(gene_id)) %>%
	(function(x) { stopifnot(with(x, all(end > start))); x }) %>%
	complete_features_all_tx(check_phase=FALSE)

with_fl_codons <- 
	completed_features_df %>%
	(function(x) data.frame(stringsAsFactors=FALSE, x, with(x, get_first_and_last_codon(type, chr, start, end, phase, strand, transcript_id, ref, mc.cores=1L))))

ensembl_tab <- mutate(
	with_fl_codons,
	cdna_sum=transcript_cumsum(
		mc.cores=1L,
		strand,
		start, 
		end, 
		type=="transcript"|type=="exon"|type=="three_prime_utr"|type=="five_prime_utr"|type=="CDS",
		factor(transcript_id) 
	),
	cds_sum=transcript_cumsum(
		mc.cores=1L,
		strand,
		start, 
		end, 
		type=="CDS",
		factor(transcript_id)
	)
)

create_database_from_rtf(
	target_db=target_db,
	feature_table=ensembl_tab,
	canonical_tx=if (!is.na(canonical_file)) readLines(canonical_file)
)
