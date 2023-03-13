# rsvr

The `rsvr` program is a set of tools for operating on RSVR IDs and building the Rareservoir database.
It can be built and installed by running `make` and `make install` in the project directory (note: tested for GCC version 11 and clang version 14).
Each individual tool can be invoked by running `rsvr <tool>` (see `rsvr -h` for usage information).

The core of a Rareservoir database comprises tables: `GENOTYPE` (containing triples RSVR ID, sample ID and genotype), `VARIANT` (RSVR IDs and corresponding annotations acquired from large source files, in principle containing records for all variants in the `GENOTYPE` table), `SAMPLE` (containing database IDs linking to the `GENOTYPE` table and optional sample metadata such as sex, ancestry, etc.), `CONSEQUENCE` (Sequence Ontology annotation for RSVR ID/transcript pairs), and complementary tables of genes, transcripts and biotypes. Users may add additional tables and columns to suit their requirements and certain applications may impose particular contraints (e.g. statistical analysis of unrelated samples requires a column in the SAMPLE table indicating membership to an unrelated set).

The procedure for building a complete Rareservoir thus depends on the target schema. The procedure for building the core rareservoir additionally depends on the format of the original genotype data, which may be in another database, a file containing variant calls for a single sample, or one containing calls for all samples in a cohort. Additionally, parameters like allele frequency threshold, inclusion of variants depending on predicted transcript consequences, availability of CPU resources, etc. maybe considered. 

The `rsvr` program and complementary scripts in this repository are designed to provide the flexibility to compose a build procedure tailored to a broad range of requirements. An example build of the core database is provided in the `scripts/example-build.sh` script, which builds a Rareservoir from a merged VCF. Supplementary tables or extensions to core tables may be added at a user's discretion.

## Building the core Rareservoir database

The core database contains the following tables: 
* `GENOTYPE`
* `VARIANT`
* `SAMPLE`
* `CONSEQUENCE`
* `GENE`
* `TX`
* `BIOTYPE`
* `FEATURE`
* `CONSEQUENCES_LIST`

See `db/rareservior-schema.sql` for complete table specifications.

### Building the `GENOTYPE` table

The method of construction of the `GENOTYPE` table depends on the format of the genotype data. Once constructed, the subsequent build procedure for the core database is the same. The format of the `GENOTYPE` table which may be used as the input of the subsequent build procedure is a SQLite database file containing one table with columns:
* `RSVR_ID`: a 64-bit integer containing the RSVR ID of a variant,
* `SAMPLE_ID`: a 32-bit integer containing the ID of a sample,
* `GT`: an 8-bit integer containing the genotype of the sample at the variant site.
The table must additionally be indexed by `RSVR_ID` (and optionally `SAMPLE_ID`, which is recommended for the scenario when the database will be used to perform sample-wise look-ups).

Variant data may be piped from the original source through programs to manipulate it into the requisite SQL commands (i.e. `INSERT INTO GENOTYPE VALUES(...); INSERT INTO...`). `rsvr enc` should be used to convert variants formatted as CHROM, POS, REF, ALT into RSVR IDs: see `rsvr -h` for other relevant programs. Output can then be directed into `sqlite3` against a target `.db` file which is then used as the input for the next stage of the build procedure. The `scripts/build-GENOTYPE.sh` script contains an example script for building the table from a merged VCF.

### Building the `VARIANT` table

The `VARIANT` table links RSVR IDs to values from large files/databases and should represent all variants present in the `GENOTYPE` table. Here we use CADD Phred scores and PMAF values (probabilistic MAF based on gnomAD). It may also optionally contain a `REF` column for storing the reference allele for each variant to facilitate simple, rapid decoding of RSVR IDs to full CHROM/POS/REF/ALT tuples. The `scripts/anno.sh` script may be used to annotate an input stream of RSVR IDs with these values for a given reference genome, CADD Phred scores tsv file (https://cadd.gs.washington.edu/download) and gnomAD VCF file. The `scripts/build-VARIANT.sh` script contains an example script for building the `VARIANT` table from a .db file containing the `GENOTYPE` table.


### Building the `CONSEQUENCE` table

The `CONSEQUENCE` table (and supplementary tables of genes, transcripts, etc.) is created by annotating the RSVR IDs from the `VARIANT` or `GENOTYPE` tables against transcript features from a GTF file of transcripts with the `rsvr seqfx` program (see `rsvr seqfx -h` and the build script `scripts/build-CONSEQUENCE.sh`). The columns of the table corresponding are:

* `RSVR_ID`: A 64-bit RSVR ID for the annotated variant.
* `GENE_ID`: An integer linking to the `ID` column of the `GENE` table of genes.
* `TX_ID`: An integer linking to the `ID` column of the `TX` table .
* `CDNA`: Position of the variant in the transcript cDNA.
* `CDS`: Position of the variant in the coding sequence.
* `WORST_CSQ_ID`: 64-bit CSQ ID for the consequences against the most severely affected transcript.
* `ALL_CSQ_ID`: 64-bit CSQ ID for any consequences against the gene.
* `CANONICAL_CSQ_ID`: 64-bit CSQ ID for the consequences against the canonical transcript (0 if canonical transcript unaffected).

The `scripts/build-CONSEQUENCE.sh` script contains an example script for building the `CONSEQUENCE` table from a .db file containing the `VARIANT` table.

## Adding the `SAMPLE` table and phenotype data

The `SAMPLE` table should link the sample IDs used in the `SAMPLE_ID` column of the `GENOTYPE` table to sample names (the `NAME` column) and sex (`SEX`). Optionally additional columns may be added. 

In order to extend the Raresevoir database to include phenotypic data, we recommend adding tables:
* `PHENOTYPE`: a table of phenotypes belonging to various standard phenotype dictionaries (e.g. the Human Phenotype Ontology - HPO, International Classification of Diseases - ICD, ...) or custom phenotype dictionaries.
* `SAMPLE_PHENOTYPE`: a table linking phenotypes in the `PHENOTYPE` table (by type/code pair) to samples in the `SAMPLE` table (by ID). Matching to the `PHENOTYPE` table in this way ensures that if the `PHENOTYPE` table is updated, the `SAMPLE_PHENOTYPE` doesn't need updating.
* `CASE_SET` (where applicable): a table assigning labels corresponding to diseases/potential aetiological disease groups to samples, linking to the `SAMPLE` table by sample ID. See `scripts/case-set.R` for example import of a set of such aetiological groups.

Templates for these tables are available in the `db/schema.sql` file. The script `scripts/samples.R` may be used to import a tab-delimited file of sample information into the DB.

## Running the example build

`example-build.sh` requires a config file passed as an argument specifying the following variables (as x=y):

* `REF` location of the reference genome fasta file.
* `VCF` location of a merged VCF file.
* `ACLIM` maximum internal allele frequency (typically 0.2% of total allele number).
* `GNOMAD` location of the gnomAD VCF file.
* `BUILD` genome build version, either v37 or v38.
* `TARGET` target SQLite database to be built.
* `CADD` location of the genome wide CADD score file for SNVs (see https://cadd.gs.washington.edu/).
* `GTF` must be set to the location of a GTF file containing the transcript features for which consequences should be predicted for each variant.

Note that the `SAMPLE` table should be added such that samples in the merged VCF are added as rows with IDs corresponding to their position in the merged VCF (typically the `NAME` values would be those given in the merged VCF).

## Building Rareservoir on a database server

A Rareservoir database maybe built from scratch on a database server, or generated as a SQLite file and then imported. Once built on a database server, convenience functions to simplify writing queries (including transforming between VCF format and RSVR ID format and converting CSQ IDs to Sequence Ontology names) maybe created using the `stored-functions.sql` script in the `db` folder.

## Rareservoir database queries

Rareservoir tables typically store RSVR IDs in a column called `RSVR_ID` hence querying these tables together requires joining on the `RSVR_ID` columns. The `GENOTYPE` table can be joined with the `SAMPLE` table by the `SAMPLE_ID` and `ID` columns respectively.

This example query demonstrates how you can select genotype data, including sample and variant identifiers and variant annotation, depending on predicted consequence and for samples with a particular phenotype (note that the table `CONSEQUENCES_LIST`/stored function `CSQ_NAME_TO_ID` can be used to find the bit in the 64-bit consequence number ('CSQ ID') corresponding to a particular Sequence Ontology term name - see here how this is used to threshold results on 'missense or more severe').

```
SELECT 
    S.NAME,
    CHROM(G.RSVR_ID) AS CHROM,
    POS(G.RSVR_ID) AS POS,
    REF(G.RSVR_ID) AS REF,
    ALT(G.RSVR_ID) AS ALT,
    G.GT,
    V.CADD_PHRED AS CADD,
    CONSEQUENCES(C.ALL_CSQ_ID) AS SEQ_ONT_CONSEQUENCES
FROM VARIANT V,
    GENOTYPE G, CONSEQUENCE C,
    SAMPLE S, SAMPLE_PHENOTYPE SP, PHENOTYPE P,
    GENE GN
WHERE G.SAMPLE_ID = S.ID 
    AND V.RSVR_ID = G.RSVR_ID
    AND G.RSVR_ID = C.RSVR_ID 
    AND GN.ID = C.GENE_ID
    AND S.ID = SP.SAMPLE_ID 
    AND SP.PHENOTYPE_TYPE = P.TYPE 
    AND SP.PHENOTYPE_CODE = P.CODE 
    AND P.NAME = 'Thrombocytopenia'
    AND C.ALL_CSQ_ID >= CSQ_NAME_TO_ID('missense_variant')
```

## Citation

If you use the `rsvr` package in your published research, please cite:
Greene et al., Genetic association analysis of 77,539 genomes reveals rare disease etiologies, Nature Medicine (2023), https://doi.org/10.1038/s41591-023-02211-z
