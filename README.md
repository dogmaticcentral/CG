# CG

CG is a set of scripts for extracting information about (localized) somatic 
mutations from 
[TCGA](https://portal.gdc.cancer.gov/)
and 
[ICGC](https://dcc.icgc.org/)
databases. The icgc branch also contains tools necessary to
merge the info from the two sources. The scripts add up to a loosely connected pipeline 
(or, rather, two pipelines).  They re-organize the data and store it in the local MySQL 
database.  While the  back end of each pipeline is rather generic, the front end is
geared toward answering particular questions for which they were originally written.

CG is not an out-of-the box solution. Rather, it is a starter kit, in case you would like to
do some cancer genomics data analysis on your own. Installing CG database(s) may take a
day or two if you are willing to go with the pipeline as-is. The (approximate) timings
for a default installation of the ICGC branch can be found in [timing.txt](icgc/timing.txt).
 With tweaks, a week is 
not an unreasonable time estimate.


Why bother with the local copy of the data? I do not know a general answer to that question.
You should check the homepage  of both databases - maybe the information you
 are looking for is already  available.
Here the original motive was to study the co-occurrence of mutations in different genes. 
These types of questions are not readily answerable using data portals - portals  tend to
agglomerate data on per-gene basis, in order to protect the privacy of sample donors. 


<!-- this is a comment -->
<!-- making TOC: https://github.com/ekalinin/github-markdown-toc -->
<!-- once installed, use with gh-md-toc README.md    -->
## Table of Contents
* [Dependencies](#dependencies)
* [TCGA](#tcga)
* [ICGC](#icgc)
     * [config file](#config-file)
     * [ICGC data download ](#icgc-data-download-00_data_download)
     * [Loading data into local version of the database](#loading-data-into-local-version-of-the-database-10_local_db_loading)
     * [Reorganizing mutation data](#reorganizing-mutation-data-20_local_db_reorganization)
     * [Merging with TCGA](#merging-with-tcga-30_tcga_merge)
     * [Housekeeping](#housekeeping-40_housekeeping)
     * [Production](#production-50_production)
* [TODO](#todo)
  
 
 
## Dependencies
 In addition to TCGA and ICGC databases themselves, CG relies on
 * MySQL
 * MySQLdb python module, installed with _sudo apt install python3-mysqldb_
 * gene symbols from HUGO gene nomenclature committee (see [here](https://www.genenames.org/download/custom/))
 * [Annovar](http://annovar.openbioinformatics.org/en/latest/) for location and functional annotation - int TCGA merge
 * CrossMap  (in TCGA merge) - (sudo pip3 install CrossMap, pyBigWig, pysam)
 * Optional: [line-profiler](https://github.com/rkern/line_profiler#line-profiler) for python; gcc compiler for c-utils
  used in postprocessing.
 
## TCGA
 The tcga branch of the icgc pipeline got obsoleted before coming to production stage. 
 It contains various blind a alleys and wrong turns. Its current use is as a prep
 step for merging with the icgc branch. It has [its own README page](tcga).
 
## ICGC
 
 A general note: throughout the pipeline, you will find scripts disabled by having exit() right on the top of the file.
 These are the scripts that drop tables and/or store without checking. Enable them by commenting the exit line.
 (The advice is to put the comment back in once the script is done.)
 
### config file
 You can set some recurring constants - such as data directories or mysql conf file(s) - 
 in the [config.py](icgc/config.py) file.
 
<!-- **********************************************************  -->
<!-- **********************************************************  -->
<!-- **********************************************************  -->
### ICGC data download ([00_data_download](icgc/00_data_download))
 
 Just like the tcga branch, this branch of the pipeline starts by downloading the data from
 the source, ICGC in this case. Note however that here you will need the access token. 
 Obtaining one is a lengthy process (as in several weeks to several months) which you can
 start investigating [here](https://icgc.org/daco).
 
 One you have the access token,  place it in the environmental variable called ICGC_TOKEN, to make
 the download scripts work.
 
 Note in particular that we are grouping some cancers under the same head-group. 
 See [06_group_cancers.py](icgc/00_data_download/06_group_cancers.py). This is because different depositors may
 use different shorthands for the same cancer (e.g. LICA == 'Liver Cancer', 
 LINC == 'Liver Cancer - NCC', LIRI == 'Liver Cancer - RIKEN'), though in some cases it 
 might not be clear which cancer the depositors refer to. Feel free to change in you version of the code
 the grouping defined in [06_group_cancers.py](icgc/00_data_download/06_group_cancers.py), or to skip it altogether.
 
 
### Loading data into local version of the database ([10_local_db_loading](icgc/10_local_db_loading))
 
 Make sure you have 
 the mysql conf file  and set its path in these two scripts, or arrange some other way to
 access the local database. The last I checked, python's MySQLdb package did not work with
 the encripted cnf files, so the only alternative is using 
 [client field in generic mysql option file](https://dev.mysql.com/doc/refman/8.0/en/option-files.html),
 like this:
 
`[client]`   
`user = blah`  
`host = localhost`  
`password = "somepasswd"`

In MySQL shell (or however you talk to your MySQL server) create the database icgc and the user _blah_, 
and give _blah_  the permissions to write to and read from _icgc_:

`create database icgc;`    
`create user 'blah'@'localhost' identified by 'somepasswd';`  
`grant all privileges on icgc.* to 'blah'@'localhost';`  
`flush privileges;`  


#### Measuring the field lengths and making MySQL tables
[05_find_max_field_length.py](icgc/10_local_db_loading/05_find_max_field_length.py) and 
[06_make_tables.py](icgc/10_local_db_loading/06_make_tables.py): 
Make sure that the fields in the mysql tables are big enough 
for each entry and create mysql tables. 
[05_find_max_field_length.py](icgc/10_local_db_loading/05_find_max_field_length.py) should 
give you an idea about the longest entries found.

#### Filling and indexing the icgc database tables
[07_write_mutations_tsv.py](icgc/10_local_db_loading/07_write_mutations_tsv.py) through 
[10_make_indices.py](icgc/10_local_db_loading/10_make_indices_on_temp_tables.py).
For large tables, rather than loading them through python, 
it turns out to be faster to create tsvs and  then load them from mysql shell 
(as in [09_load_mysql.py](icgc/10_local_db_loading/09_load_mysql.py); alternative: use mysqlimport manually) 
 to read them in wholesale. These scripts take care of that part, plus some index creating on the newly loaded tables.
 Make sure to run [10_make_indices_on_temp_tables.py](icgc/10_local_db_loading/10_make_indices_on_temp_tables.py), 
 [10_reorganize_variants.py in 20_local_db_reorganization](icgc/20_local_db_reorganization/10_reorganize_variants.py)
 pretty much does not work without it at all. 
 All index making is slow here (see [timing.txt](icgc/timing.txt)) - run overnight. 

Some checks are thrown in here that  may inform the rest of the pipeline.
 [15_assembly_check.py](icgc/10_local_db_loading/15_assembly_check.py) confirms that, as of v27, all ICGC entries
 refer to GRCh37, and [16_donor_check.py](icgc/10_local_db_loading/16_donor_check.py) highlights the fact that
 some donor ids have no somatic mutations in ICGC. This is somewhat mysterious, because some refer to
 TCGA donors with somatic mutation data available from TCGA archive.


#### Getting and storing some auxilliary data
 
In [20_hgnc_name_resolution_table.py](icgc/10_local_db_loading/20_hgnc_name_resolution_table.py) and 
 [22_ensembl_id_table.py](icgc/10_local_db_loading/22_ensembl_id_tables.py) we make and fill some tables we will use later for name resolution 
 (translating between gene and protein names used in different contexts).

The annotation across different submitters to TCGA/ICGC is not uniform In particular, for the missense mutations
sometimes it is not clear which splice they refer to. Alternatively, when the reference splice(s) is listed, it
si not clear which splice is the canonical splice. To remedy that, here we do some basic annotation of our own. For that
we will need the coding sequence of canonical transcripts.
 
The canonical transcript id is not readily available from Ensembl Mart, thus for our purposes
here you can find this info in the tables called ensembl_gene2trans_stable.tsv.bz2 
and in ensembl_deprecated2new_id.tsv.bz2 the
[hacks](icgc/hacks) directory. Decompress them (bzip2 -d) and put someplace where
[20_hgnc_name_resolution_table.py](icgc/10_local_db_loading/20_hgnc_name_resolution_table.py) can find them.

From the same place (Ensembl Mart) all human coding sequences can be downloaded in fasta format.
You can reduce that file to sequences of canonical transcript only. This is left as na exercise for the reader
 (hint: use  ensembl_gene2trans_stable.tsv and
   [blastdbcmd](https://www.ncbi.nlm.nih.gov/books/NBK279689/)
  tool; blastdbcmd has batch mode (blastdbcmd -h)). When you are happy with your fasta file, put is somewhere 
  where [23_ensembl_coding_seqs.py](icgc/10_local_db_loading/23_ensembl_coding_seqs.py) can find it.
  Note that this is optional: if you are happy wihtout knowing the position of the missense mutant on
  the canonical translation, you can move on without this step.


ICGC database does not have a complete  consensus on location annotation,  so we will be doing it ourselves.
As a prep, we download gene coordiantes from UCSC. (The coordinates are actually from Ensembl, but UCSC 
keeps it in a format that is more readily usable.) The script is
[24_ucsc_gene_coords_table.py](icgc/10_local_db_loading/24_ucsc_gene_coords_table.py). To download coordinates
from their MySQl server you will need an internet connection, mysql client, and another conf file, like this:

`[client]`  
`skip-auto-rehash`   
`user = genome`   
`host = genome-mysql.soe.ucsc.edu`   
`port = 3306`   

Again, the path to that file is expected to be defined in the [config.py](icgc/config.py) file.


 
<!-- **********************************************************  -->
<!-- **********************************************************  -->
<!-- **********************************************************  -->
### Reorganizing mutation data ([20_local_db_reorganization](icgc/20_local_db_reorganization))

This is where we depart from ICGC original database architecture - which is pretty much
nonexistent and consists of massive duplication of annotation for each occurrence of a mutation
and for each of its interpretations within various transcripts.

So instead we reorganize the database into something like this 
![db schema - schematic](illustrations/schema_schematic.png)
where *_specimen, \*\_donor, and \*\_simple_somatic tables exist for each cancer type, and mutations\_\* 
and locations\_\* tables exist for each chromosome.


<!-- to produce the schema visualization
java -jar ~/Downloads/schemaSpy_5.0.0.jar  -t mysql -host localhost  -db icgc  \
-u usrnm -p passwd -o icgc_schema  
-dp ~/Downloads/mysql-connector-java-5.1.6/mysql-connector-java-5.1.6-bin.jar 
where  icgc_schema is output dir
schemaSPy: http://schemaspy.sourceforge.net/
mysql-connector-java:  https://dev.mysql.com/downloads/connector/j/5.1.html
-->

#### Creating new tables

In [06_consequence_vocab.py](icgc/20_local_db_reorganization/06_consequence_vocab.py)
 we inspect the 'consequence' vocabulary employed by ICGC. There seems to
some confusion there about the location vs. the consequence of a mutation.  This info is used
in [13_reorganize_mutations.py](icgc/20_local_db_reorganization/13_reorganize_mutations.py)  to 
come up with the pathogenicity estimate, to be stored in the eponymous field in the
mutations\* tables.


New tables are created in [08_check_mut_tables_and_make_new_ones.py](icgc/20_local_db_reorganization/08_check_icgc_tables_and_make_new_ones.py).

In [10_reorganize_variants.py](icgc/20_local_db_reorganization/10_reorganize_variants.py),
[11_delete_variants_from_normal.py](icgc/20_local_db_reorganization/11_delete_variants_from_normal.py),
[13_reorganize_mutations.py](icgc/20_local_db_reorganization/13_reorganize_mutations.py),   and
[14_reorganize_locations.py](icgc/20_local_db_reorganization/14_reorganize_locations.py) 
 you can choose to run in parallel (the number of 'chunks' in main()). 
 
 
**Note 1:** that in [11_delete_variants_from_normal.py](icgc/20_local_db_reorganization/11_delete_variants_from_normal.py),
 we are dropping variants from normal samples. This is something you might not want to do if you are trying to
 annotate variants yourself. Though in that case you might want to go back to 
 [08_check_mut_tables_and_make_new_ones.py](icgc/20_local_db_reorganization/08_check_icgc_tables_and_make_new_ones.py)
 and keep the 'matched_icgc_sample_id' field. 
 
**Note 2:** doing things carefully leads to some interesting results. NACA (Nasopharyngeal cancer) set, for example,
consists of normal tissue samples only. 

[14_reorganize_locations.py](icgc/20_local_db_reorganization/14_reorganize_locations.py) script uses 
UCSC gene annotation to check chromosome addresses. The only information we are looking for here is the
possibility that the location falls within the splice region just outside of an exon. Mutations at these positions
are annotated as (possibly) pathogenic 
in [15_update_pathogenicity_in_mutation_tables.py](icgc/20_local_db_reorganization/15_update_pathogenicity_in_mutation_tables.py).
In [16_pathg_from_mutations_to_variants.py](icgc/20_local_db_reorganization/16_pathg_from_mutations_to_variants.py) we are adding this
same info (it is a boolean flag, no big elaboration on the table) to the variantis (*_simple_somatic) tables.

    
(Do not forget to create indices
using [10_make_indices_on_temp_tables.py](icgc/old/10_make_indices_on_temp_tables.py)) 
 
 
 #### Adding reliability info
 
 We add a couple of values to each row to later make the search for meaningful entries faster.
  we are adding mutant_allele_read_count/total_read_count ratio and pathogenicity estimate (boolean)
 to simple_somatic tables. These scripts can be run in the order indicated in the name, rather than at this point.
 In the following script,  
 [21_add_realiability_annotation_to_somatic.py](icgc/20_local_db_reorganization/21_add_reliability_annotation_to_variants.py),  
 we combine these two columns into a reliability estimate: a  somatic mutation in individual patient is considered reliable if mutant_allele_read_count>=10
 and mut_to_total_read_count_ratio>=0.2. Information about the mutation in general (mutations_chromosome tables;  
 [27_reliability_from_variants_to_mutations.py](icgc/20_local_db_reorganization/27_reliability_from_variants_to_mutations.py)) 
 is considered reliable if there is at least one patient for which it was reliably established.


#### Removing duplicates
 ICGC is rife with data duplication, coming from various sources. Some seem to be bookkeeping mistakes with the
 same patient data finding its way into the dataset through various depositors; some are the results  of the re-sampling 
 of the same  tumor, while some are completely obscure, with all identifiers being identical everywhere.
 
**All identifiers identical.**
  [18_cleanup_duplicate_entries.py](icgc/20_local_db_reorganization/18_cleanup_duplicate_entries.py):
  Some mutations  have identical tuple
 of identifiers (icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id). Note that this
 is after we have reorganized the database so that the mutation and location info sit in 
 different tables from the donor info. Not sure what this is about (the same sample analyzed independently multiple
 times?), but when a duplicate is found, this script chooses the entry with the greatest depth reported. 
 See the script for the full resolution strategy and make sure to run 
 [17_make_jumbo_index](icgc/20_local_db_reorganization/17_make_jumbo_index_on_new_tables.py) 
  beforehand.
 
**Same submitted_sample_id, different donor ids.**
 There might be further problems: See for example, mutation MU2003689, which, 
 [so the ICGC page claims](https://dcc.icgc.org/mutations/MU2003689) can be found in two distinct donors. 
 The closer  inspection of the two donors shows however that their submitter ID is the same, as is the age 
 of the 'two' women. (The tumour subtype has different description, reflecting, apparently,  the
 curator's preference.) Indeed, donors table for BRCA, at this point in the pipeline has 
 1976 distinct ICGC donor ids, and 1928 distinct submitter IDs. BRCA does turn out to be the biggest offender here,
 followed by LICA with 8 duplicated donors. It is not clear whether these duplicates refer to the same
 tumor at the same stage because even the submitter sample ids might be different
 (see [19_cleanup_multiple_donor_for_the_same_submitted_id.py](icgc/20_local_db_reorganization/19_cleanup_multiple_donor_for_the_same_submitted_id.py)). 

**Same donor with differing specimen and sample ids.**   Apparently  ICGC refers to biological replicates - i.e. 
 samples taken from different sites  - as specimens, and
 to technical replicates as samples. Perhaps they might have a role when answering different
 types of questions than what we have in mind. Here, however we do not want to have these results mistaken for recurring mutations, 
 thus we remove them in [22_cleanup_duplicate_specimens.py](icgc/20_local_db_reorganization/22_cleanup_duplicate_specimens.py) 
 and [23_cleanup_duplicate_samples.py](icgc/20_local_db_reorganization/23_cleanup_duplicate_samples.py), but not before checking
 which of the samples produced more reliable reads (see below).
 In this version of the pipeline we keep only the sample annotated as 'Primary tumour - solid tissue.' Out of
 these, if multiple refer to the same submitter id, we keep the ones with the largest reported number of
 somatic mutations. The investigation of the source of this duplication is again outside of our zone of interest.
 
**Different donor and submitter ids with overlapping variants.**
 Unfortunately, the duplicates do not stop here. See for example [DO224621](https://dcc.icgc.org/donors/DO224621)
 and [DO230968](https://dcc.icgc.org/donors/DO230968) that have different  donor *and* submitter
 ids, and mysteriously have 1703 identical variants. 
 Both donors are males diagnosed with lung squamous carcinoma. One sample turns out to be WES, while the other is WGS.
 
[25_cleanup_duplicate_donors_by_variant_overlap.py](icgc/20_local_db_reorganization/25_cleanup_duplicate_donors_by_variant_overlap.py)
attempts to detect such cases by looking for suspiciously high overlap in reported variants. When such case is detected,
 the sample with the higher average depth of sampling is retained.
 See the script for the criteria used. 
 
 #### Mutation2gene shortcut
 
 [29_index_on_mutation_tables.py](icgc/20_local_db_reorganization/29_index_on_mutation_tables.py) 
 and [30_mutation2gene_maps.py](icgc/20_local_db_reorganization/30_mutation2gene_maps.py)
 together build a table that maps icgc_mutation_id to HUGO gene
 symbol, and vice versa.
 
 #### Optional: re-annotating missense mutations
 It may be  helpful to change the annotation of missense mutations
 ([34_reannotate_missense_mutations.py](icgc/20_local_db_reorganization/34_reannotate_missense_mutations.py)) 
 to include the currently accepted  canonical transcript, according to [Ensembl](https://www.ensembl.org).
 
 After running reference [34_reannotate_missense_mutations.py](icgc/20_local_db_reorganization/34_reannotate_missense_mutations.py),
  the reference will be lost to other (non-canonical) transcripts, which can be retrieved from the locations table.
  The prerequisite is [32_load_ensembl_coord_patches.py](icgc/20_local_db_reorganization/32_load_ensembl_coord_patches.py), 
  a hacky solution to include the latest gene coordinates from Ensembl. The coordinate tables can be found in 
  [hacks/coord_patches.tar.bz2](icgc/hacks/coord_patches.tar.bz2).
 A step toward independent annotation. 
 
 #### ICGC-only  production
 It is possible to stop here and move to production scripts, if there is no interest in including TCGA.
 
 
### Merging with TCGA ([30_tcga_merge](icgc/30_tcga_merge))

#### Disaster recovery strategy
 It might be advisable to backup the newly-created tables, by storing them as an sqldump, for example.
 There is a couple of scripts in [hacks](icgc/hacks) directory to help along.
 [somatic_tables_dump.pl](hacks/somatic_tables_dump.pl) will dump them out, just make sure you move to the storage direcotry, 
 and [somatic_tables_load.pl](hacks/somatic_tables_load.pl) will load them back in if needs be.
 If you choose to do the full database dump, 
 [somatic_tables_from_dump.pl](hacks/somatic_tables_from_dump.pl) can extract only *_simple_somatic tables, 
 the process if slow, however.
 
#### Merging
 The scripts [29_index_on_mutation_tables.py](icgc/20_local_db_reorganization/29_index_on_mutation_tables.py)
 through [34_tcga_specimen_hack.py](icgc/30_tcga_merge/37_tcga_specimen_hack.py)
 concern themselves with merging TCGA info created in TCGA branch with the ICGC.
 This sub-pipe runs from preparatory indexing to data input.  
 [Duplicate data removal](#removing-duplicates) should be probably be re-applied 
 (steps [19_cleanup_duplicate_donors.py](19_cleanup_duplicate_donors.py) 
 and [22_cleanup_duplicate_specimens.py](22_cleanup_duplicate_specimens.py)).
 If everything is ok, [35_database_stats.py](icgc/40_housekeeping/37_database_stats.py) should report
 no duplicates in any of the tables.
 
### Housekeeping ([40_housekeeping](icgc/40_housekeeping))
 
### Production ([50_production](icgc/50_production))
 
 
## TODO 
* disentangle from Annovar - we have all the info we need to do own annotation here
* why does deleting normal samples take so long?
( [11_delete_variants_from_normal.py](icgc/20_local_db_reorganization/11_delete_variants_from_normal.py))
* when re-annotating we are looking only at positions already labeled as having aa_change
* note to myself: there seems to be a problem in ENSEMBL homo_sapiens_core_94_38 with exon frame/phase assignment
e.g. ENST00000435165 (ENSG00000118473), chrom 1, + strand, the translation starts at 67161031.
The first exon is [67160376, 67161176> yet the pahse in ensembl is listed as -1, which should indicate that
translation does not start in this exon at all. 
Alternatively, ENST00000435165 is listed to start in ENSE00001951262, which has phase listed as -1.
Another example: ENST00000334103	ENSG00000186094, frames listed in Ensembl vs UCSC.
I have not used the phase info here, but it should be kept in mind that it is unreliable.
* annotation in the production stage is incomplete, possibly because of the mismatch between TCAG/ICGC ref assembly and
GRCh28 that i currently the standard with Ensembl
* there might be some errors in the re-annotation \
[32_reannotate_positions.py](icgc/20_local_db_reorganization/34_reannotate_missense_mutations.py) - \
 coding sequences from biomart not correct (e.g ENST00000366645 vs the same seq on the Ensembl website
[the same seq on the Ensembl website](https://grch37.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000116903;r=1:231468499-231473598;t=ENST00000366645))


## P.S.
Do not use locks in MySQL. Locks are evil. Good luck with the rest.