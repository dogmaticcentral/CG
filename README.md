# CG

_(writeup in progress)_

CG is a set of scripts for extracting information about somatic 
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
for a default installation of the ICGC branch can be found in [timing.txt](timing.txt). With tweaks, a week is 
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
    * [ICGC data download](#icgc-data-download-00_data_download)
    * [Loading data into local version of the database](#loading-data-into-local-version-of-the-database-10_local_db_loading)
         * [Getting and storing some auxilliary data](#getting-and-storing-some-auxilliary-data)
         * [Measuring the field lengths and making MySQL tables](#measuring-the-field-lengths-and-making-mysql-tables)
         * [Filling and  indexing database tables](#filling-and--indexing-database-tables)
    * [Reorganizing mutation data](#reorganizing-mutation-data-20_local_db_reorganization)
         * [Removing duplicates](#removing-duplicates)
         * [Adding reliability info](#adding-reliability-info)
    * [Merging with TCGA](#merging-with-tcga)
    * [Some generic stats](#some-generic-stats)
    * [Project-specific stats](#project-specific-stats)
  
 
 
 ## Dependencies
 In addition to TCGA and ICGC databases themselves, CG relies on
 * MySQL
 * MySQLdb python module, installed with _sudo apt install python3-mysqldb_
 * gene symbols from HUGO gene nomenclature committee (see [here](https://www.genenames.org/download/custom/))
 * [Annovar](http://annovar.openbioinformatics.org/en/latest/) for location and functional annotation - int TCGA merge
 * CrossMap  (in TCGA merge) - (sudo pip3 install CrossMap, pyBigWig, pysam)
 * Optional: [line-profiler](https://github.com/rkern/line_profiler#line-profiler) for python
 
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
 
 ### ICGC data download ([00_data_download](icgc/00_data_download))
 
 Just like the tcga branch, this branch of the pipeline starts by downloading the data from
 the source, ICGC in this case. Note however that here you will need the access token. 
 Obtaining one is a lengthy process (as in several weeks to several months) which you can
 start investigating [here](https://icgc.org/daco).
 
 One you have the access token place it in the environmental variable called ICGC_TOKEN, to make
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
[06_find_max_field_length.py](icgc/10_local_db_loading/05_find_max_field_length.py) should 
give you an idea about the longest entries found.

#### Filling the icgc database tables
[07_write_mutations_tsv.py](icgc/10_local_db_loading/07_write_mutations_tsv.py) through 
[10_make_indices.py](icgc/old/10_make_indices_on_temp_tables.py).
For large tables, rather than loading them through python, 
it turns out to be faster to create tsvs and  then load them from mysql shell 
(as in [09_load_mysql.py](icgc/10_local_db_loading/09_load_mysql.py); alternative: use mysqlimport manually) 
 to read them in wholesale. These scripts take care of that part, plus some index creating on the newly loaded tables.
 Make sure to run [10_make_indices.py](icgc/old/10_make_indices_on_temp_tables.py), 
 [12_reorganize_mutations.py](icgc/20_local_db_reorganization/11_reorganize_variants.py)
 pretty much does not work without it at all. 
 All index making is slow here (see [timing.txt](icgc/timing.txt)) - run overnight. 

Some checks are thrown in here that  may inform the rest of the pipeline.
 [15_assembly_check.py](icgc/10_local_db_loading/15_assembly_check.py) confirms tha as of v27 all ICGC entried
 refer to GRCh37, and [16_donor_check.py](icgc/10_local_db_loading/16_donor_check.py) highlights the fact that
 some donor ids have no somatic mutations in ICGC. This is somewhat mysterious, because some refer to
 TCGA donors with somatic mutation data available from TCGA archive.


#### Getting and storing some auxilliary data
 
In [20_hgnc_name_resolution_table.py](icgc/10_local_db_loading/20_hgnc_name_resolution_table.py) and 
 [22_ensembl_id_table.py](icgc/10_local_db_loading/22_ensembl_id_table.py) we make and fill some tables we will use later for name resolution 
 (translating between gene and protein names used in different contexts).

The annotation across different submitters to TCGA/ICGC is not uniform In particular, for the missense mutations
sometimes it is not clear which splice they refer to. Alternatively, when the reference splice(s) is listed, it
si not clear which splice is the canonical splice. To remedy that, here we do some basic annotation of our own. For that
we will need the coding sequence of canonical transcripts.
 
The canonical transcript id is not readily available from Ensembl Mart, thus for our purposes
here you can find this info in the table called ensembl_gene2trans_stable.tsv.bz2 in the
[hacks](icgc/hacks) directory. Decompress it (bzip2 -d) and ut it someplace where
[20_hgnc_name_resolution_table.py](icgc/10_local_db_loading/20_hgnc_name_resolution_table.py) can find it.

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
New tables are created in [10_check_mut_tables_and_make_new_ones.py](icgc/20_local_db_reorganization/08_check_icgc_tables_and_make_new_ones.py).

Note that in [11_reorganize_variants.py](icgc/20_local_db_reorganization/11_reorganize_variants.py),
[13_reorganize_mutations.py](icgc/20_local_db_reorganization/13_reorganize_mutations.py),   and
[14_reorganize_locations.py](icgc/20_local_db_reorganization/14_reorganize_locations.py) 
 you can choose to run in parallel (the number of 'chunks' in main()). 

In [12_consequence_vocab.py](icgc/20_local_db_reorganization/12_consequence_vocab.py)
 we inspect the 'consequence' vocabulary employed by ICGC. There seems to
some confusion there about the location vs. the consequence of a mutation.  This info is used
in [13_reorganize_mutations.py](icgc/20_local_db_reorganization/13_reorganize_mutations.py)  to 
come up with the pathogenicity estimate, to be stored in the eponymous field in the
mutations\* tables.

[14_reorganize_locations.py](icgc/20_local_db_reorganization/14_reorganize_locations.py) script uses 
UCSC gene annotation to check chromosome addresses. The only information we are looking for here is the
possibility that the location falls within the splice region just outside of an exon. Mutations at these positions
are annotated as (possibly) pathogenic 
in [15_location_pathogenicity_to_variants.py](icgc/20_local_db_reorganization/15_set_pathogenicity_in_mutation_tables.py)

    
(Do not forget to create indices
using [10_make_indices_on_temp_tables.py](icgc/old/10_make_indices_on_temp_tables.py)) 
 
 
#### Removing duplicates
 ICGC is rife with data duplication, coming from various sources. Some seem to be bookkeeping mistakes with the
 same patient data finding its way into the dataset through various depositors; some are the results  of the re-sampling 
 of the same  tumor, while some are completely obscure, with all identifiers being identical everywhere.
 
 
  [18_cleanup_duplicate_entries.py](icgc/20_local_db_reorganization/18_cleanup_duplicate_entries.py):
  Some mutations  have identical tuple
 of identifiers (icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id). Note that this
 is after we have reorganized the database so that the mutation and location info sit in 
 different tables from the donor info. Not sure what this is about (the same sample analyzed independently multiple
 times?), but when a duplicate is found, this script chooses the entry with the greatest depth reported. 
 See the script for the full resolution strategy and make sure to run 
 [17_make_jumbo_index](icgc/20_local_db_reorganization/17_make_jumbo_index_on_new_tables.py) 
  beforehand.
 
 
 There might be further problems: See for example, mutation MU2003689, which, 
 [so the ICGC page claims](https://dcc.icgc.org/mutations/MU2003689) can be found in two distinct donors. 
 The closer  inspection of the two donors shows however that their submitter ID is the same, as is the age 
 of the 'two' women. (The tumour subtype has different description, reflecting, apparently,  the
 curator's preference.) Indeed, donors table for BRCA, at this point in the pipeline has 
 1976 distinct ICGC donor ids, and 1928 distinct submitter IDs. BRCA does turn out to be the biggest offender here,
 followed by LICA with 8 duplicated donors. It is not clear whether these duplicates refer to the same
 tumor at the same stage because even the submitter sample ids might be different
 (see [19_cleanup_multiple_donor_for_the_same_submitted_id.py](icgc/20_local_db_reorganization/19_cleanup_multiple_donor_for_the_same_submitted_id.py)). 

 Even after this cleanup we are still not done with the duplications problem - we might have the
 same donor with differing specimen and sample ids (apparently, not sure whether ICGC refers to
 biological replicates - i.e. samples taken from different sites  - as specimens, and
 to technical replicates as samples). Perhaps they might have a role when answering different
 types if questions, but here we do not want to have these results mistaken for recurring mutations, 
 thus we remove them in [22_cleanup_duplicate_specimens.py](icgc/20_local_db_reorganization/22_cleanup_duplicate_specimens.py) 
 and [23_cleanup_duplicate_samples.py](icgc/20_local_db_reorganization/23_cleanup_duplicate_samples.py), but not before checking
 which of the samples produced more reliable reads (see below).
 In this version of the pipeline we keep only the sample annotated as 'Primary tumour - solid tissue.' Out of
 these, if multiple refer to the same submitter id, we keep the ones with the largest reported number of
 somatic mutations. The investigation of the source of this duplication is again outside of our zone of interest.
 
 #### Adding reliability info
 
 We add a couple of values to each row to later make the search for meaningful entries faster.
  we are adding mutant_allele_read_count/total_read_count ratio and pathogenicity estimate (boolean)
 to simple_somatic tables. In the following script,  [21_add_realiability_annotation_to_somatic.py](icgc/20_local_db_reorganization/21_add_reliability_annotation_to_variants.py),  
 we combine these two columns into a reliability estimate: a  somatic mutation in individual patient is considered reliable if mutant_allele_read_count>=10
 and mut_to_total_read_count_ratio>=0.2. Information about the mutation in general (mutations_chromosome tables;  [18_copy_reliability_info_to_mutations.py](18_copy_reliability_info_to_mutations.py)) 
 is considered reliable if there is at leas one patient for which it was reliably established.
 
 ### Merging with TCGA ([30_tcga_merge](icgc/30_tcga_merge))
 
 The scripts [29_index_on_mutation_tables.py](icgc/20_local_db_reorganization/29_index_on_mutation_tables.py)
 through [34_tcga_specimen_hack.py](icgc/30_tcga_merge/37_tcga_specimen_hack.py)
 concern themselves with merging TCGA info created in TCGA branch with the ICGC.
 This sub-pipe runs from preparatory indexing to data input.  
 [Duplicate data removal](#removing-duplicates) should be probably be re-applied 
 (steps [19_cleanup_duplicate_donors.py](19_cleanup_duplicate_donors.py) 
 and [22_cleanup_duplicate_specimens.py](22_cleanup_duplicate_specimens.py)).
 If everything is ok, [35_database_stats.py](icgc/40_housekeeping/35_database_stats.py) should report
 no duplicates in any of the tables.
 
 ### Housekeeping ([]())
 
 ### Production ([]())
 
 
## TODO 
* disentangle from Annovar - we have all the info we need to do own annotation here

## P.S.

Do not use locks in MySQL. Locks are evil. Good luck with the rest.