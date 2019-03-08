# CG

_(writeup in progress)_

CG is a set of scripts for extracting information about somatic 
mutations from 
[TCGA](https://portal.gdc.cancer.gov/)
and 
[ICGC](https://dcc.icgc.org/)
databases. The ICGC branch also contains tools necessary to
merge the info from the two sources. The scripts add up to a loosely connected pipeline 
(or, rather, two pipelines).  They re-organize the data and store it in the local MySQL 
database.  While the  back end of each pipeline is rather generic, the front line is
geared toward answering particular questions for which they were originally written.

CG is not an out-of-the box solution. Rather, it is a starter kit, in case you would like to
do some cancer genomics data analysis on your own.


Why bother with the local copy of the data? I do not know a general answer to that question.
You should check the homepage  of both databases - maybe the information you
 are looking for is already  available.
Here the original motif was to study the co-occurrence of mutations in different genes. 
These types of questions are not readily answerable using data portals - they  tend to
agglomerate data on per-gene basis, in order to protect the privacy of sample donors. 


<!-- this is a comment -->
<!-- making TOC: https://github.com/ekalinin/github-markdown-toc -->
<!-- once installed, use with gh-md-toc README.md    -->
 
 ## Dependencies
 
 
 ## ICGC
 
 ### config file
 You can set some recurring constants - such as data directories or mysql conf file(s) - 
 in the [config.py](icgc/config.py) file.
 
 ### 00_data_download
 
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
 
 ### 01_hgnc_name_resolution_table and 02_ensembl_id_table
 
 Here we make and fill some tables we will use later for name resolution 
 (translating between gene and protein names used in different contexts). Make sure you have 
 the mysql conf file  and set its path in these two scripts, or arrange some other way to
 access the local database. The last I checked, python's MySQLdb package did not work with
 the encripted cnf files, so the only alternative is using 
 [client field in generic mysql option file](https://dev.mysql.com/doc/refman/8.0/en/option-files.html),
 like this:
 
`[client]`   
`user = blah`  
`host = localhost`  
`password = "somepasswd"`

The canonical transcript id is not readily available from Ensembl Mart, thus for our purposes
here you can find this info in the table called ensembl_gene2trans_stable.tsv.bz2 in the
[hacks](icgc/hacks) directory. Put it someplace where
[02_ensembl_id.py](icgc/02_ensembl_id.py) can find it.

### 03_find_max_field_length and 04_make_tables
Make sure that the fields in the mysql tables are big enough for each entry and create mysql tables.

### 05_write_mutations through 08_make_indices
For large tables, rather than loading them through python, 
it turns out to be faster to create tsvs and 
then load them from mysql shell (as in [07_load_mysql.py](icgc/07_load_mysql.py); alternative: use mysqlimport manually) 
 to read them in wholesale. These scripts take care of that part , plus some index creating on the newly loaded tables.
 Make sure to run [08_make_indices.py](icgc/08_make_indices_on_temp_tables.py) - [12_reorganize_mutations.py](icgc/12_reorganize_mutations.py)
 pretty much does not wokr without it at all. 

### 10_check_mut_etc through 17_copy_reliabilty_etc
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

Note that in [12_reorganize_mutations.py](icgc/12_reorganize_mutations.py) you can choose to
run in parallel (the number of 'chunks' in main()). It still takes a while - as in, 
leave-to-run-overnight while. This is probably the weakest part of the whole pipeline, but is unclear
whether it is worth the optimization effort. (Do not forget to create indices
 using [08_make_indices.py](icgc/08_make_indices_on_temp_tables.py))
 
 Even after we moved mutation and location ino to separate tables; 
 some entries in ICGC appear to be duplicates - the same mutation, donor, specimen, and sample id. 
 Not sure what this is about, because ICGC is not actually terribly well documented. 
 We are removing duplicates in [13_cleanup_duplicate_entries.py](icgc/14_cleanup_duplicate_entries.py).
 Even after this cleanup there might be further problems:
 
 See for example, mutation MU2003689, which, so 
 [the ICGC page claims][https://dcc.icgc.org/mutations/MU2003689] can be found in two distinct donors. 
 The closer  inspection of the two donors shows however that their submitter ID is the same, as is the age 
 of the 'two' women. (The tumour subtype has different description, reflecting, apparently,  the
 curator's preference.) Indeed, donors table for BRCA, at this point in the pipeline has 
 1976 distinct ICGC donor ids, and 1928 distinct submitter IDs. BRCA does turn out to be the biggest offender here,
 followed by LICA with 8 duplicated donors. It is not clear whether these duplicates refer to the same
 tumor at the same stage because even the submitter sample ids might be different
 (see [14_cleanup_duplicate_donors.py](icgc/15_cleanup_duplicate_donors.py)). Not sure if this is worth pursuing
 further, except for being very cautious abut making claims  about recurrent mutations, in BRCA in particular.
 
 If [12_reorganize_mutations.py](icgc/12_reorganize_mutations.py) is the weakest link in the pipeline, 
 [14_cleanup_duplicate_entries.py](icgc/14_cleanup_duplicate_entries.py) is the most likely to cover-up for a problem, 
 possibly originating in ICGC itself. Some data sets seem to have a huge number of duplicates - entries with identical tuple
 of identifiers (icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id). Note that this
 is after we have reorganized the database so that the mutation and location info sit in 
 different tables from the door info. Not sure what this is about (the same sample analyzed independently multiple
 times?), but when found, this script chooses the entry with the highest coverage if possible. See the script for the full
 resolution strategy and make sure to run [13_make_jumbo_index](icgc/13_make_jumbo_index_on_simple_somatic_tables.py) -
 [14_cleanup_duplicate_entries.py](icgc/14_cleanup_duplicate_entries.py) is useless without it.
 