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
 might not be clear which cancer the depositors refer to. Feel free to change the grouping
 defined in [06_group_cancers.py](icgc/00_data_download/06_group_cancers.py), or skip it altogether.
 
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


