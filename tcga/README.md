# tcga

TCGA pipeline: (pan)cancer stats. 

A series of python scripts to download select pieces of [TCGA](http://cancergenome.nih.gov/), 
store them as a local MySql database and extract various statistics.

Its current use is as a prep
 step for merging with icgc pipeline. The only two subdirs remaining in use are 
 [00_common_tasks](tcga/00_common_tasks) and [01_somatic_mutations](tcga/01_somatic_mutations).
 
 * [Common tasks](#common-tasks)
 * [Compiling somatic mutations](#compiling-somatic-mutations)
    * [Creating MySQL tables](#creating-mysql-tables)
    * [Reading in and cleaning up the data](#reading-in-and-cleaning-up-the-data)
    * ['Stuttering' samples](#stuttering-samples)
    * [Some basic stats](#some-basic-stats)

 
 ## Common tasks
 'Common tasks' refer to tasks needed to make a functional local subset of TCGA. The only
 non-obsolete piece remaining is [200_find_maf_files_in_GDC.py](tcga/00_common_tasks/200_find_maf_files_in_GDC.py) 
 that can be used to download somatic mutation tables from GDC - a repository of legacy TCGA data.
 
 
 ## Compiling somatic mutations
 
 ### Creating MySQL tables
 [001_drop_maf_tables]() though [002_create_maf_tables](cga/01_somatic_mutations/002_create_maf_tables.py) 
 
 ### Reading in and cleaning up the data
 [003_maf_meta](cga/01_somatic_mutations/003_maf_meta.py) through [012_drop_annotation](tcga/01_somatic_mutations/012_drop_annotation_in_remaining_conflicted.py) 
 
 
 ### 'Stuttering' samples
 Some samples in TCGA have serious problems with assembly or data interpretation. Example:
```
  broad.mit.edu_LIHC.IlluminaGA_DNASeq_automated.Level_2.1.0.0/
  An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.maf 
           273933        RPL5       Frame_Shift_Del         p.K270fs 
           273933        RPL5       Frame_Shift_Del         p.K277fs 
           273933        RPL5       Frame_Shift_Del         p.R279fs 
           273933        RPL5       Frame_Shift_Del         p.Q282fs 
```
 Such samples stand out pretty sharply and here we detect them as having two frameshift mutations within
 5 nucleotides from each other reported more than a 100 times. 
 Such samples are marked in [003_maf_meta.py](tcga/01_somatic_mutations/003_maf_meta.py).
 Later we decided to drop them in 
 [014_drop_stuttering_samples.py](tcga/01_somatic_mutations/014_drop_stuttering_samples.py) 
 
 After this point  we can move to ICGC - TCGA data will be fused into the combined dataset over there.
 
 ### Some basic stats
 ... provided by [020_db_stats.py](tcga/01_somatic_mutations/020_db_stats.py) 
 through [027_patient_freqs.py](tcga/01_somatic_mutations/027_patient_freqs.py).
 
 