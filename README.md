# tcga

TCGA pipeline: (pan)cancer stats. 

A series of python scripts to download select pieces of [TCGA](http://cancergenome.nih.gov/), 
store them as a local MySql database and extract various statistics.

Note: one of the biggest fudges I am doing -  when the normal allele fields are empty (systematically, in the whole file), I fill them from the other copy (as in  fill in allele 2 info from allele 1) or from the reference sequence field. Apprently this is what is usually done, I can see no evidence that anybody is keeping track of heterozygosity in the normal tissue.
