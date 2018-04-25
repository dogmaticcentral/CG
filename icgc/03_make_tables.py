#!/usr/bin/python

import MySQLdb
from icgc_utils.mysql   import  *

# icgc_donor_id	 icgc_specimen_id	icgc_sample_id  submitted_sample_id
# chromosome	 chromosome_start	chromosome_end	chromosome_strand	assembly_version
# mutation_type	 reference_genome_allele
# control_genotype tumour_genotype
# mutated_from_allele	 mutated_to_allele
# consequence_type	aa_mutation	 cds_mutation
# gene_affected	transcript_affected
#########################################
def make_mutations_table(cursor, db_name, mutations_table):

	switch_to_db (cursor, db_name)

	qry = "drop table " + mutations_table
	search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "  CREATE TABLE  %s (" % mutations_table
	qry += "     id INT NOT NULL, "
	qry += "  	 icgc_mutation_id VARCHAR (20) NOT NULL, "
	qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
	qry += "     icgc_specimen_id VARCHAR (20), "
	qry += "     icgc_sample_id VARCHAR (20), "
	qry += "     submitted_sample_id VARCHAR (50), "

	qry += "	 chromosome VARCHAR (20) NOT NULL, "
	qry += "	 start_position INT  NOT NULL, "
	qry += "	 end_position INT NOT NULL, "
	qry += "	 strand VARCHAR (5) NOT NULL, "
	qry += "	 assembly VARCHAR (10) NOT NULL, "
	# for mut type use ins, del, indel and sbs (=single base substitution)
	qry += "	 mutation_type VARCHAR (10), "
	qry += "	 reference_genome_allele VARCHAR (210) NOT NULL, "
	qry += "	 control_genotype VARCHAR (430) NOT NULL, "
	qry += "	 tumor_genotype VARCHAR (430) NOT NULL, "
	qry += "	 mutated_from_allele VARCHAR (210) NOT NULL, "
	qry += "	 mutated_to_allele VARCHAR (210) NOT NULL, "

	qry += "     consequence_type  VARCHAR (100), "
	# the longest entry for aa_mutation I could find was 59
	qry += "     aa_mutation  VARCHAR (100), "
	# the longest entry for cds_mutation I could find was 11
	qry += "     cds_mutation  VARCHAR (50), "
	# it looks like these are always ENS identifiers - nice
	qry += "     gene_affected  VARCHAR (20), "
	qry += "     transcript_affected  VARCHAR (20), "


	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print qry
	print rows



#########################################
#########################################
def main():

	homedir = "/data/icgc"
	cancer_types = []
	for name in os.listdir(homedir):
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql()
	cursor = db.cursor()

	db_name =  "icgc"
	for ct in cancer_types:
		mutations_table = ct + "_simple_somatic"
		make_mutations_table(cursor, db_name, mutations_table)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()


'''
icgc fields
1   icgc_mutation_id
2   icgc_donor_id
3   project_code
4   icgc_specimen_id
5   icgc_sample_id
6   matched_icgc_sample_id
7   submitted_sample_id
8   submitted_matched_sample_id
9   chromosome
10   chromosome_start
11   chromosome_end
12   chromosome_strand
13   assembly_version
14   mutation_type
15   reference_genome_allele
16   control_genotype
17   tumour_genotype
18   expressed_allele
19   mutated_from_allele
20   mutated_to_allele
21   quality_score
22   probability
23   total_read_count
24   mutant_allele_read_count
25   verification_status
26   verification_platform
27   biological_validation_status
28   biological_validation_platform
29   consequence_type
30   aa_mutation
31   cds_mutation
32   gene_affected
33   transcript_affected
34   gene_build_version
35   platform
36   experimental_protocol
37   sequencing_strategy
38   base_calling_algorithm
39   alignment_algorithm
40   variation_calling_algorithm
41   other_analysis_algorithm
42   seq_coverage
43   raw_data_repository
44   raw_data_accession
45   initial_data_release_date
'''