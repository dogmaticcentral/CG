#!/usr/bin/python3

from icgc_utils.mysql   import  *
from config import Config

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
	qry += "     total_read_count INT, "
	qry += "     mutant_allele_read_count INT, "


	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)

#########################################
# icgc_donor_id,submitted_donor_id,donor_sex,donor_diagnosis_icd10
def make_donors_table(cursor, db_name, donor_table):

	switch_to_db (cursor, db_name)

	qry = "drop table " + donor_table
	search_db(cursor, qry, verbose=True)


	qry = ""
	qry += "  CREATE TABLE  %s (" % donor_table
	qry += "     id INT NOT NULL, "
	qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
	qry += "     submitted_donor_id VARCHAR (50) NOT NULL, "
	qry += "	 donor_sex VARCHAR (10), "
	qry += "	 donor_diagnosis_icd10  VARCHAR (20), "

	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)


#########################################
# icgc_specimen_id, icgc_donor_id, specimen_type, tumour_histological_type
def make_specimen_table(cursor, db_name, specimen_table):

	switch_to_db (cursor, db_name)

	qry = "drop table " + specimen_table
	search_db(cursor, qry, verbose=True)


	qry = ""
	qry += "  CREATE TABLE  %s (" % specimen_table
	qry += "     id INT NOT NULL, "
	qry += "  	 icgc_specimen_id VARCHAR (20) NOT NULL, "
	qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
	qry += "     specimen_type VARCHAR (150), "
	qry += "     tumour_histological_type VARCHAR (150), "
	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)


#########################################
#########################################
def main():

	print("disabled bcs it drops and re-creates tables - comment out to run")
	return

	homedir = Config.data_home_local
	cancer_types = []
	for name in os.listdir(homedir): # this is not walking recursicely - it's just this dir
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name =  "icgc"
	for ct in cancer_types:
		# note 'temp' here - we'll reroganize when we get to TCGA
		mutations_table = ct + "_simple_somatic_temp"
		make_mutations_table(cursor, db_name, mutations_table)
		donors_table = ct + "_donor"
		make_donors_table(cursor, db_name, donors_table)
		specimen_table = ct + "_specimen"
		make_specimen_table(cursor, db_name, specimen_table)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()


'''
somatic fields
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

donor fields
1 icgc_donor_id
2 project_code
3 study_donor_involved_in
4 submitted_donor_id
5 donor_sex
6 donor_vital_status
7 disease_status_last_followup
8 donor_relapse_type
9 donor_age_at_diagnosis
10 donor_age_at_enrollment
11 donor_age_at_last_followup
12 donor_relapse_interval
13 donor_diagnosis_icd10
14 donor_tumour_staging_system_at_diagnosis
15 donor_tumour_stage_at_diagnosis
16 donor_tumour_stage_at_diagnosis_supplemental
17 donor_survival_time
18 donor_interval_of_last_followup
19 prior_malignancy
20 cancer_type_prior_malignancy
21 cancer_history_first_degree_relative

specimen fields
1 icgc_specimen_id
2 project_code
3 study_specimen_involved_in
4 submitted_specimen_id
5 icgc_donor_id
6 submitted_donor_id
7 specimen_type
8 specimen_type_other
9 specimen_interval
10 specimen_donor_treatment_type
11 specimen_donor_treatment_type_other
12 specimen_processing
13 specimen_processing_other
14 specimen_storage
15 specimen_storage_other
16 tumour_confirmed
17 specimen_biobank
18 specimen_biobank_id
19 specimen_available
20 tumour_histological_type
21 tumour_grading_system
22 tumour_grade
23 tumour_grade_supplemental
24 tumour_stage_system
25 tumour_stage
26 tumour_stage_supplemental
27 digital_image_of_stained_section
28 percentage_cellularity
29 level_of_cellularity


'''

