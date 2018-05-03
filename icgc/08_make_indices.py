#!/usr/bin/python

import MySQLdb
from icgc_utils.mysql   import  *

#########################################
#########################################
def main():

	homedir = "/data/icgc"
	cancer_types = []
	for name in os.listdir(homedir):
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql()
	cursor = db.cursor()

	switch_to_db(cursor,"icgc")
	for ct in cancer_types:
		mutations_table = ct + "_simple_somatic_temp"
		qry  = "create index chrom_pos_idx on %s (chromosome, start_position)" % mutations_table
		#qry  = "create index donor_mutation_idx on %s (icgc_donor_id)" % mutations_table
		search_db(cursor,qry,verbose=True)


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

