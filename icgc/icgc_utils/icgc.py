
#
# This source code is part of icgc, an ICGC processing pipeline.
# 
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
from icgc_utils.mysql import switch_to_db, search_db, check_table_exists

#########################################
def make_mutations_table(cursor, db_name, mutations_table):

	switch_to_db (cursor, db_name)

	if check_table_exists(cursor, db_name, mutations_table):
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
# submitted donor ids can be very long especially for TCGA donors
def make_donors_table(cursor, db_name, donor_table):

	switch_to_db (cursor, db_name)

	if check_table_exists(cursor, db_name, donor_table):
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

	if check_table_exists(cursor, db_name, specimen_table):
		qry = "drop table " + specimen_table
		search_db(cursor, qry, verbose=True)

	qry = ""
	qry += "CREATE TABLE  %s (" % specimen_table
	qry += "     id INT NOT NULL, "
	qry += "  	 icgc_specimen_id VARCHAR (50) NOT NULL, "
	qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
	qry += "     specimen_type VARCHAR (150), "
	qry += "     tumour_histological_type VARCHAR (150), "
	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)

