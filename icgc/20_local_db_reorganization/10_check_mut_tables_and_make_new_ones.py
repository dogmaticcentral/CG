#! /usr/bin/python3
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

# careful: mutation ids are unique for position and mutation type, not for donor!
from config import Config
from icgc_utils.mysql   import  *
from icgc_utils.icgc   import  *

#########################################
expected_chroms = set([str(i) for i in range(1,23)] + ["X", "Y", "MT"])
def sanity_checks(cursor, tables):
	# how many assemblies do I have here?
	for table in tables:
		qry = "select distinct assembly from %s "%table
		assemblies = [field[0] for field in  search_db(cursor,qry)]
		if len(assemblies) != 1  or assemblies[0] != 'GRCh37':
			print("unexpected assembly in ", table)
			print("\tdistinct assemblies:", assemblies)
			exit()
	# do I have only chromosomes or some undefined regions too?
	# do I have any non-standard chromosome annotation
	for table in tables:
		qry = "select distinct chromosome from %s "%table
		chromosomes = set([field[0] for field in  search_db(cursor,qry)])
		unexpected_chroms = chromosomes - expected_chroms
		if len(unexpected_chroms)>0:
			print("unexpected chromosomes in ", table)
			print("\t", unexpected_chroms)
			exit()

#########################################
def make_new_somatic_tables(cursor, tables):
	# variant: icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id, submitted_sample_id
	# there should be only one entry with these 4 identifiers
	# also, in the same table: total_read_count, mutant_allele_read_count

	for table in tables:
		new_table_name = table.replace("_temp","")
		make_somatic_muts_table(cursor, db_name, new_table_name)

#################
def make_mutation_tables(cursor):

	# mutation:  icgc_ mutation_id, chromosome, chrom_location, length aa_mutation cds_mutation consequence impact_estimate

	for chrom in expected_chroms:
		new_table = "mutations_chrom_{}".format(chrom)
		qry = "drop table " + new_table
		search_db(cursor, qry, verbose=True)
		qry = ""
		qry += "  CREATE TABLE  %s (" % new_table
		qry += "  	 icgc_mutation_id VARCHAR (20) NOT NULL, "
		qry += "	 start_position INT  NOT NULL, "
		qry += "	 end_position INT NOT NULL, "
		qry += "	 assembly VARCHAR (10) NOT NULL, "
		# for mut type use deletion, insertion, single and multiple
		qry += "	 mutation_type VARCHAR (10), "
		qry += "	 mutated_from_allele VARCHAR (210) NOT NULL, "
		qry += "	 mutated_to_allele VARCHAR (210) NOT NULL, "
		qry += "	 reference_genome_allele VARCHAR (210) NOT NULL, "
		# the longest entry for aa_mutation I could find was 59
		# --> yes, but I will be merging several, including ENST identifiers
		qry += "     aa_mutation  TEXT, "
		#
		qry += "     consequence TEXT, "
		qry += "     pathogenic_estimate BOOLEAN, "

		qry += "	 PRIMARY KEY (icgc_mutation_id) " # automatically indexed
		qry += ") ENGINE=MyISAM"

		rows = search_db(cursor, qry)
		print(qry)
		print(rows)

	return


#################
def make_location_tables(cursor):
	# location:  chromosome_position  gene_relative[semicolon separated list of the format ENSG:relative_location]
	# transcript_relative[semicolon separated list of the format ENST:relative_location]
	for chrom in expected_chroms:
		new_table = "locations_chrom_{}".format(chrom)
		qry = "drop table " + new_table
		search_db(cursor, qry, verbose=True)
		qry = ""
		qry += "  CREATE TABLE  %s (" % new_table
		qry += "	 position INT  NOT NULL, "
		qry += "     gene_relative TEXT, "
		qry += "     transcript_relative TEXT, "

		qry += "	 PRIMARY KEY (position) " # automatically indexed
		qry += ") ENGINE=MyISAM"

		rows = search_db(cursor, qry)
		print(qry)
		print(rows)
	return


#########################################
#########################################
def main():

	print("disabled")
	exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	switch_to_db(cursor,"icgc")
	# enable if run for the first time
	# esanity_checks(cursor, tables)

	# make new somatic mutation tables, per cancer
	make_new_somatic_tables(cursor,"icgc", tables)
	# make new mutation table, divided into chromosomes
	make_mutation_tables(cursor)
	# make new location table, divided into chromosomes
	make_location_tables(cursor)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
