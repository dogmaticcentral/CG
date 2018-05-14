#! /usr/bin/python

# careful: mutation ids are unieq for position and mutation type, not for donor!

from icgc_utils.mysql   import  *
import operator

#########################################
expected_chroms = set([str(i) for i in range(1,23)] + ["X", "Y", "MT"])
def sanity_checks(cursor, tables):
	# how many assemblies do I have here?
	for table in tables:
		qry = "select distinct assembly from %s "%table
		assemblies = [field[0] for field in  search_db(cursor,qry)]
		if len(assemblies) != 1  or assemblies[0] != 'GRCh37':
			print "unexpected assembly in ", table
			print "\tdistinct assemblies:", assemblies
			exit()
	# do I have only chromosomes or some undefined regions too?
	# do I have any non-standard chromosome annotation
	for table in tables:
		qry = "select distinct chromosome from %s "%table
		chromosomes = set([field[0] for field in  search_db(cursor,qry)])
		unexpected_chroms = chromosomes - expected_chroms
		if len(unexpected_chroms)>0:
			print "unexpected chromosomes in ", table
			print "\t", unexpected_chroms
			exit()

#########################################
def make_new_somatic_tables(cursor, tables):
	# variant: icgc_mutation_id, icgc_donor_id, icgc_specimen_id, icgc_sample_id, submitted_sample_id
	# there should be only one entry with these 4 identifiers
	# also, in the same table: total_read_count, mutant_allele_read_count

	for table in tables:
		new_table = table.replace("_temp","")
		qry = "drop table " + new_table
		search_db(cursor, qry, verbose=True)

		qry = ""
		qry += "  CREATE TABLE  %s (" % new_table
		qry += "     id INT NOT NULL AUTO_INCREMENT, "
		qry += "  	 icgc_mutation_id VARCHAR (20) NOT NULL, "
		qry += "     chromosome CHAR(2) NOT NULL,"
		qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
		qry += "     icgc_specimen_id VARCHAR (20), "
		qry += "     icgc_sample_id VARCHAR (20), "
		qry += "     submitted_sample_id VARCHAR (50), "
		qry += "	 control_genotype VARCHAR (430) NOT NULL, "
		qry += "	 tumor_genotype VARCHAR (430) NOT NULL, "
		qry += "     total_read_count INT, "
		qry += "     mutant_allele_read_count INT, "
		qry += "     mut_to_total_read_count_ratio float default 0.0,"
		qry += "     pathogenic_estimate boolean default 0,"

		qry += "	 PRIMARY KEY (id) "
		qry += ") ENGINE=MyISAM"

		rows = search_db(cursor, qry)
		print qry
		print rows

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
		# for mut type use ins, del, indel and sbs (=single base substitution)
		qry += "	 mutation_type VARCHAR (10), "
		qry += "	 mutated_from_allele VARCHAR (210) NOT NULL, "
		qry += "	 mutated_to_allele VARCHAR (210) NOT NULL, "
		qry += "	 reference_genome_allele VARCHAR (210) NOT NULL, "
		# the longest entry for aa_mutation I could find was 59
		qry += "     aa_mutation  VARCHAR (100), "
		#
		qry += "     consequence  VARCHAR (100), "
		qry += "     pathogenic_estimate BOOLEAN, "

		qry += "	 PRIMARY KEY (icgc_mutation_id) " # automatically indexed
		qry += ") ENGINE=MyISAM"

		rows = search_db(cursor, qry)
		print qry
		print rows

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
		print qry
		print rows
	return


#########################################
#########################################
def main():

	print "disabled"
	exit()

	db     = connect_to_mysql()
	cursor = db.cursor()

	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]

	switch_to_db(cursor,"icgc")
	# enable if run for the first time
	# sanity_checks(cursor, tables)

	# make new somatic mutation tables, per cancer
	make_new_somatic_tables(cursor, tables)
	# make new mutation table, divided into chromosomes
	make_mutation_tables(cursor)
	# make new location table, divided into chromosomes
	make_location_tables(cursor)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
