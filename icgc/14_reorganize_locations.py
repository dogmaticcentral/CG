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


import time

from config import Config
from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

variant_columns = ['icgc_mutation_id', 'chromosome','icgc_donor_id', 'icgc_specimen_id', 'icgc_sample_id',
				   'submitted_sample_id','control_genotype', 'tumor_genotype', 'total_read_count', 'mutant_allele_read_count']

#  'aa_mutation',  'consequence_type', and 'pathogenic_estimate'  will be filled separately
mutation_columns = ['icgc_mutation_id', 'start_position', 'end_position', 'assembly',
					'mutation_type', 'mutated_from_allele', 'mutated_to_allele', 'reference_genome_allele']

location_columns = ['position', 'gene_relative', 'transcript_relative']

################################################################
# stop_retained: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
consequence_vocab = ['stop_lost', 'synonymous', 'inframe_deletion', 'inframe_insertion', 'stop_gained',
					 '5_prime_UTR_premature_start_codon_gain',
					 'start_lost', 'frameshift', 'disruptive_inframe_deletion', 'stop_retained',
					 'exon_loss', 'disruptive_inframe_insertion', 'missense']

# location_vocab[1:4] is gene-relative
# location_vocab[4:] is transcript-relative
location_vocab = ['intergenic_region', 'intragenic', 'upstream', 'downstream',
				  'intron', 'exon',  'coding_sequence', '5_prime_UTR',  '3_prime_UTR',
				  'initiator_codon', 'splice_donor','splice_acceptor', 'splice_region'
				  ]

# this is set literal
pathogenic = {'missense', 'exon_loss',
				'stop_lost', 'stop_gained', 'start_lost',
				'5_prime_UTR_premature_start_codon_gain',
				'frameshift', 'inframe_deletion', 'inframe_insertion',
				'disruptive_inframe_deletion', 'disruptive_inframe_insertion',
				'inframe' }  # there is no way we can know at this level whether an inframe change is nondsisruptive
									# more likely it is than not

consequence2location = {
	'synonymous':'exon',
	'stop_retained':'exon',
	'missense':'exon',
	'exon_loss':'exon',
	'stop_lost':'exon',
	'stop_gained':'exon',
	'start_lost':'exon',
	'5_prime_UTR_premature_start_codon_gain':'5_prime_UTR',
	'frameshift':'exon',
	'inframe_deletion':'exon',
	'inframe_insertion':'exon',
	'disruptive_inframe_deletion':'exon',
	'disruptive_inframe_insertion':'exon',
	'inframe':'exon'
}

benign = {'synonymous', 'stop_retained'}

#########################################
def insert (cursor, table, columns, values):

	nonempty_values = []
	corresponding_columns = []
	for i in range(len(values)):
		if not values[i] or  values[i] == "": continue
		nonempty_values.append(values[i])
		corresponding_columns.append(columns[i])
	qry = "insert into %s (%s) " %(table, ",".join(corresponding_columns))
	qry += "values (%s) " % ",".join(nonempty_values)
	search_db(cursor, qry)


#########################################
def reorganize_locations(cursor, variants_table, columns):

	chromosomes = [str(i) for i in range(1,23)] + ["X", "Y", "MT"]

	for chromosome in chromosomes:
		location_table = "locations_chrom_{}".format(chromosome)

		positions = set()
		for position_significance in ['start', 'end']:
			qry =  "select distinct %s_position from %s  where chromosome='%s' " % (position_significance, variants_table, chromosome)
			qry += "and gene_affected is not null and gene_affected !='' "
			ret  = search_db (cursor, qry)
			if ret:  positions |= set([r[0] for r in ret])
		if len(positions)==0: continue


		for position in positions:

			# we trust that the first time around we got the annotation right
			if entry_exists(cursor, "icgc", location_table, "position", position): continue

			gene_relative       = set([])
			transcript_relative = set([])
			# this hinges on
			# qry  = "create index chrom_pos_idx on %s (chromosome, start_position)" % mutations_table
			qry =  "select * from %s where chromosome='%s' " % (variants_table, chromosome)
			qry += "and (start_position=%d or end_position=%d) " % (position, position)
			qry += "and gene_affected is not null and gene_affected !='' "

			ret = search_db (cursor,qry)
			if not ret:
				search_db (cursor,qry, verbose=True)
				exit()
			for fields in ret:
				named_field = dict(list(zip(columns,fields)))

				# this is not ready to be stored, because we need to work through the consequences
				gene   = named_field['gene_affected']
				tscrpt = named_field['transcript_affected']
				csq = named_field['consequence_type']
				if csq == location_vocab[0]: # intergenic
					pass
				elif csq in location_vocab[1:4]: # gene-relative
					gene_relative.add("{}:{}".format(gene,csq))
				elif csq in location_vocab[4:]: # transcript-relative
					gene_relative.add("{}:{}".format(gene,"intragenic"))
					transcript_relative.add("{}:{}".format(tscrpt,csq))
				elif csq in consequence2location:
					gene_relative.add("{}:{}".format(gene,"intragenic"))
					transcript_relative.add("{}:{}".format(tscrpt,consequence2location[csq]))

				elif csq == "":
					pass
				else:
					print("unrecognized consequence field:", csq)
					continue
				# gene and transcript are listed, but no relative position is specified:
				if len(gene_relative)==0 and gene and 'ENSG' in gene:
					gene_relative.add("{}:{}".format(gene,"unk"))
				if len(transcript_relative)==0 and tscrpt and 'ENST' in tscrpt:
					transcript_relative.add("{}:{}".format(tscrpt,"unk"))
			# now we are ready to store
			location_values = [str(position), quotify(";".join(gene_relative)), quotify(";".join(transcript_relative))]
			insert (cursor, location_table, location_columns, location_values)


#########################################
def reorganize(tables, other_args):

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,"icgc")
	for table in tables:

		# the tables should all have the same columns
		qry = "select column_name from information_schema.columns where table_name='%s'"%table
		columns = [field[0] for field in  search_db(cursor,qry)]
		# line by line: move id info into new table
		# for mutation and location check if the info exists; if not make new entry
		time0 = time.time()
		print("====================")
		print("reorganizing locations from", table, os.getpid())
		reorganize_locations(cursor, table, columns)
		time1 = time.time()
		print(("\t\t %s  done in %.3f mins" % (table, float(time1-time0)/60)), os.getpid())

	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	#########################
	# which temp somatic tables do we have
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%simple_somatic_temp'"
	tables = [field[0] for field in  search_db(cursor,qry)]
	cursor.close()
	db.close()

	number_of_chunks = 20  # myISAM does not deadlock
	parallelize(number_of_chunks, reorganize, tables, [])



#########################################
if __name__ == '__main__':
	main()
