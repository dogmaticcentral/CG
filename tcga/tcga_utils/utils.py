#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#
from mysql import *

###############################################################################################
def is_informative(string):
	if not string or len(string)==0: return False
	string = string.replace(" ", "")
	if  len(string)==0: return False
	non_info  = ['missing', '', '.', '-', '---', 'untested']
	return not string in non_info

###############################################################################################
def is_useful(fields, header):
	field_is_useful = fields != None and fields.has_key(header) and fields[header] != None
	if field_is_useful and type(fields[header]) is str: # additional check
		field_is_useful = field_is_useful and is_informative(fields[header])

	return field_is_useful

import traceback
################################################################################################
def make_named_fields (header_fields, fields, expected_fields = None):
	# type: (object, object, object) -> object

	named_fields = {}

	if len(header_fields) != len(fields) :
		print "##################################"
		print "fields length mismatch (?)" # it should have been solved by this point
		print len(header_fields), len(fields)
		for line in traceback.format_stack():
			print(line.strip())
		exit(1) # header field mismatch

	for i in range(len(header_fields)):
		header = header_fields[i]
		if expected_fields and not header in expected_fields: continue
		field = fields[i]
		named_fields[header] = field

	return named_fields

################################################################################################
def process_header_line(maffile):

	inff = open(maffile, "r")
	headerline = ""
	for line in inff:
		if line.isspace(): continue
		if line[0]=='#': continue
		headerline = line.rstrip()
		break # the first line that is not the comment should be the header
	inff.close()

	if len(headerline)==0: return []
	header_fields = [x.lower() for x in headerline.split('\t')]
	# translate the header nomenclature to what we expect to see:
	for i in range(len(header_fields)):
		if header_fields[i] == 'chrom':
			header_fields[i] ='chromosome'
		elif header_fields[i] in ['amino_acid_change_wu',
							  'aachange', 'amino_acid_change',
							  'protein_change', 'hgvsp_short']:
			header_fields[i] = 'aa_change'
		elif header_fields[i] == 'tumor_sample_id':
			header_fields[i] ='tumor_sample_barcode'
		elif header_fields[i] == 'normal_sample_id':
			header_fields[i] ='normal_sample_barcode'
		elif header_fields[i] in ['cdna_change', 'chromchange', 'c_position_wu', 'c_position', 'hgvsc']:
			header_fields[i] = 'cdna_change'

	return header_fields

################################################################################################
# Required_fields are the absolute minimum we need to 
# reconstruct the mutation - that they are missing  should not happen at all
def get_required_fields ():
	return ['start_position',  'reference_allele', 'tumor_seq_allele1', 'tumor_seq_allele2',
			'tumor_sample_barcode', 'match_norm_seq_allele1',  'match_norm_seq_allele2']
 
################################################################################################
# The fields we have predeclared in the database are what we expect to see in the file
# The rest of the code should handle the missing fields somehow.
def get_expected_fields(cursor, db_name, table):

	if not check_table_exists(cursor,db_name,table):
		print "table %s not found in %s" % (table, db_name)
		exit()
	qry  = "describe %s " % table
	rows = search_db (cursor, qry)
	if not rows:
		print "%s not found in  %s" % (table, db_name)
		exit(1)
	# the first header field is 'id' - this is entry id in the table
	expected_fields = [row[0].lower() for row in rows[1:]]
	return expected_fields

################################################################################################
def hugo2ensembl (cursor, hugo_id):

	ensembl_id = ""
	switch_to_db(cursor, 'baseline')
	qry = "select ensembl_gene_id from hgnc_id_translation where approved_symbol='%s'" %  hugo_id
	rows  = search_db (cursor, qry)
	if rows and rows[0]:
		return rows[0][0]

	for alternative_column in ['previous_symbols','synonyms']:
		qry  = "select ensembl_gene_id, %s " %  alternative_column
		qry += "from hgnc_id_translation where ensembl_gene_id is not null and  ensembl_gene_id != ''"
		qry += "and %s like '%%%s%%'  " %  (alternative_column, hugo_id)
		rows  = search_db (cursor, qry)

		if rows:
			for row in rows:
				possible_ensembl  = row[0]
				alter_symbols  = row[1].replace (' ', '')
				fields = alter_symbols.split( ',')
				for field in fields:
					if field == hugo_id:
						return possible_ensembl
	return ""


#########################################
def read_cancer_names ():
	full_name= {}
	inf = open ("/home/ivana/pypeworks/cg/cancer_names.tsv", "r")
	for line in inf:
		line   = line.rstrip()
		field = line.split ("\t")
		if field[0] == 'READ':
			field[0] = 'REA'
		full_name[field[0]] = field[1]
	inf.close()

	return full_name

