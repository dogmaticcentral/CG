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

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config
from random import sample
from math import sqrt

####################################################
def find_symbol_in_ensembl(cursor, approved_symbol, verbose=False):
	# note that some of these might not be protein coding,
	# but this is not the issue we are resolving here
	switch_to_db(cursor, 'homo_sapiens_core_94_38')
	qry  = "select x.display_label, x.description, g.stable_id "
	qry += "from xref x, object_xref o, gene g   "
	qry += "where  x.dbprimary_acc like 'HGNC%' "
	qry += "and x.display_label='%s'  " % approved_symbol
	qry += "and o.xref_id=x.xref_id and o.ensembl_object_type='Gene' "
	qry += "and o.ensembl_id=g.gene_id"
	ret = error_intolerant_search(cursor,qry)
	switch_to_db(cursor, 'icgc')
	if not ret or len(ret)==0: return None

	if verbose and len(ret)>1:
		print("Multiple ens ids found for:", approved_symbol)
		print(ret)
	ensid = ret[0][2]
	if ensid[:5] != "ENSG0":
		if verbose: print("error retrieving ensid for {}, found {}".format(approved_symbol, ensid))
		return None
	return ensid




####################################################
def main():


	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	switch_to_db(cursor, 'icgc')

	qry= "select id, approved_symbol, ensembl_gene_id, ensembl_gene_id_by_hgnc from hgnc"
	for  id, approved_symbol, ensembl_gene_id , ensembl_gene_id_by_hgnc in hard_landing_search(cursor,qry):

		# nothing we *need*  to do if ensid present and can be trnslated to ENST
		if ensembl_gene_id=='':
			# if neither ensembl column is filled try to it in ensembl homosapiens
			if ensembl_gene_id_by_hgnc=='':
				ensembl_gene_id = find_symbol_in_ensembl(cursor, approved_symbol)
			else:
				ensembl_gene_id = ensembl_gene_id_by_hgnc

			if ensembl_gene_id and len(ensembl_gene_id)>0 and ensembl_gene_id[:5] == "ENSG0":
				qry = "update hgnc set ensembl_gene_id='%s' where id=%d" %(ensembl_gene_id, id)
				error_intolerant_search(cursor,qry)
			else: # we cannot find the ensembl id for this one and we move on
				continue

		# if we are here, we have a non-empty candidate for ensembl_gene id
		qry  = "select  distinct(canonical_transcript) from icgc.ensembl_ids where  gene ='%s' " % ensembl_gene_id
		canonical_transcript = error_intolerant_search(cursor,qry)
		if canonical_transcript: continue

		# we have gene id but no transcript id
		new_id = attempt_resolve_deprecated(cursor, ensembl_gene_id, verbose=False)

		if not new_id: continue # cannot be resolved

		qry = "update hgnc set ensembl_gene_id='%s' where id=%d" %(new_id, id)
		error_intolerant_search(cursor,qry)


		qry  = "select  distinct(canonical_transcript) from icgc.ensembl_ids where  gene ='%s' " % new_id
		canonical_transcript = error_intolerant_search(cursor,qry)
		print(ensembl_gene_id, new_id, canonical_transcript)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
