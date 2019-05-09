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
from icgc_utils.utils import *
from icgc_utils.mysql import *

def resolve_id_class(cursor, id_class):
	# class can be gene or transcript
	if id_class not in ['gene', 'transcript']: return None
	id_root = 'ENS' + id_class[0].upper()

	# which ids are current
	qry = "select stable_id from %s where biotype='protein_coding'" % id_class
	ret =  search_db(cursor, qry)
	if not ret:
		search_db(cursor, qry, verbose=True)
		exit()
	current_ids = set([r[0] for r in ret])

	# get all old->new mapping pairs
	qry  = "select old_stable_id, new_stable_id from stable_id_event "
	qry += "where old_stable_id like '%s%%' and old_stable_id!=new_stable_id" % id_root
	id_pairs = search_db(cursor, qry)
	if not id_pairs:
		search_db(cursor, qry, verbose=True)
		exit()
	# rather than trying to follow the path, just lump them all ...
	# (we can do that because we are interested only on the endpoint - the current identifier
	clusters = find_clusters(id_pairs)
	print(id_class, len(id_pairs), len(clusters))

	mapped_pairs = []
	for cluster in clusters:
		idset = set(cluster)
		new_ids = current_ids.intersection(idset)
		# len(new_ids) can be 0 if the ids in question are not protein coding
		# for example, transcribed_unprocessed_pseudogene
		# in other words, the mapping exists, but we are not interested
		if len(new_ids)==0: continue
		# the intersection *can*  consist of mulitple ids (len(new_ids)>0)
		# if one old gene was mapped to 3 new (smaller) ones
		# see for example 'ENSG00000007816' --> 'ENSG00000226916', 'ENSG00000227057', 'ENSG00000236222'
		old_ids = cluster.difference(new_ids)
		# the mapping is one-to-many or many-to-one it is straightforward
		if len(old_ids)==1 or len(new_ids)==1:
			for o in old_ids:
				for n in new_ids:
					mapped_pairs.append([o,n])
		else: #otherwise we need to take a closer look
			for o in old_ids:
				for n in new_ids:
					if [o,n] in id_pairs: mapped_pairs.append([o,n])
	return mapped_pairs

#########################################
def main():
	# for this to work we need access to ensemblo hom_sapeins_core database
	db      = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor  = db.cursor()
	db_name = "homo_sapiens_core_94_38"
	switch_to_db(cursor, db_name)

	# deprecated id mapping
	mapped_pairs = resolve_id_class(cursor, 'gene')
	mapped_pairs.extend(resolve_id_class(cursor, 'transcript'))
	outf = open ("ensembl_deprecated2new_id.tsv", "w")
	idx = 0
	for mapped_pair in mapped_pairs:
		idx += 1
		outf.write("\t".join([str(idx)] + mapped_pair)+"\n")
	outf.close()

	# gene 2 transcript 2 canonical transcript
	qry  = "select t.stable_id, g.stable_id, g.canonical_transcript_id from transcript as t "
	qry += "left join gene as g on t.gene_id=g.gene_id "
	ret = search_db(cursor, qry)
	if ret:
		outf = open ("ensembl_gene2trans_stable.tsv", "w")
		for [transcript_stable_id, gene_stable_id,  canonical_transcript_id] in ret:
			ret2 =  search_db(cursor, "select stable_id from transcript where transcript_id=%s"%canonical_transcript_id)
			if not ret2: continue
			canonical_stable = ret2[0][0]
			outf.write("\t".join([transcript_stable_id, gene_stable_id, canonical_stable])+"\n")
		outf.close()
	else:
		print("No ret for", qry)



	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()
