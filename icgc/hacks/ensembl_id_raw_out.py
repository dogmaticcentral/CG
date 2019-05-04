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

from icgc_utils.mysql import *

#########################################
def main():

	db      = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor  = db.cursor()
	db_name = "homo_sapiens_core_94_38"
	switch_to_db(cursor, db_name)

	# deprecated id mapping
	outf = open ("ensembl_deprecated2new_id.tsv", "w")

	qry  = "select s.old_stable_id, s.new_stable_id from stable_id_event as s, gene as g "
	qry += "where s.old_stable_id like 'ENSG%' "
	qry += "and  s.new_stable_id=g.stable_id and s.old_stable_id!=s.new_stable_id"
	ret = search_db(cursor, qry)
	id = 0
	if ret:
		for line in ret:
			id += 1
			outf.write("\t".join([str(id)]+line)+"\n")
	else:
		print("No ret for", qry)
	qry  = "select s.old_stable_id, s.new_stable_id from stable_id_event as s, transcript as t "
	qry += "where s.old_stable_id like 'ENST%' "
	qry += "and  s.new_stable_id=t.stable_id and s.old_stable_id!=s.new_stable_id"
	ret = search_db(cursor, qry)
	if ret:
		for line in ret:
			id += 1
			outf.write("\t".join([str(id)]+line)+"\n")
	else:
		print("No ret for", qry)
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
