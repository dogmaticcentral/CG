#! /usr/bin/python
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

from icgc_utils.common_queries import *
import random

def main():

	tumor_short = "COCA"

	db     = connect_to_mysql()
	cursor = db.cursor()


	switch_to_db(cursor,'icgc')
	snvs = {}
	silent = {}
	other  = {}
	silent_ratio = {}
	qry = "select icgc_mutation_id, chromosome from  %s_simple_somatic " % tumor_short
	qry += "where reliability_estimate=1"
	for mut_id, chrom in search_db(cursor, qry):
		#silent
		qry  = "select  g.gene_symbol, m.consequence from mutation2gene g, mutations_chrom_%s  m " % chrom
		qry += "where  m.icgc_mutation_id=g.icgc_mutation_id and m.icgc_mutation_id='%s' " % mut_id
		ret = search_db(cursor,qry)
		if not ret: continue # not sure how that could happen

		for gene, consq in ret:
			if not consq in ['synonymous', 'missense', 'stop_gained']: continue
			#print "*{}*  *{}*".format(gene, consq)
			if not silent.has_key(gene):
				silent[gene] = 0
				other[gene] = 0
			if consq=='synonymous':
				silent[gene] += 1
			else:
				other[gene] += 1

	for gene, silent_count in silent.iteritems():
		other_count = other[gene]
		if silent_count==0 and other_count==0: continue # paranoia
		if (silent_count+other_count)<3: continue
		silent_ratio[gene] = float(silent_count)/(silent_count+other_count)
		snvs[gene] = (silent_count+other_count)

	gene_ratio_list =[(gene,ratio) for gene,ratio in sorted(silent_ratio.iteritems(), key=lambda (k,v): v)]

	for gene,ratio in gene_ratio_list[:100]:
		print "%15s  %5.2f  %4d " % (gene, ratio, snvs[gene])

	outf = open ("silent_ratio.%s.tsv"%tumor_short,"w")
	for gene,ratio in gene_ratio_list:
		outf.write("\t".join([gene, "%.2f"%ratio, "%d"%snvs[gene]])+"\n")
	outf.close()

	cursor.close()
	db.close()



#########################################
if __name__ == '__main__':
	main()
