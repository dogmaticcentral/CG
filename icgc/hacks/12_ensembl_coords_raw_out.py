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


#########################################
# assumes value is string
def entry_exists(cursor, db, table, column, value):
	qry = "select * from %s.%s where %s='%s' limit 1" % (db, table, column, value)
	ret =  search_db(cursor, qry)
	if not ret: return False
	if type(ret[0][0])==str and 'error' in ret[0][0].lower():
		search_db(cursor, qry, verbose=True)
		exit()
	return True


#########################################
def main():

	patch = True # we want only the coords for the trancripts we do not have already
	ref_assembly = 'hg19'

	# for this to work we need access to ensemblo hom_sapeins_core database
	db      = connect_to_mysql("/home/ivana/.tcga_conf")
	cursor  = db.cursor()
	db_name = "homo_sapiens_core_94_38"
	switch_to_db(cursor, db_name)


	# first, find database identifiers for plain old chromosomes
	chromosomes = [str(i) for i in range(1,23)] + ['X','Y']
	for chrom in chromosomes:

		outf = open("chrom_{}_coords.patch.tsv".format(chrom),"w")

		# grab the coord_system_id too so we can keep track of the assembly
		qry = "select seq_region_id, coord_system_id from seq_region where name='%s'" % chrom
		seq_region_id,coord_system_id = hard_landing_search(cursor, qry)[0]
		qry = "select version from coord_system where coord_system_id=%d" % coord_system_id
		assembly = hard_landing_search(cursor, qry)[0][0]
		print("\n==================\n", chrom, seq_region_id, assembly)

		# find protein coding genes on that gene
		qry = "select gene_id, stable_id, seq_region_strand from gene where biotype='protein_coding' "
		qry += "and seq_region_id=%d" % seq_region_id
		gene_ids = hard_landing_search(cursor,qry)
		print("number of protein coding genes:", len(gene_ids))

		total = 0
		missing = 0
		# for each gene, find transcripts, and for each transcript its start and end
		#gene_ids = [[284211, 'ENSG00000118473',  1]]
		# gene_ids = [[296055, 'ENSG00000186094',  -1]]
		for gene_id,gene_stable_id,strand in gene_ids:
			strand = "+" if strand>0 else "-"
			qry  = "select transcript_id,stable_id, seq_region_start,seq_region_end from transcript "
			qry += "where gene_id=%d " % gene_id
			transcripts = hard_landing_search(cursor,qry)
			# for each transcript
			#   get seq_start/seq_end pairs
			# for each transcript we will write positions to a file, red them in with crossmap and translate
			# this is a clunker, but not worth optimization at this point
			positions = []
			for transcript_id,transcript_stable_id,tx_start,tx_end in transcripts:
				positions.extend([tx_start,tx_end])
				if patch:
					total += 1
					if entry_exists(cursor, 'icgc', "coords_chrom_%s"%chrom, 'ens_transcript_id', transcript_stable_id): continue
					missing += 1

				#   find exons
				qry = "select exon_id, rank from exon_transcript where transcript_id=%d"%transcript_id
				exons = sorted(hard_landing_search(cursor,qry),key=lambda x: x[1],reverse=(strand=='-'))
				#   for each exon find to-from coords and frame
				exon_coords = []
				exon_coords_by_exon_id = {}
				for exon_id, rank in exons:
					qry = "select seq_region_start,seq_region_end,phase from exon where exon_id=%d" % exon_id
					exon_start, exon_end, phase = hard_landing_search(cursor,qry)[0]
					exon_coords.append([exon_start, exon_end, phase])
					positions.extend([exon_start, exon_end])
					exon_coords_by_exon_id[exon_id] = [exon_start, exon_end] # will be needed to recalculate cds_start/end

				#   find translation seq_start/seq end - they might be inverted if strand is neg
				qry = "select seq_start,start_exon_id,seq_end,end_exon_id from translation where transcript_id=%d" % transcript_id
				ret = search_db(cursor,qry)
				if not ret or len(ret)==0: continue # no translation here (shrug)
				if len(ret)>1:
					print("multiple translations (?!)")
					exit()
				seq_start,start_exon_id,seq_end,end_exon_id = ret[0]
				# seq_start is 1-based offset into the relative coordinate system of start_exon_id
				# seq_end is 1-based offset into the relative coordinate system of end_exon_id
				if strand=='+':
					cds_start = exon_coords_by_exon_id[start_exon_id][0]+seq_start-1
					cds_end   = exon_coords_by_exon_id[end_exon_id][0]+seq_end-1
				else:
					cds_start = exon_coords_by_exon_id[end_exon_id][1]-seq_end+1
					cds_end   = exon_coords_by_exon_id[start_exon_id][1]-seq_start+1

				positions.extend([cds_start, cds_end])

				# translate to GRCh37  - thr translation here is not from DNA to protein
				# but from one coord system to another
				pos_translated = translate_positions(positions, chrom, assembly, ref_assembly)
				if not pos_translated or len(pos_translated)==0: continue # translation failed; move on

				# construct ucsc-like input and compare with ucscs input that we have
				# are we counting from 1 or from 0 here?
				outfields  = [transcript_stable_id, gene_stable_id, strand]
				# 'from' gets -1 bc ucsc counts from 0; 'to' does not because ucsc uses halfopen intervals and ensembl does not
				outfields += [pos_translated[tx_start]-1, pos_translated[tx_end], pos_translated[cds_start]-1, pos_translated[cds_end]]
				exon_starts = ",".join([str(pos_translated[e[0]]-1) for e in exon_coords])
				exon_ends   = ",".join([str(pos_translated[e[1]]) for e in exon_coords])
				exon_frames = ",".join([str(e[2]) for e in exon_coords])

				outfields += [len(exons), exon_starts, exon_ends, exon_frames]
				outf.write("\t".join([str(of) for of in outfields])+"\n")

		outf.close()

		if patch: print("found in icgc:", total, " missing:", missing)

	cursor.close()
	db.close()

	return True


#########################################
if __name__ == '__main__':
	main()
