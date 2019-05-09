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

# The larget part of the reannotation seems to be innocuous - as in the
# canonical transcript has changes without changing the position on the protein seq.
# There are more serious examples, as in
#'MU113501933', 153748170, 153748170, 'GRCh37', G->A, chrom 1 + strand
# where the original annotation was 'ENST00000368661:S113N;ENST00000271857:S194N'
# and the new one is ENST00000624995:W113*
# Uniprot agrees that there is no S on position 113 nor on 194.


import time

from icgc_utils.common_queries import *
from icgc_utils.processes import *
from config import Config
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

#########################################
def genome2cds_pos(coord, mutation_position):

	#print (coord)
	# it looks like cds_start is always < cds_end
	# regardless of the strand (nice of UCSC)
	cds_start = coord['cds_start']
	cds_end   = coord['cds_end']

	exon_count  = coord['exon_count']
	exon_starts = [int(s) for s in coord['exon_starts'].strip(',').split(',')]
	exon_ends   = [int(s) for s in coord['exon_ends'].strip(',').split(',')]
	# cds start may not coincide with exon start
	length = []
	for i in range(exon_count):
		if exon_ends[i]<cds_start or cds_end<exon_starts[i]:
			length.append(0)
			continue
		start = max(exon_starts[i], cds_start)
		end   = min(cds_end, exon_ends[i])
		# the [start,end> interval is half open though
		# (so the length is end-start, rather than end-start+1
		length.append(end-start)
	# sanity checking
	if sum(length)%3>0:
		# print("CDS length not divisible by 3,", coord['ens_transcript_id'])
		# there are cases like that; in ENsmebl they are flagged as CDS 5' or 3' incomplete
		# for en example see ENST00000431696
		# there is nothing we can do except hope that it is not the case with the canonical sequence
		return None

	# UCSC is counting from 0, while other people are not
	mutation_position -= 1
	cds_position = None
	for i in range(exon_count):
		if exon_ends[i]<cds_start: continue
		if cds_end<exon_starts[i]: continue
		start = max(exon_starts[i], cds_start)
		end   = min(cds_end, exon_ends[i])
		if mutation_position<start: continue
		if mutation_position>end: break
		cds_position = sum(length[:i])+mutation_position-start
		break
	# !! coding sequence already refers to the strand that the gene is encoded in!
	if cds_position  and coord['strand']=='-': cds_position =  sum(length)-cds_position-1  # we count from 0
	#print(cds_start, cds_end, exon_starts, exon_ends)
	#f = "<<<<<<<  {}  tot_length {}  strand {}   modulo 3 {}    mut_pos-cds_start {}   cds pos {}"
	#print(f.format(coord['ens_transcript_id'], sum(length), coord['strand'],  sum(length)%3, mutation_position-cds_start, cds_position))

	return cds_position

#########################################
def find_aa_change(coding_seq, strand, cds_position, nt):

	debug = False
	# sanity:
	for label in ['from', 'to']:
		# did we end up here by mistake?
		# this is not a single nucleotide variant
		if not (len(nt[label])==1 and nt[label] in 'ATCG'): return None

	# ! coding sequence is already given  on the relevant strand
	# but nt from and nt to are not)
	pep_pos   = int(cds_position/3)
	within_codon_position = cds_position%3
	codon_seq = [coding_seq[i:i+3] for i in range(0, len(coding_seq),3)]
	if pep_pos>=len(codon_seq): return None # what's this?
	codon= {'from':codon_seq[pep_pos]}

	if debug: print("%%%%%%%%%%%", codon_seq[:12])

	# sanity: check that nt_from is what we think it is
	if strand=='-':
		complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
		nt_from = complement[nt['from']]
		nt_to   = complement[nt['to']]
	else:
		nt_from = nt['from']
		nt_to   = nt['to']

	if nt_from != codon['from'][within_codon_position]:
		if debug: print("codon mismatch", cds_position, nt_from)
		return None

	# mutate nt_from to nt_to, anc see what the aa is now
	codon['to'] = "".join([nt_to if  i==within_codon_position else codon['from'][i] for i in range(3)])

	aa = {}
	for label in ['from', 'to']:
		dna = Seq(codon[label], generic_dna)  # ! coding sequence is already given  on the relevant strand
		aa[label] = str(dna.translate())
		if debug: print(">>>", label, codon[label], aa[label])

	#if strand=='-'
	pep_pos += 1 # back to counting positions from 1
	change_string= "{}{}{}".format(aa['from'], pep_pos, aa['to'])
	return change_string

#########################################
def process_aa_change_line(cursor, chrom, line):

	#line = [... test data here ...]

	column_names = get_column_names(cursor, 'icgc', "mutations_chrom_%s"%chrom)
	named_field = dict(zip(column_names, line))

	icgc_mutation_id = named_field['icgc_mutation_id']
	start_position = named_field['start_position']
	end_position = named_field['end_position']
	aa_mutation = named_field['aa_mutation']

	# some basic sanity checking
	if start_position != end_position:
		#print("start and end position not identical for missene (?) in %s " % icgc_mutation_id)
		pass
	mutation_position = start_position

	# get (gene) location for this position
	qry = "select * from locations_chrom_%s where position=%d" % (chrom, mutation_position)
	ret = search_db(cursor, qry, verbose=False)
	if not ret:
		#print("no gene location found for %s (? it should be a missense)" % icgc_mutation_id)
		return None
	[position, gene_relative, transcript_relative] = ret[0]

	# which transcript is canonical
	canonical_transcripts = {}
	for ens_gene_id in gene_relative.split(";"):
		canonical_transcript = gene_stable_id_2_canonical_transcript_id(cursor, ens_gene_id)
		#print (ens_gene_id, get_approved_symbol(cursor, ens_gene_id), canonical_transcript)
		if canonical_transcript: canonical_transcripts[ens_gene_id] = canonical_transcript

	missing = []
	for ens_trancript_id in canonical_transcripts.values():
		if not ens_trancript_id in aa_mutation: missing.append(ens_trancript_id)

	# easy way out - nothing more to do here
	if len(missing)==0: return None


	# get transcript coordinates
	coords_table = "coords_chrom_%s" % chrom
	column_names = get_column_names(cursor, 'icgc', coords_table)
	coords = {}
	for canonical_transcript in canonical_transcripts.values():
		qry  = "select * from %s " % coords_table
		qry += "where ens_transcript_id='%s' " % canonical_transcript
		ret = search_db(cursor, qry, verbose=False)
		if not ret:
			#print("no transcript coords found for %s, %s (?)" % (icgc_mutation_id, canonical_transcript))
			continue
		coords[canonical_transcript] = dict(zip(column_names, ret[0]))
	# we cannot re-annotate because the coords are missing
	if len(coords)==0:
		#print("no coords found for", list(canonical_transcripts.values()))
		return None

	# translate genomic position to coding sequence position
	cds_positions = {}
	for enst, coord in coords.items():
		# keeo in mind here that cds_pos can be None
		# if the mutation does not hit the coding sequence
		cds_positions[enst] = genome2cds_pos(coord, mutation_position)


	# translate mutation to codon change to amino acid change
	new_annotation = []
	new_consequence = []
	new_pathogenicity = 0
	present = 0
	for enst, cds_position in cds_positions.items():
		if not cds_position:
			# we might have had problem reconstructing that gene
			# or it might be that the position does not fall within any exon
			# this is still ok because another gene might be encoded in the same region
			# and the 'missense' annotation refers to that other gene
			continue
		# find coding sequence
		qry = "select sequence from ensembl_coding_seqs where transcript_id='%s'" %  enst
		ret = search_db(cursor, qry, verbose=False)
		if not ret:
			print ("sequence not found for transcript_id='%s'" %  enst)
			continue
		coding_seq = ret[0][0]
		# note we have evaluated cds position from the left - on the "+" strand
		nt = {'from':named_field['reference_genome_allele'], 'to':named_field['mutated_to_allele']}
		change_string = find_aa_change(coding_seq, coords[enst]['strand'], cds_position, nt)
		if not change_string: continue
		annotation = "{}:{}".format(enst, change_string)
		if annotation in aa_mutation:  present += 1
		new_annotation.append(annotation)
		if change_string[0]==change_string[-1]:
			new_consequence.append('synonymous')
		elif change_string[0]=='*':
			new_consequence.append('stop_lost')
			new_pathogenicity = 1
		elif change_string[-1]=='*':
			new_consequence.append('stop_gained')
			new_pathogenicity = 1
		else:
			new_consequence.append('missense')
			new_pathogenicity = 1

	if present == len(new_annotation): return None # we already have all relevant entries

	# compare with deposited value
	annotstring = ";".join(new_annotation)
	consequence_string = ";".join(new_consequence)
	if False and len(annotstring)>0 and not annotstring in aa_mutation:
		print('\n++++++++++++++++++++++++')
		print(line)
		print("location:", position, gene_relative, transcript_relative)
		print("cds positions:", cds_positions)
		print("missing", missing)
		print("have annotation string", annotstring)
		print('++++++++++++++++++++++++\n')
	#exit()
	return annotstring, consequence_string, new_pathogenicity

#########################################
def re_annotate(chromosomes, other_args):
	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, "icgc")
	for chrom in chromosomes:
		time0 = time.time()
		print("====================")
		print("re-annotating  aa change for  ", chrom, "pid:", os.getpid())
		mutations_table = "mutations_chrom_%s" % chrom
		qry  = "select  * from %s " % mutations_table
		# note that here we trust [whoever annotated this] to have gotten at least this part right
		# (that the mutation results in aa change)
		# TODO: change this into our own independent annotation
		qry += "where aa_mutation is not null"
		ret  = search_db(cursor,qry)
		if not ret:
			print("no aa mutation entries for chrom %s (?) " % chrom)
			continue
		if 'error' in ret[0][0].lower():
			search_db(cursor,qry, verbose=True)
			exit()

		total_updates = 0
		ct = 0
		for line in ret:
			ct += 1
			if ct%1000==0: print("\t\tchrom %s  %2d%% done" % (chrom, int(float(100*ct)/len(ret))) )
			new_values =  process_aa_change_line(cursor, chrom, line)
			if not new_values: continue
			new_annotation, new_consequence, new_pathogenicity = new_values
			if len(new_annotation)==0 : continue
			# update table set aa_mutation to new_annotation
			qry  = "update %s set " % mutations_table
			qry += "aa_mutation='%s', consequence='%s', pathogenicity_estimate=%d " \
					%(new_annotation, new_consequence, new_pathogenicity)
			qry += "where icgc_mutation_id='%s' " %  line[0]
			#print(line)
			search_db(cursor, qry, verbose=False)
			total_updates += 1
			#print()


		time1 = time.time()
		print("chrom ", chrom, "done in %.3f mins, total updates %d" % (float(time1-time0)/60, total_updates))

	cursor.close()
	db.close()

	return


#########################################
#########################################
def main():

	chromosomes = [str(i) for i in range(1,13)] + ["Y"] + [str(i) for i in range(22,12,-1)] + ["X"]
	number_of_chunks = 12

	parallelize (number_of_chunks, re_annotate, chromosomes, [], round_robin=True)


	return

#########################################
if __name__ == '__main__':
	main()
