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
###################
from config import Config
from icgc_utils.common_queries import *

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

####################################################
def silent_nonsilent_count(cdna):

	mutated_codons = []
	codons = []
	# some seqs seem to miss a couple of nucleotides at the end
	# this is not really a clean solution, but we'll have to live with it for now
	# (until it is fixed in Ensembl)
	working_length = 3*int(len(cdna)/3)
	cdna = cdna[:working_length]
	standard_nucleotides = {'A','C','T','G'}
	# try to get rid of the non-standard nucleotide annotation where possible
	#  this is done in hope that there are not so many it would mess up our silent/nonsilent estimate
	ok_codon_count = 0
	bad_codon_count = 0
	for i in range(0,working_length,3):
		codon_set = set(list(cdna[i:i+3]))
		bad_codon = len(codon_set.difference(standard_nucleotides))>0
		if bad_codon:
			bad_codon_count += 1
			continue
		ok_codon_count += 1
		codon = cdna[i:i+3]
		codons.append(codon)
		for pos in range(3):
			for subs in standard_nucleotides:
				if subs == codon[pos]: continue
				mutated_codon = ''.join([subs if j==pos else codon[j] for j in range(3)])
				mutated_codons.append(mutated_codon)
	# here we bail out if there are too many bad codons
	if ok_codon_count==0 or (bad_codon_count>1 and float(bad_codon_count)/ok_codon_count>0.01): return -1, -1

	codonstring = ''.join(codons)
	mutated_codonstring = ''.join(mutated_codons)
	mutation_seq = str(Seq(mutated_codonstring, generic_dna).translate())

	aa_seq = str(Seq(codonstring, generic_dna).translate())

	if not aa_seq or len(aa_seq)==0: return -2,-2
	silent = 0
	for i in range(len(mutated_codons)):
		if int(i/9)>=len(aa_seq): return -3, -3
		if aa_seq[int(i/9)]==mutation_seq[i]: silent+=1

	return silent, len(mutated_codons)-silent

##############
def add_silent_nonsilent_columns(cursor):
	for column_name in ['silent', 'nonsilent']:
		if not column_exists(cursor, 'icgc', 'ensembl_coding_seqs', column_name):
			qry ="alter table icgc.ensembl_coding_seqs add column %s int default -1" % column_name
			error_intolerant_search(cursor,qry)
	return

####################################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	add_silent_nonsilent_columns(cursor)
	switch_to_db(cursor, 'icgc')
	errmsg = "SEQUENCEUNAVAILABLE"
	for trs_id, cdna in hard_landing_search(cursor, "select transcript_id, sequence from ensembl_coding_seqs"):
		if cdna[:len(errmsg)]==errmsg: continue
		silent, nonsilent = silent_nonsilent_count(cdna)
		if silent<0:
			print (trs_id, silent, nonsilent)
			exit()
		qry  = "update ensembl_coding_seqs set silent=%d, nonsilent=%d " % (silent, nonsilent)
		qry += "where transcript_id='%s'" % trs_id
		error_intolerant_search(cursor,qry)
	cursor.close()
	db.close()

	return


#########################################
if __name__ == '__main__':
	main()

