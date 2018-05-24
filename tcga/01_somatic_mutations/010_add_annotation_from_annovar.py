#!/usr/bin/python

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


import re

from tcga_utils.mysql import *
from tcga_utils.processes import *
from subprocess import call

#########################################
mutation_annot_pattern = re.compile('(\D+)(\-*\d+)(\D+)')
#########################################
def parse_mutation (mutation):

	match_return = re.match(mutation_annot_pattern, mutation)
	if not match_return: return [None, None, None]
	mut_from = match_return.group(1)
	mut_to   = match_return.group(3)
	mut_position = int (match_return.group(2))
	return [mut_position, mut_from, mut_to]

##################################
def output_annovar_input_file(cursor, db_name, table_name, assemblies):
	meta_table_name = table_name.split("_")[0] + "_mutations_meta"
	qry  = "select s.chromosome, s.start_position, s.end_position, "
	qry += "s.reference_allele, s.tumor_seq_allele1, s.tumor_seq_allele2, m.assembly  "
	qry += "from %s s, %s m " % (table_name, meta_table_name)
	qry += "where variant_classification='missense_mutation' "
	qry += "and (aa_change is null or aa_change='') "
	qry += "and s.meta_info_id=m.id"

	rows = search_db(cursor, qry)
	outf = {}
	outfname = {}
	for assembly in assemblies:
		outfname[assembly] = "%s.%s.avinput" % (table_name,assembly)
		outf[assembly] = open(outfname[assembly], 'w')
	for row in rows:
		(chromosome, start_position, end_position,reference_allele,
			tumor_seq_allele1, tumor_seq_allele2, assembly ) = row
		differing_allele = tumor_seq_allele1
		if differing_allele==reference_allele: differing_allele = tumor_seq_allele2
		print >> outf[assembly], "%s\t%d\t%d\t%s\t%s" % \
			(chromosome, start_position, end_position, reference_allele, differing_allele)
	for assembly in assemblies:
		outf[assembly].close()
	return outfname

##################################
def run_annovar(avinput, table_name):
	avoutname = {}
	for assembly, avinfile in avinput.iteritems():
		cmd  = "/home/ivana/third/annovar/table_annovar.pl %s " % avinfile
		cmd += "/home/ivana/third/annovar/humandb/ -buildver %s -out %s " % (assembly, table_name)
		cmd += " -protocol refGene  -operation g  -nastring ."
		call(cmd, shell=True)
		avoutname[assembly] = "%s.%s_multianno.txt" % (table_name, assembly)
		# clean the junk
		cmd = "rm %s.refGene.variant_function " % table_name
		cmd +="%s.refGene.exonic_variant_function %s.refGene.log" % (table_name, table_name)
		call ( cmd, shell=True)
	return avoutname

##################################
def store_annotation(cursor, table_name, avoutput):
	print ">>>>>> storing annotation for", table_name
	for assembly, avfile in avoutput.iteritems():
		inf = open (avfile, "r")
		for line in inf:
			if line[:3]=="Chr": continue
			fields = line.rstrip().split('\t')
			[chrom, start, end] = fields[:3]
			fields = fields[-1].split(',')[0].split(':')
			# in some cases annovar believes this is not exonic change at all
			# I am not sure wha to do in such case, and I am so sick and tired
			# of this godawful data set
			if len(fields)<2: continue
			[cdna_change_position, val1, val2] =  parse_mutation(fields[-2])
			if not cdna_change_position: continue
			[aa_change_position, val1, val2] =  parse_mutation(fields[-1].replace('p.','').replace(' ', ''))
			if not aa_change_position: continue
			aa_change = val1 + str(aa_change_position) + val2
			if val1==val2:
				classf = "silent"
			else:
				classf = "missense_mutation"
			qry = 'update %s set ' % table_name
			qry += 'variant_classification="%s",  ' % classf
			qry += 'aa_change="%s", ' % aa_change
			qry += 'cdna_change="%d" ' % cdna_change_position
			qry += 'where chromosome="%s"  ' % chrom
			qry += 'and start_position="%s" and end_position="%s"  ' % (start,end)
			search_db(cursor,qry)
		inf.close()
		print
		print
	return


##################################
##################################
def annotate_tables(tables,other_args):

	home = os.getcwd()

	db     = connect_to_mysql()
	cursor = db.cursor()

	db_name = "tcga"
	switch_to_db(cursor, db_name)

	# how many cases do we have affected?
	for somatic_mutations_table in tables:
		print
		print "========================="
		print somatic_mutations_table, os.getpid()

		tumor_short = somatic_mutations_table.split("_")[0]

		# see if there is annotation missing
		qry = "select count(1) from %s where variant_classification='missense_mutation' " % somatic_mutations_table
		qry += " and (aa_change is null or aa_change='')";
		rows = search_db (cursor, qry)
		if rows and rows[0][0] != 0:
			per_db_cases = rows[0][0]
		else:
			print ">>>>>> annotation complete"
			continue # no annotation missing - move on

		qry = "select count(1) from %s where variant_classification='missense_mutation' " % somatic_mutations_table
		out_of = search_db(cursor, qry)[0][0]
		print ">>>>>> number of missing protein annotation cases: ", per_db_cases, "out of", out_of

		continue
		# what is the assembly used here
		meta_table = tumor_short +"_mutations_meta"
		qry = "select distinct(meta_info_id) from %s  " % somatic_mutations_table
		qry += "where variant_classification='missense_mutation' "
		qry += "and (aa_change is null or aa_change='')"
		rows = search_db (cursor, qry)
		meta_ids  = [str(row[0]) for row in rows]

		qry  = "select distinct(assembly) from %s where id in (%s)" % (meta_table, ",".join(meta_ids))
		rows = search_db (cursor, qry)
		assemblies = [row[0] for row in rows]
		if len(rows) > 1:
			print ">>>>>> note: more than one assembly"
		print ">>>>>> "+ "  ".join(assemblies)

		if set(assemblies).difference({"hg19", "hg18"}):
			print "unexpected assembly/ies:", set(assemblies).difference({"hg19", "hg18"})
			exit(2)

		###### make a workdir and move there
		workdir  = tumor_short + "_annovar"
		workpath = "{}/annovar/{}".format(home,workdir)
		if not os.path.exists(workpath): os.makedirs(workpath)
		os.chdir(workpath)
		avinput  = output_annovar_input_file (cursor, db_name, somatic_mutations_table, assemblies)
		for assembly, avinfile  in avinput.iteritems():
			print "{}/{}".format(workpath,avinfile)
		avoutput = run_annovar(avinput, somatic_mutations_table)
		store_annotation (cursor, somatic_mutations_table, avoutput)

	cursor.close()
	db.close()


##################################
def main():
	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='tcga' and table_name like '%_somatic_mutations'"
	tables = [field[0] for field in search_db(cursor,qry)]

	cursor.close()
	db.close()

	number_of_chunks = 1  # myISAM does not deadlock
	parallelize(number_of_chunks, annotate_tables, tables, [])


#########################################
if __name__ == '__main__':
	main()
