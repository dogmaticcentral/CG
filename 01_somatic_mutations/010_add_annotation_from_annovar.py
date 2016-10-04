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

import MySQLdb
import re
from sets import Set
from tcga_utils.mysql   import  *
from subprocess import call

#########################################
mutation_annot_pattern = re.compile('(\D+)(\-*\d+)(\D+)')
#########################################
def parse_mutation (mutation):

    match_return = re.match(mutation_annot_pattern, mutation)
    mut_from = match_return.group(1)
    mut_to   = match_return.group(3)
    mut_position = int (match_return.group(2))
    return [mut_position, mut_from, mut_to]

##################################
def output_annovar_input_file(db_name, cursor):
    qry = "select chromosome, start_position, end_position, "
    qry += "reference_allele, tumor_seq_allele1 , tumor_seq_allele2  "
    qry += "from somatic_mutations where variant_classification='missense_mutation' "
    qry += " and (aa_change is null or aa_change='')"
    rows = search_db(cursor, qry)
    outfname = "%s.avinput" % db_name
    outf = open(outfname, 'w')
    for row in rows[:3]:
        (chromosome, start_position, end_position,
            reference_allele, tumor_seq_allele1, tumor_seq_allele2 ) = row
        differing_allele = tumor_seq_allele1
        if differing_allele==reference_allele: differing_allele = tumor_seq_allele2
        print >> outf, "%s\t%d\t%d\t%s\t%s" % \
            (chromosome, start_position, end_position, reference_allele, differing_allele )
    outf.close()
    return outfname

##################################
def run_annovar(avinput, assembly, db_name):
    cmd = "/home/ivana/third/annovar/table_annovar.pl %s " % avinput
    cmd += "/home/ivana/third/annovar/humandb/ -buildver %s -out %s " % (assembly,db_name)
    cmd += " -protocol refGene  -operation g  -nastring ."
    call(cmd, shell=True)
    avoutname = "%s.%s_multianno.txt" % (db_name, assembly)
    # clean the junk
    cmd = "rm %s.refGene.variant_function " % db_name
    cmd +="%s.refGene.exonic_variant_function %s.refGene.log" % (db_name, db_name)
    call ( cmd, shell=True)
    return avoutname

##################################
def store_annotation(cursor, db_name, avoutput):
    inf = open (avoutput, "r")
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
        [aa_change_position, val1, val2] =  parse_mutation(fields[-1].replace('p.','').replace(' ', ''))
        aa_change = val1 + str(aa_change_position) + val2
        if val1==val2:
            classf = "silent"
        else:
            classf = "missense_mutation"
        qry = 'update somatic_mutations set '
        qry += 'variant_clasfication="%s",  ' % classf
        qry += 'aa_change="%s" ' % aa_change
        qry += 'cdna_change="%d" ' % cdna_change_position
        qry += 'where chromosome="%s"  ' % chrom
        qry += 'and start_position="%s" and end_position="%s"  ' % (start,end)
        search_db(cursor,qry, verbose=True)
    inf.close()
    print
    print
    return


##################################
def main():

#######
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV",   "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


    # how many cases do we have affected?
    for db_name in db_names:
        switch_to_db (cursor, db_name)
        per_db_cases = 0
        qry = "select count(1) from somatic_mutations where variant_classification='missense_mutation' "
        qry += " and (aa_change is null or aa_change='')";
        rows = search_db (cursor, qry)
        if rows and rows[0][0] != 0:
            print ">>>>>> ",  rows[0][0]
            per_db_cases = rows[0][0]

        out_of = 0
        qry = "select count(1) from somatic_mutations where variant_classification='missense_mutation' "
        rows = search_db (cursor, qry)
        if rows and rows[0][0] != 0:
            out_of = rows[0][0]

        if per_db_cases == 0: continue

        print " ** ", db_name
        print "number of missing protein annotation cases: ", per_db_cases, "out of", out_of
        qry = "select distinct(meta_info_index) from somatic_mutations where variant_classification='missense_mutation' "
        qry += " and (aa_change is null or aa_change='')"
        print "meta info for missing info cases:"
        rows = search_db (cursor, qry)
        meta_ids  = [str(row[0]) for row in rows]
        qry  = "select distinct(assembly) from mutations_meta where id in (%s)" % ",".join(meta_ids)
        rows = search_db (cursor, qry)
        if len(rows) > 1:
            print "more than one assembly - unseen at the time of writing of this script"
            exit(1)
        assembly =  rows[0][0]
        if not assembly in ["hg19", "hg18"]:
            print "unexpected assembly:", assembly
            exit(2)
        avinput = output_annovar_input_file (db_name, cursor)
        avoutput = run_annovar (avinput, assembly, db_name)
        store_annotation (cursor, db_name, avoutput)
        exit(1)






    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
