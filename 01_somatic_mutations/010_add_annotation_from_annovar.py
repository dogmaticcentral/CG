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
from sets import Set
from   tcga_utils.mysql   import  *
from subprocess import call

##################################
def output_annovar_input_file(db_name, cursor):
    qry = "select chromosome, start_position, end_position, reference_allele, tumor_seq_allele1 "
    qry += "from somatic_mutations where variant_classification='missense_mutation' "
    qry += " and (aa_change is null or aa_change='')"
    rows = search_db(cursor, qry)
    outfname = "%s.avinput" % db_name
    outf = open(outfname, 'w')
    for row in rows[:10]:
        print >> outf, "%s\t%d\t%d\t%s\t%s" % row
    outf.close()
    return outfname

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
        if assembly != "hg19":
            print "unexpected assembly:", assembly
            exit(2)
        avinput = output_annovar_input_file (db_name, cursor)
        cmd = "/home/ivana/third/annovar/table_annovar.pl %s " % avinput
        cmd += "/home/ivana/third/annovar/humandb/ -buildver hg19 -out %s " % db_name
        cmd += " -protocol refGene  -operation g  -nastring ."
        print cmd
        call(cmd, shell=True)
        exit(1)






    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
