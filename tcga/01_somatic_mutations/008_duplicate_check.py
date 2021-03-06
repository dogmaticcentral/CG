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
import commands


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    sample_type = "metastatic"

    if sample_type == "primary":
        table_name = 'somatic_mutations'
    elif sample_type == "metastatic":
        table_name = 'metastatic_mutations'
    else:
        print "I don't know how to hadndle ", sample_type, " sample types"
        exit(1) # unknown sample type

    db_names = ["ACC", "BLCA", "BRCA", "CESC", "CHOL",  "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",  "MESO", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


    # how many cases do we have affected?
    total_cases = 0
    total_cases_conflict = 0
    for db_name in db_names:
        print " ** ", db_name
        switch_to_db (cursor, db_name)

        qry  = "select distinct  sample_barcode_short from %s " % table_name
        rows = search_db(cursor, qry)
        per_db_cases = 0
        per_db_cases_conflict = 0
        if not rows: continue
        for  row in rows:
            short_barcode = row[0]
            # this is wrong - not taking different chromosomes into account
            qry =  'select chromosome, start_position, count(*) from %s ' % table_name
            qry += "where conflict is null and sample_barcode_short = '%s' " % short_barcode
            qry += "group by chromosome, start_position having count(*) > 1 "
            rows2 = search_db (cursor, qry)
            if rows2:
                for row2 in rows2:
                    [chrom, start, count] = row2
                    print ">>>>>> ", short_barcode, chrom, start, count
                    per_db_cases += 1

            qry =  'select chromosome, start_position, count(*) from %s ' % table_name
            qry += "where conflict is not null and sample_barcode_short = '%s' " % short_barcode
            qry += "group by chromosome, start_position having count(*) > 1 "
            rows2 = search_db (cursor, qry)
            if rows2:
                for row2 in rows2:
                    [chrom, start, count] = row2
                    #print ">>>>>> ", short_barcode, chrom, start, count
                    per_db_cases_conflict+= 1



        print "number of duplicates without conflict annotation: ", per_db_cases
        print "number of duplicates with    conflict annotation: ", per_db_cases_conflict
        print
    
        total_cases += per_db_cases
        total_cases_conflict += per_db_cases_conflict


    print
    print "total number of duplicate cases without conflict annotation: ", total_cases
    print "total number of duplicate cases with    conflict annotation: ", total_cases_conflict
    print

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
