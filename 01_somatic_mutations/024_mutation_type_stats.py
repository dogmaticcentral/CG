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

import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *

#########################################
def main():

    extra_genes = ["TP53", "RPL5", "RPL11", "MDM2"]
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    table = 'somatic_mutations'

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


    full_name = read_cancer_names ()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    # the order in which we want the variants output:
    variant_order = ["missense_mutation", "silent", "nonsense_mutation", "rna", "splice_site", "frame_fhift_ins",
                     "frame_shift_del", "in_frame_del", "in_frame_ins", "translation_start_site", "nonstop_mutation",
                     "5utr", "3utr", "igr", "intron", "5flank"]
    grand_total = {}
    for variant in variant_order:
        grand_total[variant] = 0

    grand_total_extra = {}
    for gene in extra_genes: 
        grand_total_extra[gene] = {}
        for variant in variant_order:
            grand_total_extra[gene][variant] = 0
        

    total_extra    = {}
    total_patients = 0
    for db_name in db_names:

        total = 0
        for gene in extra_genes: 
            total_extra[gene] = 0
        ############################
        switch_to_db (cursor, db_name)
        print 
        print "################################"
        print db_name, full_name[db_name]
  
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        total = int (rows[0][0])
        print "number of entries:", total

        qry = "select distinct sample_barcode_short from somatic_mutations"
        rows = search_db(cursor, qry)
        number_of_patients = len(rows)
        total_patients    += number_of_patients
        print "number of patients:", number_of_patients
        ############################

        for gene in extra_genes: 
            qry = "select count(1) from " + table
            qry += " where hugo_symbol = '%s' " % gene
            rows = search_db(cursor, qry)
            total_extra[gene] = int (rows[0][0])
            print "\tnumber of entries for", gene, ":",  total_extra[gene]
            


        ############################
        print "variant classification (# cases):  %20s" %  "overall",
        for gene in extra_genes:
            print " %15s " % gene,
        print
        qry = "select distinct(variant_classification) from somatic_mutations"
        rows = search_db(cursor, qry)
        variants = [row[0] for row in  rows]
        for variant in variant_order:
            if not variant in variants: continue
            qry = "select count(1) from somatic_mutations where variant_classification='%s'" % variant
            rows = search_db(cursor, qry)
            print "\t %30s    %4d (%4.1f%%)" %  (variant, rows[0][0], float(rows[0][0])/total*100 ),
            grand_total[variant] += rows[0][0]
            for gene in extra_genes:
                if not total_extra[gene]:
                    print "\t %4d (%4.1f%%)" % (0, 0.0),
                else:
                    qry = "select count(1) from somatic_mutations "
                    qry += " where hugo_symbol = '%s' " % gene
                    qry += " and variant_classification='%s'" % variant
                    rows = search_db(cursor, qry)
                    print "\t %4d (%4.1f%%)" % (rows[0][0], float(rows[0][0])/total_extra[gene]*100 ),
                    grand_total_extra[gene][variant] += rows[0][0]
            print

    grand_grand_total = 0
    for variant in variant_order:
        grand_grand_total += grand_total[variant] 

    grand_grand_total_extra = {}
    for gene in extra_genes:
        grand_grand_total_extra[gene] = 0
        for variant in variant_order:
            grand_grand_total_extra[gene] += grand_total_extra[gene][variant] 
       

    print 
    print "################################"
    print 'pan-cancer'
    print "number of entries:", grand_grand_total
    print "number of patients:", total_patients
    print "variant classification (# cases):  %20s" %  "overall",
    for gene in extra_genes:
        print " %15s " %  gene,
    print
    for variant in variant_order:
        print "\t %30s   %8d (%4.1f%%)" %  (variant,  grand_total[variant], float(grand_total[variant] )/grand_grand_total*100 ),
        for gene in extra_genes:
            if not grand_grand_total_extra[gene]:
                print "\t %4d (%4.1f%%)" %  (0, 0.0),
            else:
                #qry = "select count(1) from somatic_mutations "
                print "\t %4d (%4.1f%%)" %  (grand_total_extra[gene][variant],
                                               float(  grand_total_extra[gene][variant] )/grand_grand_total_extra[gene]*100 ),
        print
        
  

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

