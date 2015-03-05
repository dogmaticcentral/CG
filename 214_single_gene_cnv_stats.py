#!/usr/bin/python -u

import os
from os import walk
import MySQLdb
from tcga_utils.mysql   import  *
from tcga_utils.ensembl  import  *
from time      import  time

import collections

#########################################
def read_cancer_names ():
    full_name= {}
    inf = open ("/Users/ivana/pypeworks/tcga/cancer_names.txt", "r")
    for line in inf:
        line   = line.rstrip() 
        field = line.split ("\t")
        if field[0] == 'READ':
            field[0] = 'REA'
        full_name[field[0]] = field[1]
    inf.close()

    return full_name



    


#########################################
def main():
    
    tcga_dir = '/Users/ivana/databases/TCGA'

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];
    full_name = read_cancer_names ()

    db          = connect_to_mysql()
    cursor      = db.cursor()
    db2         = connect_to_mysql()
    cursor_tcga = db2.cursor()

    db_name = 'homo_sapiens_core_75_37'
    switch_to_db (cursor, db_name)

    seqregion_id2name = get_seq_region_names(cursor)

    if 0:
        qry = "select gene_id, stable_id, description from gene where  biotype='protein_coding' "
        qry += " and description like '%ribosomal protein L%' "
        rows = search_db (cursor, qry)
        if not rows:
            print "qry failed: "
            print qry
            exit(1)

        qry = "select gene_id, stable_id, description from gene where  biotype='protein_coding' "
        qry += " and description like '%ribosomal protein L%' "
        rows2 = search_db (cursor, qry)
        if not rows2:
            print "qry failed: "
            print qry
            exit(1)
    else:
        qry = "select gene_id, stable_id, description from gene where  biotype='protein_coding' "
        qry += " and description like '%breast cancer 1%' "
        rows = search_db (cursor, qry)
        if not rows:
            print "qry failed: "
            print qry
            exit(1)

    for row in rows:
        output =  "\n\n######################################\n"
        [gene_id, stable_id, description] = row
        output += " %s, %s, %s \n" % ( gene_id, stable_id, description)

        qry = "select seq_region_id, seq_region_start, seq_region_end from gene where gene_id = %d " % gene_id
        rows = search_db (cursor, qry)
        [seq_region_id, seq_region_start, seq_region_end] = rows[0]
        gene_length = seq_region_end-seq_region_start

        qry = "select length from seq_region where seq_region_id = %d" % seq_region_id
        rows = search_db (cursor, qry)
        chrom_length = rows[0][0]

        output += " start:%d   end:%d   length:%d \n" % (seq_region_start, seq_region_end, seq_region_end-seq_region_start)
        output += " chromosome: %s    chrom. length: %d " %  (seqregion_id2name[seq_region_id], chrom_length)

        hit = False
        for db_name in db_names:
            switch_to_db (cursor_tcga, db_name)
            qry = "select count(1) from cnv_source_ids"
            rows = search_db (cursor_tcga, qry)
            count = rows[0][0]

            qry = "select * from cnv_snp where gene_stable_id = '%s' " % stable_id
            rows = search_db (cursor_tcga, qry)
            if rows:
                hit = True
                output += "\n%s, %s   (tot samples used: %d)\n" % (db_name, full_name[db_name], count)
                output += "sample_id   log2fold   num_probes    start        end         length     length/gene_length    length/chrom_length\n"
                for row in rows:
                    (db_id, stable, src, folds, num_probes, rngs) = row
                    sources = src.rsplit(';')
                    fold_changes =  folds.rsplit(';')
                    num_pbs      =  num_probes.rsplit(';')
                    ranges       =  rngs.rsplit(';')
                    for i in range (len(sources)):
                        start, end = ranges[i].rsplit("-")
                        start = int(start)
                        end   = int(end)
                        length = end-start
                       
                        output +=  "     %4d      %5s     %5d   %10d  %10d   " %  ( int(sources[i]),  fold_changes[i],  
                                                                                        int(num_pbs[i]),   start,  end,  )
                        output +=  " %10d         %8.2f            %8.2f\n" %  ( length, float(length)/gene_length, float(length)/chrom_length )
                                                                                        
        if hit: print output
 
    cursor.close()
    cursor_tcga.close()
    db.close()
    db2.close()


#########################################
if __name__ == '__main__':
    main()
