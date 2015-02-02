#!/usr/bin/python -u
#

import sys, os
import MySQLdb
import subprocess
from   tcga_utils.mysql   import  *
from random import randrange, random
from math import log10


#########################################
def read_cancer_names ():
    full_name= {}
    inf = open ("/Users/ivana/Dropbox/Sinisa/ribosomal/html/cancer_names.txt", "r")
    for line in inf:
        line   = line.rstrip() 
        field = line.split ("\t")
        if field[0] == 'READ':
            field[0] = 'REA'
        full_name[field[0]] = field[1]
    inf.close()

    return full_name
    

#########################################
def resolve_name(cursor, gene):

    name = ""
    # assume we are using 'baseline' databse

    if 'ENSG0' in gene:
        qry = "select locus_type from hgnc_id_translation where "
        qry += " ensembl_gene_id  = '%s'" % gene
    else:
        qry = "select locus_type from hgnc_id_translation where "
        qry += " approved_symbol = '%s'" % gene
        
    rows = search_db (cursor, qry)
    if rows:        
        locus_type = rows[0][0]
        if not locus_type == 'gene with protein product': return ""
        name = gene

    return name

#########################################
def main():

  
    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    
    
    db     = connect_to_mysql()
    cursor = db.cursor()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];

  
    table = 'somatic_mutations'
    # find all genes that can be found in the databases
    tmp_gene_list = set([])
    for db_name in db_names:
        switch_to_db (cursor, db_name)
        qry = "select distinct(hugo_symbol) from somatic_mutations whenr"
        rows = search_db(cursor, qry)
        tmp_gene_list |= set([row[0] for row in  rows])

    # cleanup
    if True:
        gene_list = []
        switch_to_db (cursor, 'baseline')
        for gene in tmp_gene_list:        
            name = resolve_name (cursor, gene)
            if not name: continue
            #if random() < 0.01:
            gene_list.append(gene)
    else:
        gene_list = ['TP53']
    
    
    cancers_in_which_gene_appears ={}
    for gene in gene_list:
        cancers_in_which_gene_appears[gene] = []

    pancan_samples = 0
    pancan_ct = {}
    pancan_coappearance = {}
    for  gene in  gene_list:
        pancan_ct[gene] = 0
        pancan_coappearance[gene] = 0

    number_of_patients = {}


    mut_codon_nonsilent = {}
    mut_codon_silent    = {}
    for gene  in gene_list:
        mut_codon_nonsilent[gene] = 0
        mut_codon_silent[gene] = 0


    for db_name in db_names:
        #print "######################################"
        #print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)


        ############################
        uniq_patients = []
        tbarcodes_per_patient = {}
        qry  = "select distinct  tumor_sample_barcode from somatic_mutations "
        rows = search_db(cursor, qry)
        for  row in rows:
            tbarcode = row[0]
            # the fields are 
            # project - tissue source site (TSS)  - participant -
            # sample.vial - portion.analyte  - plate - (sequencing or characterization center)
            fields = tbarcode.split('-')
            patient = '-'.join(fields[1:3])
            if not patient in  uniq_patients:
                uniq_patients.append(patient)
                tbarcodes_per_patient[patient] = []
            tbarcodes_per_patient[patient].append(tbarcode)

            
        ############################
        # one more round: get rid of patients that have multiple samples
        # there are relatively few of those and I don't know what they mean 
        ok_patients = []
        for patient in uniq_patients:
            uninterpretable = False
            sample = ""
            for tbarcode in tbarcodes_per_patient[patient]:
                fields = tbarcode.split('-')
                new_sample = fields[3][:2]
                if not sample:
                    sample = new_sample
                elif sample != new_sample:
                    uninterpretable = True
                    break
            if uninterpretable: continue
            ok_patients.append(patient)

        number_of_patients[db_name] = len(ok_patients)
        #print "number of different patients:", number_of_patients[db_name]
 
        ############################
        for patient in ok_patients:
            tbarcodes = tbarcodes_per_patient[patient]
            mutations_found = []
            for tbarcode in tbarcodes:
                ############################
                qry = "select hugo_symbol, variant_classification, aa_change "
                qry += " from somatic_mutations"
                qry += " where tumor_sample_barcode  = '%s' " % tbarcode
                #qry += " and not  variant_classification like '%s' " % "silent"

                rows = search_db (cursor, qry)
                if not rows: 
                    continue
                
                for row in rows:
                    [ hugo_symbol, variant_classification, aa_change] = row
                    if hugo_symbol in gene_list:
                        if variant_classification == 'Silent':
                            mut_codon_silent[hugo_symbol] += 1
                        elif variant_classification == 'Missense_Mutation':
                            mut_codon_nonsilent[hugo_symbol] += 1
                        elif variant_classification == 'Nonsense_Mutation':
                             mut_codon_nonsilent[hugo_symbol] += 1
                    
    print " %15s    %7s    %8s   %8s " %   ("name", "tot_muts", "nonsilent_fract", "silent_fract")                     
    for gene in gene_list:
        tot =  float (mut_codon_nonsilent[gene] +  mut_codon_silent[gene]) 
        if not tot: continue
        print " %15s    %7d    %8.3f   %8.3f  " % ( gene, tot, mut_codon_nonsilent[gene]/tot, mut_codon_silent[gene]/tot)
    
    
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

