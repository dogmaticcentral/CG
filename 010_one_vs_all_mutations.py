#!/usr/bin/python -u
#

import sys, os
import MySQLdb
import subprocess
from   tcga_utils.mysql   import  *
from random import randrange, random
from math import log10

#########################################
def simulation (M, Nr, Nb, l, number_of_iterations):
    
    avg_number_of_double_labeled = 0
    pval_le = 0.0  # probabilty of being less-or_equal-to
    pval_ge = 0.0  # probabilty of being greater-or-equal-to

    if not number_of_iterations > 0:
        return  [avg_number_of_double_labeled, pval_le, pval_ge]


    for i in range(number_of_iterations):
        #if not i%1000: print "it", i
        #####
        slots = []
        for s in range(M):
            slots.append({"r":0, "b":0})

        
        for j in range(Nr):
            random_slot = randrange(M)
            slots[random_slot]["r"] += 1
        for j in range(Nb):
            random_slot = randrange(M)
            slots[random_slot]["b"] += 1
        number_of_double_labeled = 0
        for s in range(M):
            if slots[s]["r"]>0  and  slots[s]["b"]>0:
                #print " %3d   %2d  %2d " %  (s, slots[s]["r"] ,  slots[s]["b"])
                number_of_double_labeled += 1

        #####
        avg_number_of_double_labeled += number_of_double_labeled
        if ( number_of_double_labeled <= l ): pval_le += 1.0
        if ( number_of_double_labeled >= l ): pval_ge += 1.0

    ##################################
    avg_number_of_double_labeled /= float(number_of_iterations)
    pval_le /= float(number_of_iterations)
    pval_ge /= float(number_of_iterations)

    return [avg_number_of_double_labeled, pval_le, pval_ge]

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
    

def expected (a, b, n):
    expected = 0
    # probability that there are  no mutations of  type a
    p_a = 1.0
    for i in range(a): 
        p = (1-1.0/n)
        p_a *= p
    # probability that there  are  no mutations of  type b
    p_b = 1.0
    for i in range(b): 
        p = (1-1.0/n)
        p_b *= p
    
    # expected number of co-ocurrences of a and b
    expected = (1-p_a)*(1-p_b)*n
    #if a > 0 and b > 0 :
    #    print ">>>>>>  %3d %3d %3d   %5.2f  %5.2f  %5.2f  " % ( a, b, n, 1-p_a, 1-p_b, expected)
    
    return expected

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

    if len(sys.argv) < 2:
        print  "usage: %s <gene symbol 1> " % sys.argv[0]
        exit(1)

    gene_qry = sys.argv[1].upper()
    print gene_qry
  
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
    
    if not gene_qry in gene_list: gene_list.append(gene_qry)

    # ?? how many?
    print len(gene_list)
    print (gene_qry in gene_list)
    # for  gene in  gene_list:
    #     print gene
    #exit(1)
    
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

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)

        ############################
        print "total number of entries:", 
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print  rows[0][0]

        if not rows[0][0]: continue

        ############################
        uniq_patients = {}
        tbarcodes_per_patient = {}
        qry  = "select distinct  tumor_sample_barcode from somatic_mutations "
        rows = search_db(cursor, qry)
        for  row in rows:
            tbarcode = row[0]
            # the fields are 
            # project - tissue source site (TSS)  - participant -
            # source.vial - portion.analyte  - plate - (sequencing or charcterization center)
            fields = tbarcode.split('-')
            patient = '-'.join(fields[1:3])
            if not  uniq_patients.has_key(patient):
                uniq_patients[patient] = []
                uniq_patients[patient].append('-'.join(fields[3:]))
                tbarcodes_per_patient[patient] = []
            tbarcodes_per_patient[patient].append(tbarcode)

        number_of_patients[db_name] = len(uniq_patients)
        print "number of different patients:", number_of_patients[db_name]
 
        co_appearance = {}
        mut_ct = {}
        for gene  in gene_list:
            mut_ct[gene] = 0
        total_muts = 0
        ############################
        for patient in uniq_patients:
            tbarcodes = tbarcodes_per_patient[patient]
            mutations_found = []
            for tbarcode in tbarcodes:
                ############################
                qry = "select  hugo_symbol, variant_classification, aa_change "
                qry += " from somatic_mutations"
                qry += " where tumor_sample_barcode  = '%s' " % tbarcode
                qry += " and not  variant_classification like '%s' " % "silent"

                rows = search_db (cursor, qry)
                if not rows: 
                    print " >>>>>> "
                    continue
                
                for row in rows:
                    [ hugo_symbol, variant_classification, aa_change] = row
                    if hugo_symbol in gene_list:
                        mutations_found.append (row)
                        mut_ct[hugo_symbol] += 1
                        
            ############################
            if mutations_found:
                total_muts += len(rows)
                all_mutated_genes_from_the_list = []
                # find genes that are mutated, once or twice doesn't matter
                for mut in mutations_found:
                    [ hugo_symbol, variant_classification, aa_change] = mut
                    if not hugo_symbol in all_mutated_genes_from_the_list:
                        all_mutated_genes_from_the_list.append(hugo_symbol)

                for gene2  in all_mutated_genes_from_the_list:
                    if db_name not in cancers_in_which_gene_appears[gene2]:
                         cancers_in_which_gene_appears[gene2].append(db_name)


                # now disregard the triple and up mutants, and count them as a bunch of doubles
                if gene_qry in all_mutated_genes_from_the_list:
                    for gene2  in all_mutated_genes_from_the_list:
                        if gene2 ==  gene_qry: continue
                        mut_key = gene2
                        if not  mut_key in co_appearance.keys():
                            co_appearance[mut_key] = 0
                        co_appearance[mut_key] += 1

        pancan_samples += number_of_patients[db_name]
        #print "total", total_muts
        #print " %8s   %4s   %8s  %4s    %15s    %s" %  ("gene1", "#muts1", "gene2", "#muts2", "co-appearance", "expected_no_of_co-appearances")
        
        for gene2 in gene_list:
            pancan_ct[gene2] += mut_ct [gene2]
            if gene2 ==  gene_qry: continue
            mut_key = gene2
            appears_together = 0
            if mut_key in co_appearance.keys():
                appears_together = co_appearance[mut_key]
                pancan_coappearance[mut_key] += appears_together
                ct2 = mut_ct [gene2]
                #print " %8s   %4d   %8s  %4d    %15d    %15.2f" %  ( gene_qry, mut_ct [gene_qry], gene2, ct2,  appears_together, expected (mut_ct [gene_qry], ct2, number_of_patients))
        
            
    print "######################################"
    print "pan-cancer"
    print "number of samples:", pancan_samples
    print " %8s   %4s   %8s  %4s  %15s  %15s  %15s  %15s  %15s" %  ("gene1", "#muts1", "gene2", "#muts2", 
                                                                "co-appearance", "expected no", "expected no", "-log10(pval <=)", "-log10(pval >=)")
    print " %8s   %4s   %8s  %4s  %15s        %15s  %15s  %15s  %15s" %  ("", "", "", "", "", "of co-app (expr)", 
                                                                      "of co-app (sim)", "", "")


    gene1 = gene_qry
    ct = 0
    under_pval = 0
    for gene2 in gene_list:
        if gene2 ==  gene_qry: continue
        ct += 1
        mut_key = gene2
        appears_together = pancan_coappearance[mut_key]
        ct1 = pancan_ct[gene1] 
        ct2 = pancan_ct[gene2]

        if ct2 < 20: continue

        effective_number_of_samples = 0
        for db_name in db_names:
            if db_name in cancers_in_which_gene_appears[gene_qry] and db_name in cancers_in_which_gene_appears[gene2] :
                effective_number_of_samples += number_of_patients[db_name]

        if effective_number_of_samples == 0: continue
        

        expctd =  expected (ct1, ct2, effective_number_of_samples)
        #print " %8s   %4d   %8s  %4d  %15d  %15.2f    %15.2f " %  ( gene1, ct1, gene2, ct2, appears_together, expctd, float(expctd+1)/(appears_together+1))
        #print "\t\t\t 0.5 < float(expctd+1)/(appears_together+1) < 1.3 ?", 0.5 < float(expctd+1)/(appears_together+1) < 1.3
        #if  0.5 < float(expctd+1)/(appears_together+1) < 1.3: continue
        number_of_iterations = 5*effective_number_of_samples
        #[avg, pval_le, pval_ge]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
        #print avg, pval_le, pval_ge
        # use the c implementation instaed
        cmd = "meta/ovlp_sim  %d  %d %d  %d  %d " %  (effective_number_of_samples,  ct1,   ct2,  appears_together, number_of_iterations)
        [avg, eval_le, eval_ge] = [float(x) for x in subprocess.check_output(cmd, shell=True).split()]

        #if (pval_le > 0.02 and pval_ge > 0.001): continue
        if (eval_le <3 and eval_ge <4): continue
        under_pval += 1
        print "  %8s   %4d   %8s  %4d  %15d  %15.2f  %15.2f  %15.4f  %15.4f   (%d out of %d , %d) " %  \
            ( gene1, ct1, gene2, ct2, appears_together, expctd, avg, eval_le, eval_ge,  under_pval,  ct, effective_number_of_samples)
    
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

