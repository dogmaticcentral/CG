#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);

# extracting mutations ('catch' is the output from this script)
# grep Missense catch | grep -v found | grep RPL5 | grep -v Silent | awk '{print $5}' | sed 's/p\.//g' | sed 's/[A-Z]//g' | awk '{printf "%d+", $1}' && echo


import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from random import randrange

#########################################
def simulation (M, Nr, Nb, l, number_of_iterations):
    
    avg_number_of_double_labeled = 0
    pval = 0.0

    if not number_of_iterations > 0:
        return  [avg_number_of_double_labeled, pval]


    for i in range(number_of_iterations):
        #####
        slots = []
        for s in range(M):
            slots.append({"r":0, "b":0})
        number_of_double_labeled = 0
        for j in range(Nr):
            random_slot = randrange(M)
            slots[random_slot]["r"] += 1
        for j in range(Nb):
            random_slot = randrange(M)
            slots[random_slot]["b"] += 1

        for s in range(M):
            if slots[s]["r"]>0  and  slots[s]["b"]>0:
                #print " %3d   %2d  %2d " %  (s, slots[s]["r"] ,  slots[s]["b"])
                number_of_double_labeled += 1

        #####
        avg_number_of_double_labeled += number_of_double_labeled
        if ( number_of_double_labeled <= l ): pval += 1.0

    ##################################
    avg_number_of_double_labeled /= float(number_of_iterations)
    pval /= float(number_of_iterations)

    return [avg_number_of_double_labeled, pval]

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
def main():

    if len(sys.argv) < 2:
        print  "usage: %s <gene symbol 1>  <gene symbol 2> ..." % sys.argv[0]
        exit(1)

    gene_list = [ x.upper() for x in sys.argv[1:] ]
    print gene_list
  
    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "COAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]

   
    table = 'somatic_mutations'

    pancan_samples = 0
    pancan_ct = {}
    pancan_coappearance = {}
    for i in range (len(gene_list)):
        gene1 = gene_list[i] 
        pancan_ct[gene1] = 0
        for j in range (i+1,len(gene_list)):
            gene2 = gene_list[j]
            mut_key = gene1 + "_"  + gene2
            pancan_coappearance[mut_key] = 0

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

        number_of_patients =  len(uniq_patients)
        print "number of different patients:", number_of_patients
 
        mut_breakdown = {}
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
                # make sure the key is always in the same order
                for gene  in gene_list:
                    for mut in mutations_found:
                        [ hugo_symbol, variant_classification, aa_change] = mut
                        if hugo_symbol == gene and not gene in all_mutated_genes_from_the_list:
                            all_mutated_genes_from_the_list.append(hugo_symbol)

                if len(all_mutated_genes_from_the_list)==1:
                    mut_key = all_mutated_genes_from_the_list[0]
                    if not  mut_key in mut_breakdown.keys():
                        mut_breakdown[mut_key] = 0
                    mut_breakdown[mut_key] += 1

                else:   # now disregard the triple and up mutants, and count them as a bunch of doubles

                    for i in range (len(all_mutated_genes_from_the_list)):
                        gene1 = all_mutated_genes_from_the_list[i]
                        for j in range (i+1, len(all_mutated_genes_from_the_list)):
                            gene2 = all_mutated_genes_from_the_list[j]
                            mut_key = gene1 + "_" + gene2
                            if not  mut_key in mut_breakdown.keys():
                                mut_breakdown[mut_key] = 0
                            mut_breakdown[mut_key] += 1

        pancan_samples += number_of_patients
        print "total", total_muts
        print " %8s   %4s   %8s  %4s    %15s    %s" %  ("gene1", "#muts1", "gene2", "#muts2", "co-appearance", "expected_no_of_co-appearances")
        for i in range (len(gene_list)):
            gene1 = gene_list[i] 
            ct1 = mut_ct [gene1]
            pancan_ct[gene1] += ct1
            for j in range (i+1,len(gene_list)):
                gene2 = gene_list[j]
                mut_key = gene1 + "_"  + gene2
                appears_together = 0
                if mut_key in mut_breakdown.keys():
                     appears_together = mut_breakdown[mut_key]
                pancan_coappearance[mut_key] += appears_together
                ct2 = mut_ct [gene2]
                print " %8s   %4d   %8s  %4d    %15d    %15.2f" %  ( gene1, ct1, gene2, ct2,  appears_together, expected (ct1, ct2, number_of_patients))
     
            
    print "######################################"
    print "pan-cancer"
    print "number of samples:", pancan_samples
    print " %8s   %4s   %8s  %4s  %15s  %15s  %15s  %15s  %15s" %  ("gene1", "#muts1", "gene2", "#muts2", 
                                                                "co-appearance", "expected no", "expected no", "pval", "1-pval")
    print " %8s   %4s   %8s  %4s  %15s  %15s  %15s  %15s  %15s" %  ("", "", "", "", "", "of co-app (expr)", 
                                                                      "of co-app (sim)", "", "")
    for i in range (len(gene_list)):
        for j in range (i+1,len(gene_list)):
            gene1 = gene_list[i] 
            gene2 = gene_list[j]
            mut_key = gene1 + "_"  + gene2
            appears_together = pancan_coappearance[mut_key]
            ct1 = pancan_ct[gene1] 
            ct2 = pancan_ct[gene2]
            number_of_iterations = 2*pancan_samples
            [avg, prob]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            print " %8s   %4d   %8s  %4d  %15d  %15.2f  %15.2f  %15.4f  %15.4f" %  ( gene1, ct1, gene2, ct2, 
                                                                               appears_together, expected (ct1, ct2, pancan_samples),
                                                                               avg, prob, 1-prob )
    
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

