#!/usr/bin/python -u
#

import os
from time import time
import commands
from   tcga_utils.mysql   import  *
from random import randrange, sample

#########################################
def simulation (M, Nr, Nb, l, number_of_iterations):
    
    avg_number_of_double_labeled = 0
    pval_le = 0.0  # probabilty of being less-or_equal-to
    pval_ge = 0.0  # probabilty of being greater-or-equal-to

    if not number_of_iterations > 0:
        return  [avg_number_of_double_labeled, pval_le, pval_ge]


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
def mkey (gene1, gene2):
    mut_key = ""
    if gene1 < gene2: # alphabetical
        mut_key = gene1 + "_" + gene2
    else:
        mut_key = gene2 + "_" + gene1

    return mut_key


#########################################
def coappearance_stats (cursor, db_list, gene_list, full_name):

    pancan_samples = 0
    pancan_ct = {}
    pancan_coappearance = {}

    for i in range (len(gene_list)):
        gene1 = gene_list[i]
        pancan_ct[gene1] = 0
        for j in range (i+1,len(gene_list)):
            gene2 = gene_list[j]
            mut_key = mkey(gene1, gene2)
            pancan_coappearance[mut_key] = 0

    for db_name in db_list:
        print "######################################"
        print db_name, full_name[db_name]
        start = time()
        outf = sys.stdout
        switch_to_db (cursor, db_name)

        ############################
        print "total number of entries:",
        qry = "select count(1) from  somatic_mutations"

        rows = search_db(cursor, qry)
        db_entries =  rows[0][0]
        print db_entries
        if not rows[0][0]: continue

        ############################
        short_barcodes = []
        tbarcodes_per_patient = {}
        qry  = "select distinct  sample_barcode_short from somatic_mutations "
        rows = search_db(cursor, qry)

        number_of_patients =  len(rows)

        for row in rows:
            short_barcodes.append(row[0])

        co_appearance = {}
        mut_ct = {}
        patients_per_gene = {}
        for gene  in gene_list:
            mut_ct[gene] = 0
            patients_per_gene[gene] = 0

        for i in range (len(gene_list)):
            gene1 = gene_list[i]
            for j in range (i+1, len(gene_list)):
                gene2   = gene_list[j]
                mut_key = mkey (gene1, gene2)
                co_appearance[mut_key] = 0

        total_muts = 0

        ############################
        for sample_barcode_short in short_barcodes:
            ############################
            qry = "select  hugo_symbol, variant_classification, aa_change "
            qry += " from somatic_mutations"
            qry += " where sample_barcode_short  = '%s' " %  sample_barcode_short
            qry += " and not  variant_classification like '%s' " % "silent"
            qry += " and not  variant_classification like '%s' " % "RNA"

            rows = search_db (cursor, qry)
            if not rows: continue

            mutations_found = {}
            for row in rows:
                [ hugo_symbol, variant_classification, aa_change] = row
                if hugo_symbol in gene_list + ['TP53']:
                    # find genes that are mutated, once or twice, doesn't matter
                    mutations_found[hugo_symbol] = True
                    # here keep track of the actual number of mutations
                    mut_ct[hugo_symbol] += 1

            ############################
            if mutations_found:

                total_muts += len(rows)
                for hugo_symbol in mutations_found.keys():
                    patients_per_gene[hugo_symbol] += 1

                # make sure the key is always in the same order
                all_mutated_genes_from_the_list = mutations_found.keys();
                for i in range (len(all_mutated_genes_from_the_list)):
                    gene1 = all_mutated_genes_from_the_list[i]
                    for j in range (i+1, len(all_mutated_genes_from_the_list)):
                        gene2 = all_mutated_genes_from_the_list[j]
                        mut_key = mkey (gene1, gene2)
                        co_appearance[mut_key] += 1

        pancan_samples += number_of_patients
        print >> outf, "number of different patients:", number_of_patients
        print >> outf, "total number of entries:", db_entries
        print >> outf, "number of functional mutations (not silent and not 'RNA')", total_muts
        print >> outf, " %8s   %4s   %4s   %8s  %4s   %4s    %15s    %s" %  ("gene1", "#pts1", "#muts1",
	      	       	  "gene2", "#pts2", "#muts2", "co-appearance", "expected_no_of_co-appearances")

        for i in range (len(gene_list)):
            gene1 = gene_list[i]
            ct1 = mut_ct [gene1]
            pt1 = patients_per_gene[gene1]
            pancan_ct[gene1] += ct1
            for j in range (i+1,len(gene_list)):
                gene2   = gene_list[j]
                mut_key = mkey (gene1, gene2)
                pancan_coappearance[mut_key] += co_appearance[mut_key]
                ct2 = mut_ct [gene2]
                pt2 = patients_per_gene[gene2]
                print >> outf,  "%8s   %4d   %4d   %8s  %4d    %4d  " %  (gene1, pt1, ct1, gene2, pt2, ct2),
                print >> outf,  "%15d    %15.2f" %  ( co_appearance[mut_key], expected (ct1, ct2, number_of_patients))

        print "db done in %8.2f min" % ( (time() - start)/60 )

    print "######################################"
    print "pan-cancer"
    print >> outf, "number of samples:", pancan_samples
    print >> outf, " %8s   %4s   %8s  %4s  %15s  %15s  %15s  %15s  %15s" %  ("gene1", "#muts1", "gene2", "#muts2",
                                                                "co-appearance", "expected no", "expected no", "pval of <=", "pval of >=")
    print >> outf, " %8s   %4s   %8s  %4s  %15s        %15s  %15s  %15s  %15s" %  ("", "", "", "", "", "of co-app (expr)","of co-app (sim)", "", "")

    for i in range (len(gene_list)):
        gene1 = gene_list[i]
        ct1 = pancan_ct[gene1]
        for j in range (i+1,len(gene_list)):
            gene2 = gene_list[j]
            mut_key = mkey (gene1, gene2)
            appears_together = pancan_coappearance[mut_key]
            ct2 = pancan_ct[gene2]
            number_of_iterations = 2*pancan_samples
            #[avg, pval_le, pval_ge]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            cmd = "coapp_sim  %d  %d    %d  %d   %d   "  % (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            [avg, pval_le, pval_ge] = [float(x) for x in commands.getoutput(cmd).split()]
            print >> outf,  " %8s   %4d   %8s  %4d  %15d  %15.2f  %15.2f  %15.4f  %15.4f" %  ( gene1, ct1, gene2, ct2,
                                                                                               appears_together, expected (ct1, ct2, pancan_samples),
                                                                                               avg, pval_le, pval_ge )


#########################################
def mut_freqs (cursor, db_name, gene):
    non_silent_ct = 0
    silent_ct     = 0

    switch_to_db (cursor, db_name)
    qry  = "select count(1) from somatic_mutations "
    qry += "where hugo_symbol='%s' " % gene
    qry += "and variant_classification in ('Missense_Mutation', 'Nonstop_Mutation', 'Nonsense_Mutation')"
    rows = search_db(cursor, qry)
    if not rows:
        non_silent_ct += 0
    else:
        non_silent_ct += rows[0][0]

    qry  = "select count(1) from somatic_mutations "
    qry += "where hugo_symbol='%s' " % gene
    qry += "and variant_classification='silent'"
    rows = search_db(cursor, qry)
    if not rows:
        silent_ct += 0
    else:
        silent_ct += rows[0][0]

    qry  = "select distinct  sample_barcode_short from somatic_mutations "
    rows = search_db(cursor, qry)
    number_of_patients =  len(rows)

    return [silent_ct, non_silent_ct, number_of_patients]

#########################################
def main():

    if len(sys.argv) != 3:
        print "please give me exactly 2 gene names"
        exit(1)
  
    full_name = read_cancer_names ()

    # unbuffered output
    #sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    
    db     = connect_to_mysql()
    cursor = db.cursor()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    gene_list = [ x.upper() for x in sys.argv[1:] ]

    non_silent_freq = {}
    for gene in gene_list:
        non_silent_freq[gene] = {}
        for db_name in db_names:
            [silent_ct, non_silent_ct,  number_of_patients] =  mut_freqs(cursor, db_name, gene)
            non_silent_freq[gene][db_name] = float(non_silent_ct)/number_of_patients
            print  "%5s  %5s  %5d  %4d  %4d  %5.3f " % (gene,  db_name, number_of_patients, non_silent_ct, number_of_patients, non_silent_freq[gene][db_name] )

    # filter databases to have the requisite  minimal frequency of functional mutations
    # in both genes
    [gene1, gene2] = gene_list
    #for min_freq in [ 0.019, 0.009, 0.01, 0.0,  -1.0]:
    for min_freq in [ 0.0]:
        filtered_db_list = [db_name for db_name in db_names if
                            non_silent_freq[gene1][db_name] > min_freq and  non_silent_freq[gene2][db_name] > min_freq]
        print
        print "==================================================================="
        print min_freq, filtered_db_list
        coappearance_stats (cursor, filtered_db_list, gene_list, full_name)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

