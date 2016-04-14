#!/usr/bin/python -u
#

import os
from time import time
import commands
from   tcga_utils.mysql   import  *
from random import randrange, sample
from scipy  import stats

use_metastatic =  False
verbose        =  True
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
    

#########################################
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
def coappearance_stats (cursor, db_list, primary_samples, metastatic_samples, gene_list, full_name):

    pancan_samples = 0
    pancan_ct = {}
    pancan_pt = {}
    pancan_coappearance = {}

    for i in range (len(gene_list)):
        gene1 = gene_list[i]
        pancan_ct[gene1] = 0
        pancan_pt[gene1] = 0
        for j in range (i+1,len(gene_list)):
            gene2 = gene_list[j]
            mut_key = mkey(gene1, gene2)
            pancan_coappearance[mut_key] = 0

    outf = sys.stdout

    for db_name in db_list:

        if verbose:
            print "######################################"
            print db_name, full_name[db_name]
        start = time()
        switch_to_db (cursor, db_name)

        ############################
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
        all_samples = primary_samples[db_name][:]
        if use_metastatic:
            all_samples += metastatic_samples[db_name]
        number_of_patients = len(all_samples)

        for sample_barcode_short in all_samples:
            ############################
            qry = "select  hugo_symbol, variant_classification, aa_change "
            if use_metastatic and sample_barcode_short in metastatic_samples[db_name]:
                qry += "from metastatic_mutations"
            else:
                qry += " from somatic_mutations"
            qry += " where sample_barcode_short  = '%s' " %  sample_barcode_short
            qry += " and not variant_classification in ('Silent', 'RNA') "

            rows = search_db (cursor, qry)
            if not rows: continue

            mutations_found = {}
            for row in rows:
                [ hugo_symbol, variant_classification, aa_change] = row
                if hugo_symbol in gene_list:
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

        ok = True
        for i in range (len(gene_list)):
            gene1 = gene_list[i]
            pt1 = patients_per_gene[gene1]
            if pt1==0:
                ok = False
                break

        if not ok: continue


        if verbose:
            print >> outf, "number of different patients:", number_of_patients
            print >> outf, "number of functional mutations in codons (not silent and not 'RNA')", total_muts
            print >> outf, " %8s   %4s   %4s   %8s  %4s   %4s    %15s    %s" %  ("gene1", "#pts1", "#muts1",
                          "gene2", "#pts2", "#muts2", "co-appearance", "expected_no_of_co-appearances")

        for i in range (len(gene_list)):
            gene1 = gene_list[i]
            ct1 = mut_ct [gene1]
            pt1 = patients_per_gene[gene1]
            pancan_ct[gene1] += ct1
            pancan_pt[gene1] += pt1
            for j in range (i+1,len(gene_list)):
                gene2   = gene_list[j]
                mut_key = mkey (gene1, gene2)
                pancan_coappearance[mut_key] += co_appearance[mut_key]
                ct2 = mut_ct [gene2]
                pt2 = patients_per_gene[gene2]
                if verbose:
                    print >> outf,  "%8s   %4d   %4d   %8s  %4d    %4d  " %  (gene1, pt1, ct1, gene2, pt2, ct2),
                    print >> outf,  "%15d    %15.2f" %  (co_appearance[mut_key], expected (ct1, ct2, number_of_patients))

        if verbose: print "db done in %8.2f min" % ( (time() - start)/60 )

    if pancan_samples==0: return

    if verbose:
        print "######################################"
        print "pan-cancer"
    print >> outf, "number of samples:", pancan_samples
    print >> outf, " %8s   %4s  %4s   %8s  %4s  %4s   %15s  %15s  %15s  %15s " %  ("gene1", "#muts1", "#pts1", "gene2", "#muts2", "#pts2",
                                                                "co-appearance", "expected",  "pval of <=", "pval of >=")
    #print >> outf, " %8s   %4s  %4s   %8s  %4s  %4s   %15s        %15s  %15s  %15s  %15s  %15s" %  ("", "", "", "", "", "", "", "of co-app (expr)","of co-app (sim)", "", "", "")

    for i in range (len(gene_list)):
        gene1 = gene_list[i]
        ct1 = pancan_ct[gene1]
        pt1 = pancan_pt[gene1]
        for j in range (i+1,len(gene_list)):
            gene2 = gene_list[j]
            mut_key = mkey (gene1, gene2)
            appears_together = pancan_coappearance[mut_key]
            ct2 = pancan_ct[gene2]
            pt2 = pancan_pt[gene2]
            #number_of_iterations = 2*pancan_samples
            #[avg, pval_le, pval_ge]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            #cmd = "coapp_sim  %d  %d    %d  %d   %d   "  % (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            #[avg, pval_le, pval_ge] = [float(x) for x in commands.getoutput(cmd).split()]
            a = pt2 -  appears_together                     # rpl5 mutated and p53 wt
            b = pancan_samples - pt1 - pt2 + appears_together # rpl5 wt and p53 wt (in pt1 andp2 we subtracted the overlap twice
            c = appears_together                        # rpl5 mutated and p53 mutated
            d = pt1 - appears_together                    # rpl5 wt and p53  mutated
            [odds,pval_fisher_le] = stats.fisher_exact([[a, b], [c, d]], "greater")
            [odds,pval_fisher_ge] = stats.fisher_exact([[a, b], [c, d]], "less")
            print >> outf,  " %8s   %4d  %4d   %8s  %4d %4d  %15d  %15.2f  %15.2f  %15.4f  "  % ( gene1, ct1, pt1, gene2, ct2, pt2,
                                                                                               appears_together, float(pt1)*pt2/pancan_samples,
                                                                                               pval_fisher_le, pval_fisher_ge)


#########################################
def find_uniq_metastatic (cursor, db_name):
    switch_to_db (cursor, db_name)

    primary_samples = []
    meta_samples    = []
    duplicates      = []
    qry  = "select distinct  sample_barcode_short from somatic_mutations "
    rows = search_db(cursor, qry)
    primary_samples = [row[0] for row in rows]

    primary_sample_patients = [x[:-3] for x in primary_samples]
    if use_metastatic:
        # find patients from metastatic that do not have the primary tissue counterpart
        qry  = "select distinct  sample_barcode_short from metastatic_mutations "
        rows = search_db(cursor, qry)
        meta_samples = []
        if rows:
            meta_samples = [row[0] for row in rows if row[0][:-3] not in primary_sample_patients]
            duplicates   = [row[0] for row in rows if row[0][:-3]  in primary_sample_patients]

    return [primary_samples, meta_samples, duplicates]

#########################################
def get_count(cursor, qry):

    count = 0
    rows = search_db(cursor, qry)
    if not rows:
        count = 0
    else:
        count = rows[0][0]

    return count

#########################################
def mut_freqs (cursor, db_name, metastatic_samples, duplicate_primary_meta_samples, gene):
    non_silent_ct = 0
    silent_ct     = 0
    switch_to_db (cursor, db_name)

    qry  = "select count(1) from somatic_mutations "
    qry += "where hugo_symbol='%s' " % gene

    qry_non    =  qry + "and not variant_classification in ('Silent', 'RNA') "
    non_silent_ct += get_count(cursor,qry_non)

    qry_silent =  qry + "and variant_classification in ('Silent', 'RNA')"
    silent_ct += get_count(cursor,qry_silent)

    if use_metastatic and metastatic_samples:
        qry  = "select count(1) from metastatic_mutations "
        qry += "where hugo_symbol='%s' " % gene
        if duplicate_primary_meta_samples:
            duplicates_string = ', '.join( duplicate_primary_meta_samples)
            qry += "and not find_in_set(sample_barcode_short, '%s') " % duplicates_string

        qry_non    =  qry + "and not variant_classification in ('Silent', 'RNA') "
        non_silent_ct += get_count(cursor,qry_non)

        qry_silent =  qry + "and variant_classification in ('Silent', 'RNA')"
        silent_ct += get_count(cursor,qry_silent)

    return [silent_ct, non_silent_ct]

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

    primary_samples = {}
    meta_samples    = {}
    duplicates      = {}
    for db_name in db_names:
        [primary_samples[db_name], meta_samples[db_name], duplicates[db_name]] = find_uniq_metastatic(cursor, db_name)


    non_silent_freq = {}
    for gene in gene_list:
        non_silent_freq[gene] = {}
        for db_name in db_names:
            [silent_ct, non_silent_ct] =  mut_freqs(cursor, db_name, meta_samples[db_name], duplicates[db_name], gene)
            number_of_patients = len(primary_samples[db_name])
            if use_metastatic: number_of_patients += len(meta_samples)
            non_silent_freq[gene][db_name] = float(non_silent_ct)/number_of_patients
            #if verbose: print  "%5s  %5s  %5d  %4d  %4d  %5.3f " % (gene,  db_name, number_of_patients, non_silent_ct, number_of_patients, non_silent_freq[gene][db_name] )

    # filter databases to have the requisite  minimal frequency of functional mutations
    # in both genes
    [gene1, gene2] = gene_list
    for min_freq in [ 0.02, 0.01, 0.001, 0.0,  -1.0]:
    #for min_freq in [0.0,  -1.0]:
        filtered_db_list = [db_name for db_name in db_names if
                            non_silent_freq[gene1][db_name] > min_freq and  non_silent_freq[gene2][db_name] > min_freq]
        print
        print "==================================================================="
        if min_freq<0:
            print "min freq of nonsilent mutations: any"
        else:
            print "min freq of nonsilent mutations > %.1f%%" %  (min_freq*100)
        print "tumors: ", ", ".join(filtered_db_list)
        coappearance_stats (cursor, filtered_db_list, primary_samples, meta_samples, gene_list, full_name)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

