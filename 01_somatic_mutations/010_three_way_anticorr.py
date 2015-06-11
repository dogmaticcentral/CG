#!/usr/bin/python -u
#
# we'll rely on the output from 009 - proteins anticorr with tp53


import os, sys
from time import time
import commands
from   tcga_utils.mysql   import  *
from random import randrange, sample


#########################################
def simulation (number_of_slots, number_of_marbles_of_type, l, number_of_iterations):

    avg_number_of_single_labeled = 0
    avg_number_of_double_labeled = 0
    pval_le = 0.0  # probability of being less-or_equal-to
    pval_ge = 0.0  # probability of being greater-or-equal-to

    if not number_of_iterations > 0:
        return  [avg_number_of_single_labeled, pval_le, pval_ge]

    number_of_marble_types = len(number_of_marbles_of_type)

    for i in range(number_of_iterations):
        if not i%100: print i
        slot = []
        for s in range(number_of_slots):
            slot.append([0]*number_of_marble_types)
        #####
        for a in range(number_of_marble_types):
            for j in  range(number_of_marbles_of_type[a]):
                random_slot = randrange(number_of_slots)
                slot[random_slot][a] += 1
        #####
        number_of_single_labeled = 0
        for s in range(number_of_slots):
            number_of_types_in_the_slot = 0
            for a in range(number_of_marble_types):
                if slot[s][a] > 0: number_of_types_in_the_slot += 1
            if number_of_types_in_the_slot==1:
                number_of_single_labeled += 1

        #####
        avg_number_of_single_labeled += number_of_single_labeled
        if ( number_of_single_labeled <= l ): pval_le += 1.0
        if ( number_of_single_labeled >= l ): pval_ge += 1.0

    ##################################
    avg_number_of_single_labeled /= float(number_of_iterations)
    avg_number_of_double_labeled /= float(number_of_iterations)
    pval_le /= float(number_of_iterations)
    pval_ge /= float(number_of_iterations)


    return [avg_number_of_single_labeled, pval_le, pval_ge]

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
def main():

    print simulation(8554, [3692,230], 3271, 5000)
    exit(1)

    if len(sys.argv) == 1:
        tp53_mode = True
    else:
        tp53_mode = False
  
    full_name = read_cancer_names ()

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    table = 'somatic_mutations'

    pancan_samples = 0
    pancan_ct = {}
    pancan_avoidance = {}
    
    special_list = ['TP53', 'RPL5']
    special_list = ['TP53']

    cmd =  "awk '$1==\"TP53\" && $8<0.066 {print $3}' coapp_tables/pancan_tp53_coapps.table"
    gene_list  = [x for x in  commands.getoutput(cmd).splitlines() if x not in special_list ]


    print "number of different genes considered:", len(gene_list)
    for gene in gene_list + special_list:
        pancan_ct[gene] = 0
        pancan_avoidance[gene] = 0

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        start = time()
        outf = sys.stdout
        switch_to_db (cursor, db_name)

        ############################
        print "total number of entries:", 
        qry = "select count(1) from " + table
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
        print "number of different patients:", number_of_patients


        for row in rows:
            short_barcodes.append(row[0])

        avoidance = {}
        mut_ct = {}
        for gene  in gene_list+ special_list:
            avoidance[gene] = 0
            mut_ct[gene] = 0

        total_muts = 0
        ############################
        for sample_barcode_short in short_barcodes:

            ############################
            qry = "select  hugo_symbol, variant_classification, aa_change "
            qry += " from somatic_mutations"
            qry += " where sample_barcode_short  = '%s' " % sample_barcode_short
            qry += " and not  variant_classification like '%s'" % "silent"
            qry += " and not  variant_classification like '%s'" % "RNA"

            rows = search_db (cursor, qry)
            if not rows: continue

            mutations_found = {}
            for row in rows:
                [ hugo_symbol, variant_classification, aa_change] = row
                if hugo_symbol in gene_list + special_list:
                    # find genes that are mutated, once or twice, doesn't matter
                    mutations_found[hugo_symbol] = True
                    # here keep track of the actual number of mutations
                    mut_ct[hugo_symbol] += 1
            if not mutations_found: continue

            ############################
            total_muts += len(rows)
            # look for cases where exactly one of the three genes appear
            for gene in gene_list:
                exclusion_list = [gene] + special_list
                for i in range(len(exclusion_list)):
                    cyclic_perm = exclusion_list[i:] + exclusion_list[:i]
                    exclusive = mutations_found.has_key(cyclic_perm[0])
                    for j in range (1, len(exclusion_list)):
                        exclusive = exclusive and not mutations_found.has_key(cyclic_perm[j])
                    if exclusive:
                        avoidance[gene] += 1
                        break

        pancan_samples += number_of_patients
        for gene  in gene_list + special_list:
            pancan_ct[gene] +=  mut_ct[gene]
            pancan_avoidance[gene] += avoidance[gene]

        print >> outf, "number of different patients:", number_of_patients
        print >> outf, "total number of entries:", db_entries
        print >> outf, "number of functional mutations (not silent and not 'RNA')", total_muts
        #print >> outf, " %8s   %4s   %8s  %4s   %8s  %4s " %  ("gene1", "#muts1", "gene2", "#muts2", "gene2", "#muts3"),
        #print >> outf, "  %10s    %s" %  ("avoidance", "expected_no_of_avoidance_cases")

        for gene in  avoidance.keys():
            all_mutated = True
            for g in [gene] + special_list:
                if not mut_ct[g]:  all_mutated=False
            if not all_mutated: continue
            for g in special_list + [gene]:
                print >> outf,  " %8s   %4d " % (g, mut_ct[g]),
            print >> outf,  "  %10d " %  ( avoidance[gene])

        print "db done in %8.2f min" % ( (time() - start)/60 )


    outf  = sys.stdout
    print "######################################"
    print "pan-cancer"
    print >> outf, "number of samples:", pancan_samples
    #print >> outf, " %8s   %4s   %8s  %4s   %8s  %4s " %  ("gene1", "#muts1", "gene2", "#muts2", "gene2", "#muts3"),
    #print >> outf, "  %10s   %15s  %15s  %15s " %  ("avoidance",   "expected no", "pval of <=", "pval of >=")

    for gene in  pancan_avoidance.keys():

        all_mutated = True
        for g in [gene] + special_list:
            if not pancan_ct[g]:  all_mutated=False
        if not all_mutated: continue
        counts = []
        for g in special_list + [gene]:
            counts.append(pancan_ct[g])
        [avg_number_of_single_labeled, pval_le, pval_ge] = simulation (pancan_samples, counts,  pancan_avoidance[gene], 300)

        for g in special_list + [gene]:
            print >> outf,  " %10s %4d " % (g, pancan_ct[g]),
        
        print >> outf,  "  %10d   %15.2f  %15.4f  %15.4f " %  (pancan_avoidance[gene], avg_number_of_single_labeled, pval_le, pval_ge)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

