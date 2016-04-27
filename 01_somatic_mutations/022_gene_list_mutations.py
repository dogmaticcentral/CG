#!/usr/bin/python -u
#
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
#

import os
from time import time
import commands
from   tcga_utils.mysql   import  *
from random import randrange, sample
from scipy import stats
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
def main():

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

    #db_names = ["ACC", "UCS"]

    table = 'somatic_mutations'

    pancan_samples = 0
    pancan_ct_gene1 = {}
    pancan_ct = {}
    pancan_pt_gene1 = {}
    pancan_pt = {}
    pancan_coappearance = {}
    number_of_samples_in_dbs_in_which_both_appear = {}
    if tp53_mode:

        gene_list = ['RPL11', 'RPL5', 'MDM2']

        #print "number of different genes:"
        switch_to_db(cursor, 'baseline')
        qry = "select distinct approved_symbol from hgnc_id_translation where locus_type = 'gene with protein product' "
        rows = search_db(cursor, qry)
        gene_list = [row[0] for row in  rows if row[0] != "TP53"]
        print "total genes:", len(gene_list)

        #gene_list = gene_list[:30]
        #gene_list = sample(gene_list, 500)
        gene1 = 'TP53'
        pancan_ct[gene1] = 0
        pancan_pt[gene1] = 0
        for j in range (len(gene_list)):
            gene2 = gene_list[j]
            pancan_ct_gene1[gene2] = 0
            pancan_pt_gene1[gene2] = 0
            pancan_ct[gene2] = 0
            pancan_pt[gene2] = 0
            mut_key = mkey(gene1, gene2)
            pancan_coappearance[mut_key] = 0
            number_of_samples_in_dbs_in_which_both_appear[gene2] = 0
    else:
        gene_list = [ x.upper() for x in sys.argv[1:] ]
        #print gene_list
        for i in range (len(gene_list)):
            gene1 = gene_list[i] 
            pancan_ct[gene1] = 0
            pancan_pt[gene1] = 0
            for j in range (i+1,len(gene_list)):
                gene2 = gene_list[j]
                mut_key = mkey(gene1, gene2)
                pancan_coappearance[mut_key] = 0

    for db_name in db_names:
        header = "\n"

        header += "######################################" + "\n"
        header += " %s  %s " % (db_name, full_name[db_name])
        header += "\n"

        start = time()
        #outf = sys.stdout
        if tp53_mode:
            outf = open ("coapp_tables/%s_tp53_coapps.table" % db_name, "w")
        else:
            outf = open ("coapp_tables/%s_coapps.table" % db_name, "w")
        switch_to_db (cursor, db_name)

        ############################
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        db_entries =  rows[0][0]
        if not rows[0][0]: continue

        ############################
        short_barcodes = []
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

        if tp53_mode:
            mut_ct['TP53'] = 0
            patients_per_gene['TP53'] = 0

        if tp53_mode:
            gene1 = 'TP53'
            for j in range (len(gene_list)):
                gene2   = gene_list[j]
                mut_key = mkey(gene1, gene2)
                co_appearance[mut_key] = 0
        else:
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
                if tp53_mode:
                    if mutations_found.has_key('TP53'):
                        for gene2 in  mutations_found.keys():
                            if gene2=='TP53': continue
                            mut_key = mkey (gene1, gene2)
                            co_appearance[mut_key] += 1

                else:
                    all_mutated_genes_from_the_list = mutations_found.keys();
                    for i in range (len(all_mutated_genes_from_the_list)):
                        gene1 = all_mutated_genes_from_the_list[i]
                        for j in range (i+1, len(all_mutated_genes_from_the_list)):
                            gene2 = all_mutated_genes_from_the_list[j]
                            mut_key = mkey (gene1, gene2)
                            co_appearance[mut_key] += 1

        pancan_samples += number_of_patients
        header +=  "number of different patients: " + str(number_of_patients)+ "\n"
        header +=  "total number of entries: " + str(db_entries)+ "\n"
        header +=  "number of functional mutations (not silent and not 'RNA') " + str(total_muts)+ "\n"
        header +=  " %8s   %4s   %4s %8s  %4s   %4s  %15s  %10s  %10s %10s  " %  ("gene1", "#pts1", "#muts1", "gene2",
                                                                                     "#pts2", "#muts2", "co-appearance",
                                                                                     "expct_co-app", "pval <=", "pval >=")
        #header += "\n"
        outstr = ""
        if tp53_mode:
            gene1 = 'TP53'
            ct1 = mut_ct [gene1]
            pt1 = patients_per_gene[gene1]
            if not pt1: continue
            if float(pt1)/number_of_patients < 0.001: continue

            for j in range (len(gene_list)):
                gene2   = gene_list[j]
                mut_key = mkey (gene1, gene2)
                pancan_coappearance[mut_key] += co_appearance[mut_key]
                appears_together = co_appearance[mut_key]
                ct2 = mut_ct [gene2]
                if not ct2: continue
                pt2 = patients_per_gene[gene2]
                if float(pt2)/number_of_patients < 0.001: continue

                # the number of times gene1 appears in tumors in which both gene1 and gene2 appear
                pancan_ct_gene1[gene2] += ct1
                pancan_pt_gene1[gene2] += pt1

                pancan_ct[gene2] += ct2
                pancan_pt[gene2] += pt2

                number_of_samples_in_dbs_in_which_both_appear [gene2] += number_of_patients

                expctd = float(pt1)/pancan_samples*pt2
                if abs((expctd - appears_together)/expctd) < 0.1: continue
                a = pt2 -  appears_together                     # rpl5 mutated and p53 wt
                b = pancan_samples - pt1 - pt2 + appears_together # rpl5 wt and p53 wt (in pt1 andp2 we subtracted the overlap twice
                c = appears_together                        # rpl5 mutated and p53 mutated
                d = pt1 - appears_together                    # rpl5 wt and p53  mutated
                if expctd > appears_together:
                    # pval to have smaller overlap than expected - that is greater overalp with wt p53 type
                    [odds,pval] = stats.fisher_exact([[a, b], [c, d]], "greater")
                    pval_lt = pval
                    pval_gt = 1.0
                else:
                    [odds,pval] = stats.fisher_exact([[a, b], [c, d]], "less")
                    pval_lt = 1.0
                    pval_gt = pval

                outstr +=  "%8s   %4d   %4d   %8s  %4d    %4d  " %  (gene1, pt1, ct1, gene2, pt2, ct2)
                outstr +=  "%15d    %10.2f   %10.4f  %10.4f" %  ( co_appearance[mut_key], expctd, pval_lt, pval_gt)
                outstr += "\n"

        else:
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
                    ovlp = co_appearance[mut_key]
                    number_of_iterations = 2*number_of_patients
                    ct2 = mut_ct [gene2]
                    if not ct2: continue
                    pt2 = patients_per_gene[gene2]
                    # cmd = "coapp_sim  %d  %d    %d  %d   %d   "  % (number_of_patients, pt1, pt2, ovlp, number_of_iterations)
                    # [avg, pval_le, pval_ge] = [float(x) for x in commands.getoutput(cmd).split()]#
                    a = pt2 - ovlp                     # rpl5 mutated and p53 wt
                    b = number_of_patients - pt1 - pt2 + ovlp # rpl5 wt and p53 wt (in pt1 andp2 we subtracted the overlap twice
                    c = ovlp                           # rpl5 mutated and p53 mutated
                    d = pt1 - ovlp                     # rpl5 wt and p53  mutated
                    [odds,pval] = stats.fisher_exact([[a, b], [c, d]], "greater")
                    outstr += "%8s   %4d   %4d    %8s  %4d  %4d" %  (gene1, pt1, ct1, gene2, pt2, ct2)
                    outstr += "%15d    %10.2f    %10.2f" %  ( co_appearance[mut_key], float(pt1)/number_of_patients*pt2, pval)
                    outstr += "\n"
        if outstr:
            print >> outf, header
            print >> outf, outstr
        outf.close()
        print db_name, "done in %8.2f min" % ( (time() - start)/60 )


    #outf = sys.stdout
    if tp53_mode:
        outf = open ("coapp_tables/pancan_tp53_coapps.table", "w")
    else:
        outf = open ("coapp_tables/pancan_coapps.table", "w")

    print >> outf, "######################################"
    print >> outf, "pan-cancer"
    print >> outf, " %8s   %4s   %4s %8s  %4s   %4s  %15s %15s %10s %10s  " %  ("gene1", "#pts1", "#muts1", "gene2",
                                                                                "#pts2", "#muts2", "co-appearance",
                                                                                 "expct_co-app", "pval <=", "pval >=")

    if tp53_mode:
        gene11 = 'TP53'
        ct1 = pancan_ct[gene1]
        pt1 = pancan_pt[gene1]

        start = time()
        for j in range (len(gene_list)):

            gene2 = gene_list[j]
            mut_key = mkey (gene1, gene2)
            appears_together = pancan_coappearance[mut_key]

            number_of_patients = number_of_samples_in_dbs_in_which_both_appear[gene2]

            ct1 = pancan_ct_gene1[gene2]
            pt1 = pancan_pt_gene1[gene2]
            if not pt1  or float(pt1)/number_of_patients < 0.001: continue

            ct2 = pancan_ct[gene2]
            pt2 = pancan_pt[gene2]
            if not pt2 or float(pt2)/number_of_patients < 0.001: continue

            expctd = float(pt1)/number_of_patients*pt2
            if abs((expctd - appears_together)/expctd) < 0.1: continue
            a = pt2 -  appears_together                     # rpl5 mutated and p53 wt
            b = number_of_patients - pt1 - pt2 + appears_together # rpl5 wt and p53 wt (in pt1 andp2 we subtracted the overlap twice
            c = appears_together                        # rpl5 mutated and p53 mutated
            d = pt1 - appears_together                    # rpl5 wt and p53  mutated
            if expctd > appears_together:
                [odds,pval] = stats.fisher_exact([[a, b], [c, d]], "greater")
                pval_lt = pval
                pval_gt = 1.0
            else:
                [odds,pval] = stats.fisher_exact([[a, b], [c, d]], "less")
                pval_lt = 1.0
                pval_gt = pval
            print >> outf,  "%8s   %4d   %4d   %8s  %4d    %4d  " %  (gene1, pt1, ct1, gene2, pt2, ct2),
            print >> outf,  "%15d    %10.2f   %10.4f  %10.4f" %  ( appears_together, expctd, pval_lt, pval_gt)
            #if not j%10:
            #    print " %4d  time:  %8.2f min" % (j, (time()-start)/60 )
            #    start = time()
    else:
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
                number_of_iterations = 4*pancan_samples
                #[avg, pval_le, pval_ge]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
                #cmd = "coapp_sim  %d  %d    %d  %d   %d   "  % (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
                #[avg, pval_le, pval_ge] = [float(x) for x in commands.getoutput(cmd).split()]
                a = pt2 -  appears_together                     # rpl5 mutated and p53 wt
                b = pancan_samples - pt1 - pt2 + appears_together # rpl5 wt and p53 wt (in pt1 andp2 we subtracted the overlap twice
                c = appears_together                        # rpl5 mutated and p53 mutated
                d = pt1 - appears_together                    # rpl5 wt and p53  mutated
                [odds,pval] = stats.fisher_exact([[a, b], [c, d]], "greater")

                print >> outf,  "%8s   %4d   %4d   %8s  %4d    %4d  " %  (gene1, pt1, ct1, gene2, pt2, ct2),
                print >> outf,  "%15d    %10.2f   %10.2f" %  (appears_together, float(pt1)/pancan_samples*pt2, pval)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

