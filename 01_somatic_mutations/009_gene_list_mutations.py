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
def main():

    if len(sys.argv) == 1:
        tp53_mode = True
    else:
        tp53_mode = False
  
    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    
    db     = connect_to_mysql()
    cursor = db.cursor()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    table = 'somatic_mutations'

    pancan_samples = 0
    pancan_ct = {}
    pancan_coappearance = {}
    
    if tp53_mode:

        print "number of different genes:"
        switch_to_db(cursor, 'baseline')
        qry = "select distinct approved_symbol from hgnc_id_translation where locus_type = 'gene with protein product' "
        rows = search_db(cursor, qry)
        gene_list = [row[0] for row in  rows if row[0] != "TP53"]
        print "\t db_name", len(gene_list)

        #gene_list = gene_list[:30]
        #gene_list = sample(gene_list, 500)
        gene1 = 'TP53'
        pancan_ct[gene1] = 0
        for j in range (len(gene_list)):
            gene2 = gene_list[j]
            pancan_ct[gene2] = 0
            mut_key = gene1 + "_"  + gene2
            pancan_coappearance[mut_key] = 0

    else:
        gene_list = [ x.upper() for x in sys.argv[1:] ]
        print gene_list
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
        start = time()
        if tp53_mode:
            outf = open ("coapp_tables/%s_tp53_coapps.table" % db_name, "w")
        else:
            outf = open ("coapp_tables/%s_coapps.table" % db_name, "w")
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

        co_appearance = {}
        mut_ct = {}
        for gene  in gene_list:
            mut_ct[gene] = 0
        if tp53_mode: mut_ct['TP53'] = 0

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

                # make sure the key is always in the same order
                if tp53_mode:
                    if mutations_found.has_key('TP53'):
                        for gene2 in  mutations_found.keys():
                            mut_key = gene1 + "_" + gene2
                            if not  mut_key in co_appearance.keys():
                                co_appearance[mut_key] = 0
                            co_appearance[mut_key] += 1

                else:
                    all_mutated_genes_from_the_list = mutations_found.keys();
                    for i in range (len(all_mutated_genes_from_the_list)):
                        gene1 = all_mutated_genes_from_the_list[i]
                        for j in range (i+1, len(all_mutated_genes_from_the_list)):
                            gene2 = all_mutated_genes_from_the_list[j]
                            mut_key = gene1 + "_" + gene2
                            if not  mut_key in co_appearance.keys():
                                co_appearance[mut_key] = 0
                            co_appearance[mut_key] += 1

        pancan_samples += number_of_patients
        print >> outf, "number of different patients:", number_of_patients
        print >> outf, "total number of entries:", db_entries
        print >> outf, "number of functional mutations (not silent and not 'RNA')", total_muts
        print >> outf, " %8s   %4s   %8s  %4s    %15s    %s" %  ("gene1", "#muts1", "gene2", "#muts2", "co-appearance", "expected_no_of_co-appearances")
        if tp53_mode:
            gene1 = 'TP53'
            ct1 = mut_ct [gene1]
            pancan_ct[gene1] += ct1
            for j in range (len(gene_list)):
                gene2 = gene_list[j]
                mut_key = gene1 + "_" + gene2
                appears_together = 0
                if mut_key in co_appearance.keys():
                    appears_together = co_appearance[mut_key]
                pancan_coappearance[mut_key] += appears_together
                ct2 = mut_ct [gene2]
                pancan_ct[gene2] += ct2
                print >> outf,  " %8s   %4d   %8s  %4d    %15d    %15.2f" %  ( gene1, ct1, gene2, ct2,  appears_together, expected (ct1, ct2, number_of_patients))

        else:
            for i in range (len(gene_list)):
                gene1 = gene_list[i]
                ct1 = mut_ct [gene1]
                pancan_ct[gene1] += ct1
                for j in range (i+1,len(gene_list)):
                    gene2   = gene_list[j]
                    mut_key = gene1 + "_" + gene2
                    appears_together = 0
                    if mut_key in co_appearance.keys():
                         appears_together = co_appearance[mut_key]
                    pancan_coappearance[mut_key] += appears_together
                    ct2 = mut_ct [gene2]
                    print >> outf,  " %8s   %4d   %8s  %4d    %15d    %15.2f" %  ( gene1, ct1, gene2, ct2,  appears_together, expected (ct1, ct2, number_of_patients))

        outf.close()
        print "db done in %8.2f min" % ( (time() - start)/60 )


    if tp53_mode:
        outf = open ("coapp_tables/pancan_tp53_coapps.table", "w")
    else:
        outf = open ("coapp_tables/pancan_coapps.table", "w")

    print "######################################"
    print "pan-cancer"
    print >> outf, "number of samples:", pancan_samples
    print >> outf, " %8s   %4s   %8s  %4s  %15s  %15s  %15s  %15s  %15s" %  ("gene1", "#muts1", "gene2", "#muts2",
                                                                "co-appearance", "expected no", "expected no", "pval of <=", "pval of >=")
    print >> outf, " %8s   %4s   %8s  %4s  %15s        %15s  %15s  %15s  %15s" %  ("", "", "", "", "", "of co-app (expr)",
                                                                      "of co-app (sim)", "", "")
    if tp53_mode:
        gene11 = 'TP53'
        ct1 = pancan_ct[gene1]
        start = time()
        for j in range (len(gene_list)):
            gene2 = gene_list[j]
            mut_key = gene1 + "_"  + gene2
            appears_together = pancan_coappearance[mut_key]
            ct2 = pancan_ct[gene2]
            number_of_iterations = 2*pancan_samples
            #[avg, pval_le, pval_ge]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            cmd = "coapp_sim  %d  %d    %d  %d   %d   "  % (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
            [avg, pval_le, pval_ge] = [float(x) for x in commands.getoutput(cmd).split()]
            print >> outf, " %8s   %4d   %8s  %4d  %15d  %15.2f  %15.2f  %15.4f  %15.4f" %  ( gene1, ct1, gene2, ct2,
                                                                                     appears_together, expected (ct1, ct2, pancan_samples),
                                                                                     avg, pval_le, pval_ge )
            if not j%10:
                print " %4d  time:  %8.2f min" % (j, (time()-start)/60 )
                start = time()
    else:
        for i in range (len(gene_list)):
            gene1 = gene_list[i]
            ct1 = pancan_ct[gene1]
            for j in range (i+1,len(gene_list)):
                gene2 = gene_list[j]
                mut_key = gene1 + "_"  + gene2
                appears_together = pancan_coappearance[mut_key]
                ct2 = pancan_ct[gene2]
                number_of_iterations = 2*pancan_samples
                #[avg, pval_le, pval_ge]  = simulation (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
                cmd = "coapp_sim  %d  %d    %d  %d   %d   "  % (pancan_samples, ct1, ct2, appears_together, number_of_iterations)
                [avg, pval_le, pval_ge] = [float(x) for x in commands.getoutput(cmd).split()]
                print >> outf,  " %8s   %4d   %8s  %4d  %15d  %15.2f  %15.2f  %15.4f  %15.4f" %  ( gene1, ct1, gene2, ct2,
                                                                                         appears_together, expected (ct1, ct2, pancan_samples),
                                                                                         avg, pval_le, pval_ge )
    
    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

