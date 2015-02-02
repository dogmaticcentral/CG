#!/usr/bin/python -u
#

import sys, os
import MySQLdb
import subprocess
from   tcga_utils.mysql   import  *
from random import randrange, random
from math   import log10
from time   import  time
#########################################
def read_intersting_genes ():
    genes = []
    #inf = open ("/Users/ivana/pypeworks/tcga/test.txt", "r")
    inf = open ("/Users/ivana/pypeworks/tcga/pct_silent.txt", "r")
    for line in inf:
        line   = line.rstrip() 
        if line[0] == '%': continue
        (name, tot, non_silent, silent) = line.split ()
        if int(tot) <= 30 : continue
        if float(non_silent) <= 0.8: continue
        #if int(tot) <= 50 : continue
        #if float(non_silent) <= 0.85: continue
        genes.append(name)
    inf.close()

    return genes
    


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


    full_name = read_cancer_names ()

    db     = connect_to_mysql()
    cursor = db.cursor()


    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS", "UVM"];

  
    table = 'somatic_mutations'


    # read in  all genes that appear  the databases that we work with
    # and slect ones where the fraction of non-silent mutations is > 0.76 (see script 010_ok_genes_in_ok_sample.py)
    gene_list = read_intersting_genes()
    
    print >> sys.stderr, "number of genes:", len(gene_list)
    
    cancers_in_which_gene_appears ={}
    for gene in gene_list:
        cancers_in_which_gene_appears[gene] = []

    

    pancan_samples = 0
    pancan_ct = {}
    pancan_coappearance = {}
    for  gene in  gene_list:
        pancan_ct[gene] = 0

    for i in range(len(gene_list)):
        gene1 = gene_list[i]
        for j in range(i+1,len(gene_list)):
            gene2 = gene_list[j]
            if gene2 > gene1:
                mut_key = gene1+"_"+gene2
            else:
                mut_key = gene2+"_"+gene1

            pancan_coappearance[mut_key] = 0

    number_of_patients = {}

    for db_name in db_names:
        print >> sys.stderr, "######################################"
        print >> sys.stderr, db_name, full_name[db_name]
        start = time()
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
 
 
        co_appearance = {}
        for i in range(len(gene_list)):
            gene1 = gene_list[i]
            for j in range(i+1,len(gene_list)):
                gene2 = gene_list[j]
                if gene2 > gene1:
                    mut_key = gene1+"_"+gene2
                else:
                    mut_key = gene2+"_"+gene1
                co_appearance[mut_key] = 0

        mut_ct = {}
        for gene  in gene_list:
            mut_ct[gene] = 0
        total_muts = 0
        ############################
        for patient in ok_patients:
            tbarcodes = tbarcodes_per_patient[patient]
            mutations_found = []
            for tbarcode in tbarcodes:
                ############################
                qry = "select hugo_symbol, variant_classification, aa_change "
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
                # find genes that are mutated, once or twice doesn't matter
                for mut in mutations_found:
                    [ hugo_symbol, variant_classification, aa_change] = mut
                    if not hugo_symbol in all_mutated_genes_from_the_list:
                        all_mutated_genes_from_the_list.append(hugo_symbol)

                for gene  in all_mutated_genes_from_the_list:
                    if db_name not in cancers_in_which_gene_appears[gene]:
                         cancers_in_which_gene_appears[gene].append(db_name)

                if True:


                    for i in range(len(all_mutated_genes_from_the_list)):
                        gene1 = all_mutated_genes_from_the_list[i]
                        for j in range(i+1,len(all_mutated_genes_from_the_list)):
                            gene2 = all_mutated_genes_from_the_list[j]
                            if gene2 > gene1:
                                mut_key = gene1+"_"+gene2
                            else:
                                mut_key = gene2+"_"+gene1
                            co_appearance[mut_key] += 1

        pancan_samples += number_of_patients[db_name]
         

        for gene in gene_list:
            pancan_ct[gene] += mut_ct [gene]

        for  mut_key in co_appearance.keys():
            pancan_coappearance[mut_key] += co_appearance[mut_key]

        ##########################;
        print >> sys.stderr, "%s done in %8.2f min " % (db_name, (time()-start)/60)
        

    #print "######################################"
    #print "pan-cancer"
    #print "number of samples:", pancan_samples
    print "%%  %8s   %4s   %8s  %4s  %15s  %15s  %15s    %15s  %15s" %  ("gene1", "#muts1", "gene2", "#muts2", 
                                                                "co-appearance", "expected no", "expected no", "-log10(pval <=)", "-log10(pval >=)")
    print " %8s   %4s   %8s  %4s  %15s          %15s  %15s    %15s  %15s" %  ("", "", "", "", "", "of co-app (expr)", "of co-app (sim)", "", "")



    for i in range(len(gene_list)):
        gene1 = gene_list[i]
        for j in range(i+1,len(gene_list)):
            gene2 = gene_list[j]
            if gene2 > gene1:
                mut_key = gene1+"_"+gene2
            else:
                mut_key = gene2+"_"+gene1


            appears_together = pancan_coappearance[mut_key]
            ct1 = pancan_ct[gene1] 
            ct2 = pancan_ct[gene2]

            effective_number_of_samples = 0
            for db_name in db_names:
                if db_name in cancers_in_which_gene_appears[gene1] and db_name in cancers_in_which_gene_appears[gene2] :
                    effective_number_of_samples += number_of_patients[db_name]

            if effective_number_of_samples == 0: continue


            expctd =  expected (ct1, ct2, effective_number_of_samples)
            #print " %8s   %4d   %8s  %4d  %15d  %15.2f    %15.2f " %  ( gene1, ct1, gene2, ct2, appears_together, expctd, float(expctd+1)/(appears_together+1))
            #print "\t\t\t 0.5 < float(expctd+1)/(appears_together+1) < 1.3 ?", 0.5 < float(expctd+1)/(appears_together+1) < 1.3
            if  2.0/3.0 < float(expctd+1)/(appears_together+1) < 3.0/2.0: continue
            number_of_iterations = 5*effective_number_of_samples
            cmd = "meta/ovlp_sim  %d  %d %d  %d  %d " %  (effective_number_of_samples,  ct1,   ct2,  appears_together, number_of_iterations)
            [avg, eval_le, eval_ge] = [float(x) for x in subprocess.check_output(cmd, shell=True).split()]

            if (eval_le <3 and eval_ge <3): continue
            print "  %8s   %4d   %8s  %4d  %15d  %15.2f  %15.2f    %15.4f  %15.4f " %  \
                ( gene1, ct1, gene2, ct2, appears_together, expctd, avg, eval_le, eval_ge)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

