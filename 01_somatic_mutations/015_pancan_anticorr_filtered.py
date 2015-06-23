#!/usr/bin/python -u
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *

drop_silent = True


#########################################
def silent_proportion(cursor, db_names, gene):
    non_silent_ct = 0
    silent_ct     = 0
    for db_name in db_names:
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

    return [silent_ct, non_silent_ct]


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["STAD","ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM",  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]


    infile = open ("coapp_tables/pancan_tp53_coapps.table", "r")
    for line in infile.readlines()[3:]:
        spl = line.split()
        if len(spl) < 9: continue
        [gene1, no_muts1, gene2, no_muts2, coapps, exp_coapp_analytic, exp_coapp_sim, pval_less, pval_gt] = spl
        pval_less = float(pval_less)
        if pval_less > 0.1: continue
        [silent_ct, non_silent_ct] = silent_proportion (cursor, db_names, gene2)
        if not non_silent_ct: continue
        if float(silent_ct)/non_silent_ct > 0.15: continue
        print " %5s  %4d  %12s  %4d    %4d   %7.2f    %.4f " % (gene1, int(no_muts1), gene2, int(no_muts2),
                                                           int(coapps), float(exp_coapp_analytic), float(pval_less))


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

