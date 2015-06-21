#!/usr/bin/python -u
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *
import matplotlib.pyplot as plt

#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

    full_name = read_cancer_names ()

    for db_name in db_names:
        print "######################################"
        print db_name
        switch_to_db (cursor, db_name)

        ############################
        print db_name, "number of somatic_mutations:"
        qry = "select count(1)from somatic_mutations"
        rows = search_db(cursor, qry)
        print "\t", rows[0][0]
        print
        ############################

        ############################
        print "number of different genes:"
        qry = "select distinct(hugo_symbol) from somatic_mutations"
        rows  = search_db(cursor, qry)
        genes = [row[0] for row in  rows]
        print "\t", len(genes)

        ############################
        print "number of different patients:"
        qry = "select distinct(sample_barcode_short) from somatic_mutations"
        rows = search_db(cursor, qry)
        patients = [row[0] for row in  rows]
        total_patients = len(patients)
        print "\t", total_patients

        ############################
        print "frequencies reported per gene"
        freq_gene = {}
        for gene in genes:
            qry  = "select sample_barcode_short, count(sample_barcode_short) from somatic_mutations "
            qry += "where hugo_symbol='%s' " % gene
            qry += "and variant_classification!='silent' and variant_classification!='RNA' "
            qry += "group by sample_barcode_short"
            rows = search_db(cursor, qry)
            if rows:
                no_patients = len(rows)
            else:
                no_patients = 0
            if no_patients==0: continue
            if no_patients>total_patients/5:
                print " %10s   %4d   %6.2f%%" % (gene, no_patients, float(no_patients)/total_patients)
            freq_gene[gene] = float(no_patients)/total_patients*100

        sorted_genes =  sorted(freq_gene.keys(), key= lambda x: -freq_gene[x])

        fig, ax1 = plt.subplots(1, 1)
        title = full_name[db_name]
        ax1.set_title (title, fontsize=20)
        ax1.set_xlabel('genes, listed by their rank', fontsize = 24)
        ax1.set_ylabel('% of patients', fontsize = 24)
        if 'RPL5' in sorted_genes:
            less_mutated = len( [y for y in  freq_gene.values() if y<freq_gene['RPL5']])
            more_mutated = len( [y for y in  freq_gene.values() if y>freq_gene['RPL5']])
            rank_msg = "RPL5 rank:   %d-%d  (%.1f%% patients)" % (more_mutated, len(freq_gene)-less_mutated, freq_gene['RPL5'])
        else:
            rank_msg = "RPL5 rank:   %d  (no patients)" %  len(sorted_genes)
        ax1.set_axis_bgcolor((0, 102./255, 204./255))
        x = range(1,len(sorted_genes)+1)
        y = [freq_gene[gene] for gene in sorted_genes]
        ax1.plot (x, y, 'r-', lw=1, alpha=0.6, label=rank_msg)
        ax1.fill_between(x, y,  interpolate=True, color=(255./255,153./255,51./255))
        ax1.legend()
        plt.ylim(0,min(max(y),10))
        #plt.show()
        filename = db_name+"_somatic_freq.png"
        plt.savefig(filename)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

