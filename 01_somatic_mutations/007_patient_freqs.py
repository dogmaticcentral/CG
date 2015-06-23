#!/usr/bin/python -u
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *
import matplotlib.pyplot as plt

drop_silent = True

#########################################
def rank_message ( gene_name, freq_gene):
    rank_msg = ""
    if gene_name in freq_gene.keys():
        less_mutated = len( [y for y in  freq_gene.values() if y<freq_gene[gene_name]])
        more_mutated = len( [y for y in  freq_gene.values() if y>freq_gene[gene_name]])
        rank_msg = "%7s rank:   %4d-%4d  (%.1f%% patients)" % (gene_name, more_mutated, len(freq_gene)-less_mutated, freq_gene[gene_name])
    else:
        rank_msg = "%7s rank:   %4d  (no patients)" %  (gene_name, len(freq_gene))
    return rank_msg

#########################################
def  live_plot ( title, freq_gene, filename):

    sorted_genes =  sorted(freq_gene.keys(), key= lambda x: -freq_gene[x])
    fig, ax1 = plt.subplots(1, 1)
    ax1.set_title (title, fontsize=20)
    ax1.set_xlabel('genes, listed by their rank', fontsize = 24)
    ax1.set_ylabel('% of patients', fontsize = 24)
    rpl5_rank_msg  = rank_message('RPL5', freq_gene)
    rpl11_rank_msg = rank_message('RPL11', freq_gene)
    bg_color = (0, 102./255, 204./255)
    ax1.set_axis_bgcolor(bg_color)
    x = range(1,len(sorted_genes)+1)
    y = [freq_gene[gene] for gene in sorted_genes]
    ax1.plot (x, y, 'b-', lw=1, alpha=0.6, label=rpl5_rank_msg)
    ax1.plot (x, y, 'r-', lw=1, alpha=0.6, label=rpl11_rank_msg)
    ax1.fill_between(x, y,  interpolate=True, color=(255./255,153./255,51./255))
    ax1.legend()
    plt.ylim(0,min(max(y),10))

    if filename:
        plt.savefig(filename)
    else:
        plt.show()

#########################################
def silent_proportion(cursor, gene):


    qry  = "select count(1) from somatic_mutations "
    qry += "where hugo_symbol='%s' " % gene
    qry += "and variant_classification in ('Missense_Mutation', 'Nonstop_Mutation', 'Nonsense_Mutation')"
    rows = search_db(cursor, qry)
    if not rows:
        non_silent_ct = 0
    else:
        non_silent_ct = rows[0][0]

    qry  = "select count(1) from somatic_mutations "
    qry += "where hugo_symbol='%s' " % gene
    qry += "and variant_classification='silent'"
    rows = search_db(cursor, qry)
    if not rows:
        silent_ct = 0
    else:
        silent_ct = rows[0][0]

    return [silent_ct, non_silent_ct]


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["STAD","ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM",  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
    #db_names = ["ACC", "GBM", "UVM"]

    full_name = read_cancer_names ()
    pancan_freq = {}
    pancan_silent = {}
    pancan_non_silent = {}
    grand_total_patients = 0
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

        grand_total_patients += total_patients

        ############################
        print "frequencies reported per gene"
        freq_gene = {}
        for gene in genes:

            if drop_silent:
                [silent_ct, non_silent_ct] = silent_proportion(cursor, gene)
                if non_silent_ct==0 or float(silent_ct)/non_silent_ct>0.15: continue

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

            if no_patients>total_patients/10:
                print " %10s   %4d   %6d%%   %3d  %3d   %5.2f" % (gene, no_patients, float(no_patients)/total_patients*100,
                                                      silent_ct, non_silent_ct, float(silent_ct)/non_silent_ct)

            if not pancan_freq.has_key(gene):
                pancan_freq[gene]   = 0
                pancan_silent[gene] = 0
                pancan_non_silent[gene] = 0
            pancan_freq[gene]       += no_patients
            pancan_silent[gene]     += silent_ct
            pancan_non_silent[gene] += non_silent_ct
            freq_gene[gene] = float(no_patients)/total_patients*100

        ###################################
        # in individual tumor types:
        if False:
            filename = db_name+"_somatic_freq.png"
            title = full_name[db_name]
            live_plot (title, freq_gene, "")

    for gene in pancan_freq.keys():
        pancan_freq[gene] /= float(grand_total_patients)
        pancan_freq[gene] *= 100

    sorted_genes =  sorted(pancan_freq.keys(), key= lambda x: -pancan_freq[x])
    for gene in sorted_genes[:100]:
        print " %10s   %6d%%    %5.2f" % (gene,  pancan_freq[gene], float(pancan_silent[gene])/pancan_non_silent[gene])

    filename = "pancan_somatic_freq.filtered_for_silent_proportion.png"
    live_plot ("Pan-cancer statistics", pancan_freq, filename)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

