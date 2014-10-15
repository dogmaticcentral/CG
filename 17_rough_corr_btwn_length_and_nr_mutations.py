#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);


import sys, os, math
import MySQLdb
import commands
from   tcga_utils.mysql   import  *
from scipy import stats
from scipy.stats import binom

#########################################
def retrieve_ensembl_id (cursor, gene):
    ensembl_gene_id = ""

    if gene=='Unknown':
        qry = "select entrez_gene_id from somatic_mutations  where hugo_symbol='%s'" % gene
        rows = search_db(cursor, qry)
        if not rows:   return ""
        entrez_gene_id = rows[0][0]
        qry = "select ensembl_gene_id from baseline.hgnc_id_translation where entrez_gene_id='%s'" % str(entrez_gene_id)
        rows = search_db(cursor, qry)
        if not rows:   return ""
        ensembl_gene_id = rows[0][0]
    else: 
        qry = "select ensembl_gene_id, comment from name_resolution where tcga_hugo_symbol='%s'" % gene
        rows = search_db(cursor, qry)
        if not rows or rows[0][1]=='failure':   return ""
        ensembl_gene_id = rows[0][0]

    return ensembl_gene_id 
      
#########################################
# what is the length of the peptide this gene translates to?
def get_peptide_length (cursor, gene):
    peptide_length = 0
    ensembl_id = retrieve_ensembl_id (cursor, gene)
    if not ensembl_id:
        #print "ensembl id not found for ", gene
        return peptide_length

    qry = "select peptide from baseline.aa_sequence where ensembl_id='%s'" % ensembl_id
    rows = search_db(cursor, qry)
    if not rows or not rows[0][0]:
        #print 'no peptide found for ', gene
        return peptide_length

    peptide = rows[0][0]
    peptide_length = len(peptide)
    return peptide_length
            
#########################################
def mutations_per_gene (cursor, db_name):
    table = 'somatic_mutations'

    switch_to_db (cursor, db_name)

    ############################
    print "%% ", commands.getoutput('grep ' + db_name +' cancer_names.txt')
    qry = "select count(1) from " + table
    rows = search_db(cursor, qry)
    no_entries = rows[0][0]
    print "%% number of entries:    %7d" %  no_entries
    ############################

    if not no_entries: return

    ############################
    uniq_patients = {}
    qry  = "select tumor_sample_barcode from somatic_mutations "
    rows = search_db(cursor, qry)
    for  row in rows:
        tbarcode = row[0]
        fields = tbarcode.split('-')
        patient = '-'.join(fields[1:3])
        if not  uniq_patients.has_key(patient):
            uniq_patients[patient] = []
        uniq_patients[patient].append('-'.join(fields[3:]))
    print "%% number of patients:   %7d   (mutations per sample %8.2f)" % (len(uniq_patients), float(no_entries)/len(uniq_patients))

    ############################
    qry = "select distinct(hugo_symbol) from somatic_mutations"
    rows = search_db(cursor, qry)
    genes = [row[0] for row in  rows if (row[0] != 'abParts' and row[0] != 'Unknown')]
    print "%% number of  genes:     %7d" %  len(genes)

    ############################
    qry = "select distinct(tumor_sample_barcode) from somatic_mutations"
    rows = search_db(cursor, qry)
    samples = [row[0] for row in  rows]
    print "%% number of  samples:   %7d" %  len(samples)
    mutations_in_sample = {}
    for sample in samples:
        qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" %  sample
        rows = search_db(cursor, qry)
        if not rows: 
            mutations_in_sample[sample] = 0
        else:
            mutations_in_sample[sample] = rows[0][0]
 
    ############################
    entries_per_gene = {}
    indel_per_gene   = {}
    codon_mutations_per_gene   = {}
    codon_mutations_per_sample = {}
    codon_mutations_per_gene_per_sample = {}
    total_length_per_sample = {}
    silent_per_gene  = {}
    peptide_length   = {}
    total_length     = 0
    total_hits = 0 
    gene2sample = {}
    for gene in genes:
        switch_to_db (cursor, db_name)
        qry = "select variant_classification, variant_type, tumor_sample_barcode  from somatic_mutations where hugo_symbol='%s'" % gene
        rows = search_db(cursor, qry)
        entries_per_gene[gene] = len(rows)
        peptide_length[gene]   = get_peptide_length(cursor, gene)

        indel_per_gene[gene]  = 0
        silent_per_gene[gene] = 0
        codon_mutations_per_gene[gene] = 0
        codon_mutations_per_gene_per_sample[gene] = {}
        gene2sample[gene] = []
        samples_seen = {}
        for row in rows:
            sample = row[2]
            if  not sample in codon_mutations_per_sample.keys():
                codon_mutations_per_sample[sample]  = 0
 
            if  not sample in codon_mutations_per_gene_per_sample[gene].keys():
                codon_mutations_per_gene_per_sample[gene][sample] = 0

            gene2sample[gene].append(sample) 
            [variant_classification, variant_type] = [x.lower() for x in row[:2]]
            if variant_type=='ins' or variant_type=='del':
                indel_per_gene[gene] +=1
            elif  variant_classification in ['silent', 'missense_mutation','nonsense_mutation']:
                codon_mutations_per_gene[gene]     += 1
                codon_mutations_per_sample[sample] += 1
                codon_mutations_per_gene_per_sample[gene][sample] += 1
                if  variant_classification =='silent':
                    silent_per_gene[gene]  += 1

            if not sample in samples_seen.keys():
                samples_seen[sample] = True
                if  not sample in total_length_per_sample.keys():
                    total_length_per_sample[sample]  = 0
                total_length_per_sample[sample] += 3*peptide_length[gene] 

        total_length += 3*peptide_length[gene]
        total_hits   += codon_mutations_per_gene[gene]
        
    number_of_samples = len (codon_mutations_per_sample.keys())
    uniform_bg_mutation_rate    = float(total_hits)/total_length
    bg_mutation_rate_per_sample = {}
    for sample in codon_mutations_per_sample.keys():
        if not total_length_per_sample[sample]: 
            print ">>>",  sample, total_length_per_sample[sample]
            continue
        bg_mutation_rate_per_sample[sample] = float(codon_mutations_per_sample[sample])/total_length
        #print "\t %s  %8d  %8d   %8d  %5.3e" % ( sample, codon_mutations_per_sample[sample], mutations_in_sample[sample],
        #                                total_length_per_sample[sample], bg_mutation_rate_per_sample[sample])

    
    print  "total_hits,  total_length, uniform_bg_mutation_rate:  %8d  %8d  %5.3e" %\
        ( total_hits,  total_length, uniform_bg_mutation_rate)

    if not total_hits: return

    model1 = []
    model2 = []
    number_of_hits_1 = []
    number_of_hits_2 = []

    for gene in genes:
        pep_len = peptide_length[gene]
        if not pep_len: continue
        model1.append (3*pep_len*uniform_bg_mutation_rate)
        number_of_hits_1.append (codon_mutations_per_gene[gene])
        # print " ** ", 3*pep_len*uniform_bg_mutation_rate, codon_mutations_per_gene[gene]
        # for model2 find expected number fo mutations in gene in sample, and sum them over all samples
        expected_total_number_of_mutations = 0
        for sample, mutation_rate in bg_mutation_rate_per_sample.iteritems():
            if  sample in codon_mutations_per_gene_per_sample[gene].keys():
                expected_number_of_mutations = 3*mutation_rate*pep_len;
                model2.append (expected_number_of_mutations)
                number_of_hits_2.append (codon_mutations_per_gene_per_sample[gene][sample])

    [pearson_corr, p_val] = stats.pearsonr(model1, number_of_hits_1)
    print "pearson_corr, model1", pearson_corr
    [pearson_corr, p_val] = stats.pearsonr(model2, number_of_hits_2)
    print "pearson_corr, model2", pearson_corr

    return

    avg = {}
    avg_sq = {}
    count = {}
    for gene in sorted_genes:
        if  mutations_in_sample[gene2sample[gene]] > 50: continue
        if  not entries_per_gene.has_key(gene): continue
        pep_len =  peptide_length[gene]
        if not pep_len: continue
        pep_bin = pep_len/20
        if not avg.has_key(pep_bin):
            avg[pep_bin]    = 0
            count[pep_bin]  = 0
            avg_sq[pep_bin] = 0
        avg[pep_bin]    += entries_per_gene[gene]
        avg_sq[pep_bin] += entries_per_gene[gene]*entries_per_gene[gene]
        count[pep_bin] += 1

    for pep_bin in avg.keys():
        avg_no_muts =  float(avg[pep_bin])/count[pep_bin]
        stdev       =  float(avg_sq[pep_bin])/count[pep_bin]
        stdev       =  math.sqrt(stdev -avg_no_muts*avg_no_muts)
        print "  %8d   %8.2f   %8.2f  " %  (pep_bin*20, avg_no_muts, stdev)


#########################################
def main():
    
    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = [ "ACC", "BLCA", "BRCA", "CESC", "COAD",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]
    for db_name  in db_names:
        print '% ==================================='
        mutations_per_gene (cursor, db_name)
        print

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
