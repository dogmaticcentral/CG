#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);


import sys, os, math
import MySQLdb
import commands
import random
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *
from scipy import stats
from scipy.stats import binom

peptide_length    = {}
silent_length     = {}
non_silent_length = {}
    


#########################################
# what is the length of the peptide this gene translates to?
def set_peptide_lengths (cursor):

    position_pattern   = re.compile ('\D*(\d+)\D*')
    qry = "select ensembl_id, sequence from baseline.canonical_sequence" 
    rows = search_db(cursor, qry)
    if rows:
        for row in rows:
            [ensembl_id, canonical_sequence] = row
            fields = canonical_sequence.split(';')
            last_field    = fields[-2]
            last_position = position_pattern.match(last_field).group(1)
            peptide_length[ensembl_id] =  int(last_position)
    else:
        print "error searching for protein lengths:"
        print qry

#########################################
#
def check_number_of_mutations (cursor):
    qry = "select ensembl_id, silent, nonsense,  missense  from baseline.possible_mutations" 
    rows = search_db(cursor, qry)
    tot_length = {}
    if rows:
        # I have the mutations borken down into categories;
        # for the purposes here, sum them all up
        for row in rows:
            [ensembl_id,  silent, nonsense,  missense] = row
            if not silent_length.has_key(ensembl_id):
                tot_length[ensembl_id] = 0
            tot_length[ensembl_id] += int(silent) +  int(nonsense)  + int(missense)
    for ensembl_id, length in silent_length.iteritems():
        if  3*3*peptide_length[ensembl_id] != length:
            print ensembl_id,  3*3*peptide_length[ensembl_id],length
    exit(1)


#########################################
# what is the effective 'silent' length of the peptide = the number of silent mutations possible
def set_silent_mutation_lengths (cursor):
    qry = "select ensembl_id, silent, nonsense, missense from baseline.possible_mutations" 
    rows = search_db(cursor, qry)
    if rows:
        # I have the mutations borken down into categories;
        # for the purposes here, sum them all up
        for row in rows:
            [ensembl_id, number_of_silent_mutations, nonsense, missense] = row
            if not silent_length.has_key(ensembl_id):
                silent_length[ensembl_id] = 0
                non_silent_length[ensembl_id] = 0
            silent_length[ensembl_id] += int(number_of_silent_mutations)
            non_silent_length[ensembl_id] += int(nonsense) + int (missense)
    else:
        print "error searching for 'silent' lengths:"
        print qry


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
    else: 
        qry = "select ensembl_gene_id, comment from name_resolution where tcga_hugo_symbol='%s'" % gene
        rows = search_db(cursor, qry)
        if not rows or rows[0][1]=='failure':   return ""
        
    if rows and 'ENSG' in rows[0][0]: 
        return rows[0][0]
    else:
        return ""

            
#########################################
def overall_stats  (cursor, db_name):
    table = 'somatic_mutations'

    switch_to_db (cursor, db_name)

    ############################
    print "%% ", commands.getoutput('grep ' + db_name +' cancer_names.txt')
    qry = "select count(1) from " + table
    rows = search_db(cursor, qry)
    no_entries = rows[0][0]
    print "%% number of entries:    %7d" %  no_entries
    ############################

    if not no_entries: return []

    ############################
    uniq_patients = {}
    qry  = "select distinct tumor_sample_barcode from somatic_mutations "
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

    ############################
    mutations_in_sample = {}
    for sample in samples:
        qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" %  sample
        rows = search_db(cursor, qry)
        if not rows: 
            mutations_in_sample[sample] = 0
        else:
            mutations_in_sample[sample] = rows[0][0]

    return [genes, samples, mutations_in_sample]



############################################
def mutations_per_gene (cursor, db_name, genes, samples, mutations_in_sample):

    # qry = "select variant_classification, variant_type, tumor_sample_barcode from somatic_mutations where hugo_symbol='%s'" % gene
 
    ############################
    gene2ensembl_id     = {}
    ensembl2gene        = {}

    not_found           = 0
    peplen_not_found    = 0
    silent_len_not_found = 0

    for gene in genes:
        ensembl_id = hugo2ensembl(cursor, gene) # in utils
        if not ensembl_id:
            ensembl_id  = retrieve_ensembl_id (cursor, gene)
            if not ensembl_id:
                not_found += 1
                continue
            else:
                print  gene, "resolved in name_res_table to ", ensembl_id
        if not peptide_length.has_key(ensembl_id):
            peplen_not_found += 1
            continue
        if not silent_length.has_key(ensembl_id):
            silent_len_not_found += 1

        gene2ensembl_id[gene] = ensembl_id
        if not ensembl2gene.has_key(ensembl_id):
            ensembl2gene[ensembl_id] = []
        ensembl2gene[ensembl_id].append(gene)

    ############################
    total_length        = 0
    total_silent_length = 0
    total_non_silent_length = 0
    print "% genes entered under different names:"
    for ensembl_id, genes_w_this_ensembl in ensembl2gene.iteritems():
        #if len(genes_w_this_ensembl) >= 2: 
        #    print "%", ensembl_id, genes_w_this_ensembl
        total_silent_length     += silent_length[ensembl_id]
        total_non_silent_length += non_silent_length[ensembl_id]
    total_length = total_silent_length + total_non_silent_length
       

    print "% number of unresolved gene names:", not_found
    print "% peptide_length_not_found:", peplen_not_found
    print "% silent_length_not_found:", silent_len_not_found
    print "% total_length, total_silent_length, total_silent_length: ", total_length, total_silent_length, total_non_silent_length

    ############################
    # total and silent number of mutations
    number_of_codon_mutations_per_gene = {}
    number_of_silent_mutations_per_gene = {}
    number_of_non_silent_mutations_per_gene = {}
    for ensembl_id in ensembl2gene.keys():
        number_of_codon_mutations_per_gene[ensembl_id]  = 0
        number_of_silent_mutations_per_gene[ensembl_id] = 0
        number_of_non_silent_mutations_per_gene[ensembl_id] = 0
    total_number_of_codon_mutations  = 0;
    total_number_of_silent_mutations  = 0;
    total_number_of_non_silent_mutations  = 0;

    switch_to_db(cursor, db_name)

    for ensembl_id, genes_w_this_ensembl in ensembl2gene.iteritems():
        for gene in genes_w_this_ensembl:

            qry = "select tumor_sample_barcode, variant_classification  from somatic_mutations where hugo_symbol='%s'" % gene
            rows = search_db(cursor, qry)
            if not rows: continue

            for row in rows:
                tumor_sample_barcode = row[0]
                if mutations_in_sample[tumor_sample_barcode] < 90: continue
                variant_classification = row[1].lower()
                if  variant_classification in ['silent', 'missense_mutation','nonsense_mutation']:
                    number_of_codon_mutations_per_gene[ensembl_id] += 1
                    total_number_of_codon_mutations += 1
                    if  variant_classification == 'silent':
                        number_of_silent_mutations_per_gene[ensembl_id] += 1
                        total_number_of_silent_mutations += 1
                    else:
                        number_of_non_silent_mutations_per_gene[ensembl_id] += 1
                        total_number_of_non_silent_mutations += 1
                       

    overall_mutation_rate    = float(total_number_of_codon_mutations)/total_length;
    overall_silent_rate      = float(total_number_of_silent_mutations)/total_silent_length;
    overall_non_silent_rate  = float(total_number_of_non_silent_mutations)/total_non_silent_length;
   
    print  "%% total_number_of_codon_mutations %d,  mutation rate: %5.2e,  silent rate: %5.2e,  non_silent rate: %5.2e" %  \
        (total_number_of_codon_mutations, overall_mutation_rate, overall_silent_rate, overall_non_silent_rate)

    if not total_number_of_codon_mutations: return

    print "% correlation between expected and observed number of mutations, "
    print "% taking the rates to be uniform accross genes and across samples "
    print "% (that is one mutation rate fo each cancer type) "
    print "% interestingly enough, focusing on silent mutations only makes things worse"

    ################################################
    x = []
    y = []
    for ensembl_id in ensembl2gene.keys():
        expected_number_of_mutations = (silent_length[ensembl_id] + non_silent_length[ensembl_id])*overall_mutation_rate
        observed =  number_of_silent_mutations_per_gene[ensembl_id] + number_of_non_silent_mutations_per_gene[ensembl_id]
        x.append(expected_number_of_mutations)
        y.append(observed)

    [size, [hi,lo],  mean, variance, skew, kurtosis] = stats.describe(y)
    print "%% statistics for y:   size %d     mean %6.1f    variance:  %6.1f  (sqrt:  %6.1f)" % (size, mean, variance, math.sqrt(variance))
    [pearson_corr, p_val] = stats.pearsonr(x, y)
    print "%% pearson_corr for  codon_mutation_rate :  %6.2e,   pval  %6.1e" % (pearson_corr, p_val)


    ################################################
    
    print " silent mutations only "

    ################################################
    x = []
    y = []
    for ensembl_id in ensembl2gene.keys():
        expected_number_of_mutations = silent_length[ensembl_id]*overall_silent_rate
        observed = number_of_silent_mutations_per_gene[ensembl_id]
        x.append(expected_number_of_mutations)
        y.append(observed)

    [size, [hi,lo],  mean, variance, skew, kurtosis] = stats.describe(y)
    print "%% statistics for y:   size %d     mean %6.1f    variance:  %6.1f  (sqrt:  %6.1f)" % (size, mean, variance, math.sqrt(variance))
    [pearson_corr, p_val] = stats.pearsonr(x, y)
    print "%% pearson_corr for  silent_mutation_rate :  %6.2e,   pval  %6.1e" % (pearson_corr, p_val)

    ################################################
    
    print " non-silent mutations only "

    ################################################
    x = []
    y = []
    for ensembl_id in ensembl2gene.keys():
        expected_number_of_mutations = non_silent_length[ensembl_id]*overall_non_silent_rate
        observed = number_of_non_silent_mutations_per_gene[ensembl_id]
        x.append(expected_number_of_mutations)
        y.append(observed)

    [size, [hi,lo],  mean, variance, skew, kurtosis] = stats.describe(y)
    print "%% statistics for y:   size %d     mean %6.1f    variance:  %6.1f  (sqrt:  %6.1f)" % (size, mean, variance, math.sqrt(variance))
    [pearson_corr, p_val] = stats.pearsonr(x, y)
    print "%% pearson_corr for  non_silent_mutation_rate :  %6.2e,   pval  %6.1e" % (pearson_corr, p_val)


#########################################
def main():
    
    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = [ "ACC", "BLCA", "BRCA", "CESC", "COAD",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]

    #db_names = ["COAD"]

    # cancer that look like they are driven by radnom point  mutations
    # large number of patients, large number of mutations per gene, large chunk of the genome covered
    # as a byproduct: higher rate of mutations in 
    # even thoug I think I'll have to drop BRCA bcs they do not give me the precise mutation location
    db_names = ["BRCA", "COAD", "HNSC", "KIRC", "LIHC", "LUAD", "REA", "SKCM", "STAD"]
    
    set_peptide_lengths(cursor) # fills peptide_length dictionary
    #check_number_of_mutations(cursor) $ looks like we pass this test
    set_silent_mutation_lengths(cursor) # fills peptide_length dictionary
   
    for db_name  in db_names:
        print '% ==================================='
        ret = overall_stats (cursor, db_name)
        if not len(ret) == 3: continue
        [genes, samples, mutations_in_sample] = ret
        mutations_per_gene (cursor, db_name, genes, samples, mutations_in_sample )
        print

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
