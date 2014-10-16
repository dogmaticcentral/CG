#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);


import sys, os, math
import MySQLdb
import commands
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *
from scipy import stats
from scipy.stats import binom

peptide_length = {}
length_per_category  = {}
silent_length_per_category  = {}
    

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

    qry = "select distinct category from baseline.possible_mutations "
    rows = search_db(cursor, qry)
    if rows:
        for row in rows:
            category = row[0]
            length_per_category[category] = {}
            silent_length_per_category[category] = {}
            
    else:
        print "no categories in baseline.possible_mutations (?)"
        exit (1)


    qry = "select * from baseline.possible_mutations" 
    rows = search_db(cursor, qry)
    if rows:
        # I have the mutations borken down into categories;
        # for the purposes here, sum them all up
        for row in rows:
            [ensembl_id, category, silent, nonsense, missense ] = row

            if not silent_length_per_category[category].has_key(ensembl_id):
                silent_length_per_category[category][ensembl_id] = 0
            silent_length_per_category[category][ensembl_id] += int(silent)

            if not length_per_category[category].has_key(ensembl_id):
                length_per_category[category][ensembl_id] = 0
            length_per_category[category][ensembl_id] += int(silent) + int(nonsense) + int(missense)
    else:
        print "error searching for categorized lengths:"
        print qry
        exit(1)


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
class Failure_counter:
    def __init__ (self):
        self.no_aa_pos_match       = 0
        self.no_cdna_pos_match     = 0
        self.no_aa_from_to_match   = 0
        self.no_cdna_from_to_match = 0
        self.not_a_point_mutation  = 0
        self.not_a_snp            = 0
        self.nt_and_aa_positions_do_not_match  = 0
        self.no_pos_match_to_canonical_ensembl = 0
        self.not_aa_match_to_canonical_ensembl  = 0
        self.not_nt_match_to_canonical_ensembl  = 0

    # print
    def __str__ (self):
        printstr = ""
        for attr, value in self.__dict__.iteritems():
            if ( not value is None):
                printstr += " %-20s    %d" % (attr,  value)
            else:
                printstr += " %-20s    None"  % attr
            printstr += "\n"
        return printstr
    
        

############################################
def parse_cdna_change (cursor, db_name, ensembl_id, variant_classification, aa_change, cdna_change, failure_counter, qry0):


    #print ensembl_id, variant_classification, aa_change, cdna_change

    position_pattern   = re.compile ('[pc]\.\D*(\d+)\D*')
    from_to_pattern    = re.compile ('[pc]\.[\d\>]*(\D*)[\d\>]+(\D*)')

    aa_position_match   = position_pattern.match (aa_change)
    cdna_position_match = position_pattern.match (cdna_change)

    aa_from_to_match   = from_to_pattern.match (aa_change)
    cdna_from_to_match = from_to_pattern.match (cdna_change)

    if not aa_position_match:
        failure_counter.no_aa_pos_match += 1
        return

    if not cdna_position_match:
        failure_counter.no_cdna_pos_match += 1
        return

    if not aa_from_to_match:
        failure_counter.no_aa_from_to_match += 1
        return

    if not cdna_from_to_match:
        failure_counter.no_cdna_from_to_match += 1
        return
    # some curators believe that the actual value of aa or nt is not important
    # (like we have a reliable sequence info, so we need no double checking)
    if not  (len(aa_from_to_match.group(1))<=1  and len(aa_from_to_match.group(2)) <= 1) :
        failure_counter.not_a_point_mutation  +=1 
        #print
        #print "not a point mutation"
        #print ensembl_id, variant_classification, aa_change, cdna_change
        #print "canonical according to ensembl:", ensembl_aa, ensembl_position, ensembl_codon
        return

    if not  (len(cdna_from_to_match.group(1))<=1  and len(cdna_from_to_match.group(2)) <= 1) :
        failure_counter.not_a_snp  +=1 
        #print
        #print "not a snp"
        #print ensembl_id, variant_classification, aa_change, cdna_change
        #print "canonical according to ensembl:", ensembl_aa, ensembl_position, ensembl_codon
        return

    
    position_protein =  int (aa_position_match.group(1))
    position_dna     =  int (cdna_position_match.group(1))
    # are the protein and dna positions related as expected
    pos_computes =   position_protein*3-2 <=  position_dna  <= position_protein*3;

    if not pos_computes:
        failure_counter.nt_and_aa_positions_do_not_match += 1
        return

    aa_from =  aa_from_to_match.group(1)
    aa_to   =  aa_from_to_match.group(2)
    nt_from =  cdna_from_to_match.group(1)
    nt_to   =  cdna_from_to_match.group(2)

    position_pattern   = re.compile ('(\D*)(\d+)(\D*)')
    qry = "select sequence from baseline.canonical_sequence where ensembl_id = '%s'"  %  ensembl_id
    rows = search_db(cursor, qry)
    if not  rows or len(rows)>1:
        print  "oink ?"
        exit(1)

    canonical_sequence = rows[0][0]
    fields = canonical_sequence.split(';')
    ensembl_entry = None
    for field in fields:
        if not field: continue
        position = int(position_pattern.match(field).group(2))
        if position == position_protein:
            ensembl_entry = field
            break

    if not ensembl_entry:
        failure_counter.no_pos_match_to_canonical_ensembl += 1
        return

    ensembl_aa       =  position_pattern.match(ensembl_entry).group(1)
    ensembl_position =  int(position_pattern.match(field).group(2))
    ensembl_codon    =  position_pattern.match(ensembl_entry).group(3)

    #print "cannonical according to ensembl:", ensembl_aa, ensembl_position, ensembl_codon

    if aa_from and ensembl_aa != aa_from:
        #print " !!!! ", qry0
        failure_counter.not_aa_match_to_canonical_ensembl += 1
        return
    # 'index' goes from 0-2, position counts from 1
    tcga_index_in_codon = (position_dna+2)%3
    ensembl_nt = ensembl_codon[tcga_index_in_codon]
    if nt_from and ensembl_nt != nt_from:
        #print " !!!! ", qry0
        failure_counter.not_nt_match_to_canonical_ensembl += 1
        return
       



############################################
def mutations_per_gene (cursor, db_name, genes, samples, mutations_in_sample):

    categories = length_per_category.keys()

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
        found = False

        for category in categories:
            if silent_length_per_category[category].has_key(ensembl_id):
                found = True

        if not found:
            silent_len_not_found += 1

        gene2ensembl_id[gene] = ensembl_id
        if not ensembl2gene.has_key(ensembl_id):
            ensembl2gene[ensembl_id] = []
        ensembl2gene[ensembl_id].append(gene)


    ############################
    total_length        = 0

    total_length_per_category        = {}
    total_silent_length_per_category = {}

    for category in categories:
        total_length_per_category[category]        = 0
        total_silent_length_per_category[category] = 0

    #print "% genes entered under different names:"
    for ensembl_id, genes_w_this_ensembl in ensembl2gene.iteritems():
        #if len(genes_w_this_ensembl) >= 2: 
        #    print "%", ensembl_id, genes_w_this_ensembl
        total_length        += peptide_length[ensembl_id]
        
        for category in categories:
            total_length_per_category[category]         += length_per_category[category][ensembl_id]
            total_silent_length_per_category[category]  += silent_length_per_category[category][ensembl_id]

    print "% number of unresolved gene names:", not_found
    print "% peptide_length_not_found:", peplen_not_found
    print "% silent_length_not_found:", silent_len_not_found
    print "% total_length", total_length
    for category in categories:
        print "%%\t total_length  in category %s:  %d   silent: %d"  % (category,  total_length_per_category[category] ,  total_silent_length_per_category[category])

    ############################
    # total and silent number of mutations
    number_of_codon_mutations_per_gene = {}
    number_of_silent_mutations_per_gene = {}
    for category in categories:
        number_of_codon_mutations_per_gene[category] = {}
        number_of_silent_mutations_per_gene[category]  = {}
        

    for category in categories:
        for ensembl_id in ensembl2gene.keys():
            number_of_codon_mutations_per_gene[category][ensembl_id]  = 0
            number_of_silent_mutations_per_gene[category][ensembl_id] = 0

    total_number_of_codon_mutations   = {};
    total_number_of_silent_mutations  = {};

    for category in categories:
        total_number_of_codon_mutations[category]   = 0;
        total_number_of_silent_mutations[category]  = 0;

    switch_to_db(cursor, db_name)

    failure_counter = Failure_counter()
    ct = 0
    for ensembl_id, genes_w_this_ensembl in ensembl2gene.iteritems():
        ct += 1
        for gene in genes_w_this_ensembl:

            qry = "select aa_change, cdna_change, variant_classification  from somatic_mutations where hugo_symbol='%s'" % gene
            rows = search_db(cursor, qry)
            if not rows: continue

            for row in rows:
                variant_classification = row[2].lower()
                if  not variant_classification in ['silent', 'missense_mutation','nonsense_mutation']: continue

                aa_change   = row[0]
                cdna_change = row[1]

                parse_cdna_change (cursor, db_name, ensembl_id, variant_classification, aa_change, cdna_change, failure_counter, qry)
        #if  ct == 1000: break
    
    print failure_counter
    return
  
    overall_mutation_rate = float(total_number_of_codon_mutations)/total_length;
    overall_silent_rate   = float(total_number_of_silent_mutations)/total_silent_length;
   
    print  "%% total_number_of_codon_mutations %d,  mutation rate: %5.2e,  silent rate: %5.2e" %  \
        (total_number_of_codon_mutations, overall_mutation_rate, overall_silent_rate)

    if not total_number_of_codon_mutations: return

    print "% correlation between expected and observed number of mutations, "
    print "% taking the rates to be uniform accross genes and across samples "
    print "% (that is one mutation rate fo each cancer type) "
    print "% interestingly enough, focusing on silent mutations only makes things worse"

    x = []
    y = []
    for ensembl_id in ensembl2gene.keys():
        expected_number_of_mutations = 3*peptide_length[ensembl_id]*overall_mutation_rate
        observed = number_of_codon_mutations_per_gene[ensembl_id]
        x.append(expected_number_of_mutations)
        y.append(observed)

    [pearson_corr, p_val] = stats.pearsonr(x, y)
    print "% pearson_corr for  overall_mutation_rate", pearson_corr
       
    x = []
    y = []
    for ensembl_id in ensembl2gene.keys():
        expected_number_of_mutations = silent_length[ensembl_id]*overall_silent_rate
        observed = number_of_silent_mutations_per_gene[ensembl_id]
        x.append(expected_number_of_mutations)
        y.append(observed)

    [pearson_corr, p_val] = stats.pearsonr(x, y)
    print "% pearson_corr for  silent_mutation_rate", pearson_corr
     

    

    

#########################################
def main():
    
    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = [ "ACC", "BLCA", "BRCA", "CESC", "COAD",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]

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
