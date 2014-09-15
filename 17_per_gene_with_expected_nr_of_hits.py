#!/usr/bin/python
# needed the index on hugoSymbol for this to work with any speed:
# create index hugo_idx on somatic_mutations (hugoSymbol);




import sys, os
import MySQLdb
import commands
from   tcga_utils.mysql   import  *

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
    print commands.getoutput('grep ' + db_name +' cancer_names.txt')
    qry = "select count(1) from " + table
    rows = search_db(cursor, qry)
    print "number of entries:    %4d" % rows[0][0]
    ############################

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
    print "number of patients:   %5d" % len(uniq_patients)

    ############################
    qry = "select distinct(hugo_symbol) from somatic_mutations"
    rows = search_db(cursor, qry)
    genes = [row[0] for row in  rows if (row[0] != 'abParts' and row[0] != 'Unknown')]
    print "number of  genes:     %5d" %  len(genes)

    ############################
    qry = "select distinct(tumor_sample_barcode) from somatic_mutations"
    rows = search_db(cursor, qry)
    samples = [row[0] for row in  rows]
    print "number of  samples:     %5d" %  len(samples)
    mutations_in_sample = {}
    for sample in samples:
        qry = "select count(1) from somatic_mutations where tumor_sample_barcode = '%s'" %  sample
        rows = search_db(cursor, qry)
        if not rows: 
            mutations_in_sample[sample] = 0
        else:
            mutations_in_sample[sample] = rows[0][0]
 
    ############################
    print "mutations reported per gene"
    entries_per_gene = {}
    indel_per_gene   = {}
    codon_mutations_per_gene     = {}
    silent_per_gene  = {}
    peptide_length   = {}
    ad_hoc_score = {}
    total_length = 0
    total_hits = 0 
    for gene in genes:
        switch_to_db (cursor, db_name)
        qry = "select variant_classification, variant_type, tumor_sample_barcode  from somatic_mutations where hugo_symbol='%s'" % gene
        rows = search_db(cursor, qry)
        entries_per_gene[gene] = len(rows)

        indel_per_gene[gene]  = 0
        silent_per_gene[gene] = 0
        codon_mutations_per_gene[gene] = 0
        ad_hoc_score[gene] = 0
        for row in rows:
            sample = row[2]
            [variant_classification, variant_type] = [x.lower() for x in row[:2]]
            if variant_type=='ins' or variant_type=='del':
                indel_per_gene[gene] +=1
                if mutations_in_sample[sample]:  ad_hoc_score[gene] += float(100)/mutations_in_sample[sample]
            elif  variant_classification in ['silent', 'missense_mutation','nonsense_mutation']:
                codon_mutations_per_gene[gene]  += 1
                if  variant_classification =='silent':
                    silent_per_gene[gene]  += 1
                else:
                    if mutations_in_sample[sample]:  ad_hoc_score[gene] += float(100)/mutations_in_sample[sample]


        peptide_length[gene] = get_peptide_length(cursor, gene)
        total_length += peptide_length[gene]
        total_hits   += entries_per_gene[gene]

    sorted_genes =  sorted(genes, key= lambda x: entries_per_gene[x])


    print "  %4s  %20s  %5s   %5s    %5s   %8s    %5s    %8s  %5s  %8s" % ( 'rank', 'name', 'tot',
                                                                     'indels',   'codon',                 
                                                                      'silent',  'pep len', 'expected',  'hits/expected', 'score')

    fract_expected = {}
    genes_filtered = []
    for gene in sorted_genes:
        #if gene[:3] != 'RPL' and gene[:3] != 'RPs' : continue
        #iif gene[:3] != 'WNT': continue
        if  not entries_per_gene.has_key(gene): continue
        if codon_mutations_per_gene[gene]:
            ratio = float(silent_per_gene[gene])/codon_mutations_per_gene[gene]
        else:
            ratio = 0
        if  0.25 < ratio  < 0.50: continue
        pep_len =  peptide_length[gene]
        if not pep_len: continue
        expected = float(pep_len)/total_length*total_hits
        if expected:
            fract_expected[gene] = float(entries_per_gene[gene])/expected
        else:
            fract_expected[gene] = 0.0
        if 0.3 < fract_expected[gene] < 3.0: continue
        genes_filtered.append(gene)

    ct = 0
    #max_rank     = len(genes_filtered)
    max_rank     =100
    #sorted_genes =  sorted(genes_filtered, key= lambda x:  (fract_expected[x]*fract_expected[x] + ad_hoc_score[x]*ad_hoc_score[x] ) )
    sorted_genes =  sorted(genes_filtered, key= lambda x:   ad_hoc_score[x] )
    for gene in sorted_genes[-100:]:
        if codon_mutations_per_gene[gene]:
            ratio = float(silent_per_gene[gene])/codon_mutations_per_gene[gene]
        else:
            ratio = 0
        pep_len =  peptide_length[gene]
        if not pep_len: continue
        expected = float(pep_len)/total_length*total_hits
      
        print "  %4d  %20s  %5d   %5d    %5d   %8.2f    %5d    %8.2f       %10.2f     %10.2f  " % ( max_rank-ct, gene, entries_per_gene[gene],
                                                                                       indel_per_gene[gene],    codon_mutations_per_gene[gene], 
                                                                                       ratio,  pep_len, expected,  fract_expected[gene], ad_hoc_score[gene])
        ct += 1
    print

#########################################
def main():
    
    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = [ "ACC", "BLCA", "BRCA", "CESC", "COAD",  "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA", # READ is reseved word
                 "SKCM", "STAD", "THCA", "UCEC", "UCS"]
    db_names = ["COAD"]
    for db_name  in db_names:
        print '==================================='
        mutations_per_gene (cursor, db_name)
        print

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

