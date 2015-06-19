#!/usr/bin/python -u

#/usr/local/bin/python -u
# this version of python has scipy 0.15.1 - earlier versions had some bug in lognorm; but it dosn't have mysqlDB

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)

# the lognormal is more right skewed (CV3+3CV) than the gamma (2CV).

import MySQLdb
from   tcga_utils.mysql   import  *
from math  import sqrt, floor, ceil, exp, log
import os
from scipy import stats
from time import time
from numpy import linalg, histogram
from sympy import Matrix

output            = True
to_file           = False

use_normal_tissue = True
use_metastatic = False

store_in_db       = False

if store_in_db: use_normal_tissue = False # at least until we create a separate table
# tumor: ['01', '03', '08', '09']
# normal: ['10', '11', '12', '14']

#########################################
def store (cursor, sample_ids, scaling_values):

    for sample_id in sample_ids:

        scaling   = scaling_values[sample_id]
        fixed_fields  = {'experiment_id' : sample_id}
        update_fields = {'scaling_factor': scaling}

        ok = store_or_update (cursor, 'rnaseq_scaling', fixed_fields, update_fields)
        if not ok:
            print 'store failure:'
            print fixed_fields
            print update_fields
            exit(1)

#########################################
def process_data_set (cursor, db_name, gene_list):

    switch_to_db(cursor, db_name)

    no_match = []
    start = time()
    ct = 0

    toyset =  gene_list

    # how many different samples do we have here?
    if use_normal_tissue:
        qry = "select distinct experiment_id from rnaseq_rpkm where source_code in (10,11,12,14)"
    elif use_metastatic:
        qry = "select distinct experiment_id from rnaseq_rpkm where source_code in (6, 7)"
    else:
        qry = "select distinct experiment_id from rnaseq_rpkm where source_code in (1,3,8,9)"
    rows = search_db(cursor, qry)

    if not rows:
        print "no sample info (?!)"
        return

    no_of_samples = len(rows)
    sample_ids = [row[0] for row in rows]
    print "no_of_samples = ", no_of_samples

    rpkms = {}
    for symbol in toyset:

        ct += 1
        if not ct%1000: print "%d done in %5.2fs" % (ct, time() - start)
        if use_normal_tissue:
            qry = "select rpkm, experiment_id from rnaseq_rpkm  where symbol='%s' and source_code in (10,11,12,14)" % symbol
        else:
            qry = "select rpkm,  experiment_id  from rnaseq_rpkm  where symbol='%s' and source_code in (1,3,8,9)" % symbol
        rows = search_db(cursor, qry)
        if not rows:
            no_match.append(symbol)
            continue

        rpkms[symbol] = {}

        for row in rows:
            [rpkm, sample] = row
            rpkms[symbol][sample] = log(rpkm+0.01)
            #rpkms[symbol][sample] = rpkm

    print "not_found: ", len (no_match)
    # check that all sets are the same size - i.e. that for entry for which
    # there is a return  from each sample
    incomplete = 0
    for symbol, values_per_sample in rpkms.iteritems():
        if len(values_per_sample) != no_of_samples:
            incomplete += 1
            #print "number of points mismatch:", symbol, len (values_per_sample), no_of_samples
    print "incomplete:", incomplete

    avg   = {}
    stdev = {}
    for symbol, values_per_sample in rpkms.iteritems():
        val_array = values_per_sample.values()
        [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(val_array)
        avg[symbol] = mean
        stdev[symbol] = sqrt (variance)

    [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(avg.values())
    avg_avg   = mean
    avg_stdev = sqrt(variance)

    [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(stdev.values())
    stdev_avg   = mean
    stdev_stdev = sqrt(variance)

    well_behaved = [symbol for symbol in rpkms.keys() if
                    2 < avg[symbol] < log(1000)  and   stdev[symbol] <  stdev_avg + stdev_stdev/2]

    print
    print "well behaved:", len(well_behaved)

    # new avg
    avg   = {}
    stdev = {}
    grand_avg = 0
    for symbol in well_behaved:
        val_array = rpkms[symbol].values()
        [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(val_array)
        avg[symbol] = mean
        stdev[symbol] = sqrt (variance)
        grand_avg += avg[symbol]
    grand_avg /= len(well_behaved)

    for sample_id in sample_ids:
        val_array = [rpkms[symbol][sample_id] for symbol in well_behaved]
        [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(val_array)
        avg[sample_id] = mean

    if True: # this is very much correct for big number of samples, and is useful as a correction
        correction = {}
        for sample_id in sample_ids:
            correction[sample_id] = grand_avg - avg[sample_id]
        val_array = correction.values()
        [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(val_array)
        print "correction range:   %5.2f  %5.2f  (  %5.2f  %5.2f ) " % (min, max, exp(min), exp(max))


    print "dropping samples with scaling factor > 5 and rescaling - keep in mind when using the data"
    sample_ids_clean = [sample_id for sample_id in sample_ids if exp(correction[sample_id])<= 5]


    # this is robust for small number of samples - least square fit, scaling factor as an additive term
    # for log-transformed case
    S = len(sample_ids_clean)
    M = [[0]*S for i in range(S)]
    b = [0]*S
    for i in range(S):
        M[i][i] = 1.0 - 1./S
        for j in range(i+1,S):
            M[i][j] = -1./S
            M[j][i] = M[j][i]
        b[i] = grand_avg - avg[sample_ids_clean[i]]

    try:
        solvec = linalg.solve(M,b)
    except:
        print "linalg fail"
        exit(1)
    [nobs, [min,max], mean, variance, skewness, kurtosis] = stats.describe(solvec)
    print "correction range:   %5.2f  %5.2f  (  %5.2f  %5.2f ) " % (min, max, exp(min), exp(max))
    print

    scaling = {}
    for i in range(S):
        sample_id = sample_ids_clean[i]
        scaling[sample_id] = exp(solvec[i])

    for sample_id in sample_ids:
        if sample_id in sample_ids_clean: continue
        scaling[sample_id] = exp(correction[sample_id])

    store (cursor, sample_ids, scaling)

#########################################
def make_scaling_table(cursor, db_name):

    switch_to_db (cursor, db_name)
    qry = "select database()"
    rows = search_db(cursor, qry)
    print qry
    print rows

    qry = "";
    qry += "CREATE TABLE rnaseq_scaling ("
    qry += "	 experiment_id VARCHAR(100) NOT NULL,"
    qry += "	 scaling_factor FLOAT (10,2) DEFAULT 1.0"
    qry += ") ENGINE=MyISAM"
    rows = search_db(cursor, qry)
    print qry
    print rows


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BLCA", "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "REA", "UCEC"]
    #db_names  = ["UCEC" ]
    #db_names  = []
    # gene names
    gene_list = []
    if True:
        switch_to_db(cursor, 'baseline')
        qry = "select distinct approved_symbol from hgnc_id_translation where locus_type = 'gene with protein product' "
        rows = search_db(cursor, qry)
        gene_list = [row[0] for row in  rows]
        print "\t number of genes", len(gene_list)

    for db_name in db_names:
        print " ** ", db_name
        table = 'rnaseq_scaling'
        if ( check_table_exists (cursor, db_name, table)):
            print table, " found in ", db_name
            print "unreasonable scaling values:"
            qry = "select * from rnaseq_scaling where scaling_factor > 5"
            print search_db(cursor, qry)
        else:
            print table, " not found in ", db_name
            make_scaling_table(cursor, db_name)

        process_data_set (cursor, db_name, gene_list)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
