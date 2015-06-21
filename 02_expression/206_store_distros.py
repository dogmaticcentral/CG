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
from math  import sqrt, floor, ceil, exp
import os
from scipy import stats
from time import time
import matplotlib.pyplot as plt
import numpy as np
from random import sample



output            = True
to_file           = False
use_normal_tissue = False
store_in_db       = False

if store_in_db: use_normal_tissue = False # at least until we create a separate table
# tumor: ['01', '03', '08', '09']
# normal: ['10', '11', '12', '14']


#########################################
def bad (cursor, symbol):
    if symbol[:3] == 'LOC':  return True
    if symbol[0] == 'C' and 'ORF' in symbol: return True
    if symbol[:4]== 'KIAA': return True
    switch_to_db(cursor, 'baseline')
    qry = "select locus_type from hgnc_id_translation where approved_symbol = '%s'" % symbol
    rows = search_db(cursor, qry)
    if not rows: return  True
    if not rows[0][0] == 'gene with protein product': return True

    return False

#########################################
def  blurb (symbol, description, comment, outf, cursor):
    [nobs, [min,max], mean, variance, skewness, kurtosis, normaltest, left, right]  = description
    if output:
        if False:
            print >> outf,"===================================="
            print >> outf,"%s " % symbol
            print >> outf,"pts: %4d    min = %6.2f   max = %10.2f    "   %  ( nobs, min, max)
            print >> outf,"mean = %10.2f   stdev = %6.2f   skew = %5.2f  kurt = %8.2f"   %  (mean, sqrt(variance), skewness, kurtosis)
            print >> outf,"fract left = %5.3f     fract right  = %5.3f " %  (left, right)
            print >> outf, comment
        else:
            print >> outf,"%15s " % symbol,
            print >> outf,"%4d  %6.2f  %10.2f   "   %  ( nobs, min, max),
            print >> outf,"%10.2f   %6.2f   %5.2f   %8.2f"   %  (mean, sqrt(variance), skewness, kurtosis),
            print >> outf,"   %5.3f  %5.3f " %  (left, right)
            #print >> outf, comment

    # if we are storing, do that here too
    if store_in_db:
        fixed_fields  = {'symbol':symbol}
        update_fields = {'number_of_points':nobs, 'min':min, 'max':max, 'mean':mean,
                         'stdev':sqrt(variance), 'skewness':skewness, 'kurtosis':kurtosis, 'normaltest':normaltest}

        ok = store_or_update (cursor, 'rnaseq_distro_description', fixed_fields, update_fields)
        if not ok:
            print 'store failure:'
            print fixed_fields
            print update_fields
            exit(1)



#########################################
def is_bimodal(histogram):
    for i in range (1, len(histogram)):

        val = histogram[i]
        left_max = max(10,val)
        left_max_found = False
        for j in range(0,i):
            if histogram[j] > left_max:
                left_max = histogram[j]
                left_max_found = True
                break

        right_max = max(10,val)
        right_max_found = False
        for j in range(i+1, len(histogram)):
            if histogram[j] > left_max:
                right_max = histogram[j]
                right_max_found = True
                break

        if left_max_found  and  right_max_found:
            return True

    return False

#################################################
def  optimize_distro_domain (symbol,val_array, init_width, mean, stdev, distro, verbose=True):

    if distro=='gamma':
        model_distro = stats.gamma
    elif distro== 'lognorm':
        model_distro = stats.lognorm

    done = False
    round_ct = 0
    right_side = False
    pval_max = -1.0
    prev_pval = -1.0
    pval = -1.0
    opt_bound = [0,len(val_array)]
    opt_params = [0, 0, 0]
    left_cut  = 0
    right_cut = len(val_array)
    # take the core of the data and then work outwards from there
    while (val_array[left_cut]    < mean - stdev*init_width) and left_cut < len(val_array): left_cut += 1
    while (val_array[right_cut-1] > mean + stdev*init_width) and right_cut >1: right_cut -= 1

    left_bound  = [0,left_cut]
    right_bound = [right_cut,len(val_array)]
    init_left_cut  = left_cut
    init_right_cut = right_cut

    while not done or round_ct==len(val_array):
        round_ct += 1
        prev_pval = pval
        [shape, location, scale] =  model_distro.fit(val_array[left_cut:right_cut])
        [D, pval] = stats.kstest(val_array[left_cut:right_cut], distro, (shape, location, scale))

        test1 = (pval < prev_pval)
        #test2 = (location < 0.0)
        test2 = False
        if distro=='gamma':
            test3 = (shape < 1.0)
        else:
            test3  = False

        if pval >= 0.9 and not test2  and  not test3  and right_side==True and round_ct>1:
            if verbose: print "passed the basic citerion - out of the optimization loop with cut values", left_cut, right_cut
            done = True
            break # let's break here not to go crazy with the indentation


        if right_side:
            if round_ct>1  and  (test1 or test3 ): #getting worse, backpedal
                right_bound[1] = right_cut
                right_cut = int (floor ( (right_bound[1] + right_bound[0])/2.0 ))
                if verbose: print "backpedal: right_cut =", right_cut
                if right_bound[1] <= right_bound[0]+1:
                    done=True
                    if verbose: print "the interval closed, were done on the right", right_bound
                if right_cut == init_right_cut:
                    if verbose: print "we're back ot the iniitival value of the right cut (done)", right_cut
                    if verbose: done=True
            else: #this is an improvement, keep moving in the same direction
                if (pval > pval_max and not test2 and not test3):
                    pval_max = pval
                    opt_bound = [left_cut, right_cut]
                    opt_params = [shape, location, scale]
                # but first check if there is room to move
                if right_bound[1] <= right_bound[0]+1:
                    done = True
                    if verbose: print "we're done on the right - bounds:", right_bound
                else:
                    right_bound[0] = right_cut
                    right_cut = int (ceil ( (right_bound[1] + right_bound[0])/2.0 ))
                    if verbose: print "move forward: right_cut =", right_cut

        else: # left_side
            if round_ct>1  and  (test1 or test3 ): #getting worse, backpedal
                left_bound[0] = left_cut
                left_cut = int (ceil ((left_bound[1] + left_bound[0])/2.0))
                if verbose: print "backpedal: left_cut =", left_cut
                if left_bound[1] <= left_bound[0]+1:
                    right_side = True
                    if verbose: print "the interval closed, were done on the left", left_bound
                if  left_cut == init_left_cut:
                    if verbose: print "we're back to the initial value of the left cut (done)", left_cut
                    right_side =True
            else: #this is an improvement, keep moving in the same direction
                if (pval > pval_max and not test2 and not test3):
                    pval_max = pval
                    opt_bound = [left_cut, right_cut]
                    opt_params = [shape, location, scale]
                if left_bound[1] <= left_bound[0]+1:
                    right_side = True
                    round_ct = 0
                    if verbose: print "we're done on the left - bounds:", left_bound
                else:
                    left_bound[1] = left_cut
                    left_cut = int (floor ( (left_bound[1] + left_bound[0])/2.0 ))
                    if verbose: print "move forward: left_cut =", left_cut

    return [pval_max] + opt_bound + opt_params

#########################################
def  optimization_loop (symbol, val_array,  mean, stdev, use_gamma, verbose=False):

    opt_results = []
    pval_max = -1
    for init_width in [10, 5, 2, 1, 0.5]:

        retvals = optimize_distro_domain (symbol, val_array, init_width, mean, stdev, use_gamma, verbose)
        pval = retvals[0]
        if (pval > pval_max):
            pval_max = pval
            opt_results = retvals
        if pval > 0.9:
            break

    return opt_results

#########################################
def store (cursor, symbol, distro, description, distro_params):


    fixed_fields = {'symbol':symbol}
    update_fields = {}
    [nobs, [min,max], mean, variance, skew, kurtosis]  = description
    update_fields['number_of_points'] =  nobs
    update_fields['min']      =  min
    update_fields['max']      =  max
    update_fields['mean']     =  mean
    update_fields['stdev']    =  sqrt(variance)
    update_fields['skewness'] =  skew
    update_fields['kurtosis'] =  kurtosis
    update_fields['distro']   =  distro

    if distro_params:
        update_fields['KL_pval']   = distro_params[0]
        update_fields['left_cut']  = distro_params[1]
        update_fields['right_cut'] = distro_params[2]
        update_fields['shape']     = distro_params[3]
        update_fields['location']  = distro_params[4]
        update_fields['scale']     = distro_params[5]
        update_fields['interval_endpoints'] = distro_params[6]

    ok = store_or_update(cursor, 'rnaseq_distro_description', fixed_fields, update_fields)
    if not ok:
        print 'store failure:'
        print fixed_fields
        print update_fields
        exit(1)

#########################################
def   store_fitted (cursor, symbol, distro, description,  opt_result_lognorm):
    if distro=="lognorm":
        interval = stats.lognorm.interval
    elif distro=="gamma":
        interval = stats.gamma.interval
    else:
        print "unrecognized distro:", distro
        exit(1)
    [left_cut, right_cut] = opt_result_lognorm[1:3]
    [shape, location, scale] = opt_result_lognorm[3:]
    intervals_blob = ""
    for alpha in [99, 95, 90]:
        if intervals_blob: intervals_blob += ";"
        intervals_blob += "%.2f;%.2f" % tuple(interval(alpha/100.0, shape, location, scale))
    distro_params = opt_result_lognorm + [intervals_blob]
    store (cursor, symbol, distro, description, distro_params)


#########################################
def process_data_set (cursor, db_name, gene_list):

    switch_to_db(cursor, db_name)

    if use_normal_tissue:
        qry = "select count(1) from rnaseq_rpkm where source_code in (10,11,12,14)"
        rows = search_db(cursor, qry)
        if not rows or not rows[0][0] > 0:
            print "no normal tissue samples in ", db_name
            return
    values   = {}
    no_match = []
    start = time()
    ct = 0

    # get  the scaling factors
    qry = "select * from rnaseq_scaling"
    rows = search_db(cursor, qry)

    if not rows:
        print "no scaling info"
        exit(1)

    scaling = {}

    for row in rows:
        [sample_id, scaling_factor] = row
        scaling[sample_id] = scaling_factor

    #toyset = sample(gene_list, 100)
    toyset = ['TP53', 'RPL5', 'RPL11', "RPL22", 'WNT11', 'WLS', "PORCN", 'MDM2',
              'CDKN2A', 'ACTN', 'ACTB', 'LEP', 'GPR161', 'CDK2', 'HIF1A',
              'ERBB2', 'PTEN','CDKN2A', 'LTN1','NEMF', 'NMD3','EIF6'] # melanoma

    for symbol in toyset:
        ct += 1
        if not ct%1000: print "%d done in %5.2fs" % (ct, time() - start)
        if use_normal_tissue:
            qry = "select rpkm, experiment_id from rnaseq_rpkm  where symbol='%s' and source_code in (10,11,12,14)" % symbol
        else:
            qry = "select rpkm, experiment_id from rnaseq_rpkm  where symbol='%s' and source_code in (1,3,8,9)" % symbol
        rows = search_db(cursor, qry)
        if not rows:
            #print "no match for ", symbol
            no_match.append(symbol)
            continue
        values[symbol] = []
        for row in rows:
            [rpkm, sample_id] = row
            values[symbol].append(rpkm*scaling[sample_id])

        #print "search done in %5.2fs" % (time() - start)
        # TODO looks like some name esolution is in order here in general

        tot    = 0

        if not values.has_key(symbol): continue
        val_array = values[symbol]
        if bad (cursor, symbol): continue
        if not val_array: continue
        switch_to_db(cursor, db_name)
        tot += 1
        description = stats.describe(val_array)
        [nobs, [min,max], mean, variance, skewness, kurtosis]  = description
        stdev = sqrt(variance)

        if mean < 2:
            #print symbol, "weak"
            store (cursor, symbol, "weak", description, [])
            continue

        if nobs < 20:
            #print "too few samples for a fit - bailing out"
            store (cursor, symbol, "small", description, [])
            continue

        val_array = sorted(val_array)

        distro = 'lognorm'
        opt_result_lognorm = optimization_loop (symbol, val_array, mean, stdev, distro, verbose=False)
        lognorm_ok = (len(opt_result_lognorm)==6)

        # if we are happy with the lognorm fit, we won't do gamma fitting at all
        if lognorm_ok and opt_result_lognorm[0] > 0.9:
            #print symbol, "lognorm p: %5.2f" % opt_result_lognorm[0]
            store_fitted (cursor, symbol, "lognorm", description,  opt_result_lognorm)
            continue

        distro = 'gamma'
        opt_result_gamma = optimization_loop (symbol, val_array, mean, stdev, distro, verbose=False)
        gamma_ok = (len(opt_result_gamma)==6)

        # use gamma if it worked and lognorm failed
        if gamma_ok and not lognorm_ok:
            store_fitted (cursor, symbol, "gamma", description,  opt_result_gamma)
            #print "use gamma bcs lognorm failed"
            continue

        # one more chance for gamma if it is better fit, and covers larger interval than lognorm
        use_gamma = False
        if gamma_ok and lognorm_ok:
            [left_cut_g, right_cut_g] = opt_result_gamma[1:3]
            [left_cut_l, right_cut_l] = opt_result_lognorm[1:3]
            use_gamma = (opt_result_gamma[0] - opt_result_lognorm[0]) > 0.1
            use_gamma = use_gamma and left_cut_g < left_cut_l  and right_cut_l < right_cut_g

        if use_gamma:
            store_fitted (cursor, symbol, "gamma", description,  opt_result_gamma)
            #print "use gamma bcs it is a batter fit"
            continue
        elif lognorm_ok: # lognorm did return something though we were hoping for something better
            store_fitted (cursor, symbol, "lognorm", description,  opt_result_lognorm)
            #print "use lognorm after all %5.2f" % opt_result_lognorm[0]
            continue

        else:
            #print "no lognorm, no gamma"
            store (cursor, symbol, "fail", description, [])




#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BLCA", "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "REA", "UCEC"]
    db_names  = ["BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "REA", "UCEC"]
    # gene names
    gene_list = []
    if False:
        switch_to_db(cursor, 'baseline')
        qry = "select distinct approved_symbol from hgnc_id_translation where locus_type = 'gene with protein product' "
        rows = search_db(cursor, qry)
        gene_list = [row[0] for row in  rows]
        print "\t number of genes", len(gene_list)


    for db_name in db_names:
        print " ** ", db_name

        process_data_set (cursor, db_name, gene_list)


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
