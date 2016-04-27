#!/usr/bin/python -u

#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

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
    elif distro == 'burr':
        model_distro = stats.burr
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
        if distro == 'burr':
            [c, d, location, scale] =  model_distro.fit(val_array[left_cut:right_cut])
            [D, pval] = stats.kstest(val_array[left_cut:right_cut], distro, (c, d, location, scale))
        else:
            [shape, location, scale] =  model_distro.fit(val_array[left_cut:right_cut])
            [D, pval] = stats.kstest(val_array[left_cut:right_cut], distro, (shape, location, scale))

        if False:
            [fit_mean, fit_variance,fit_skew, fit_kurtosis] = model_distro.stats(shape, location, scale, 'mvsk')
            if verbose:
                print
                print symbol, "   round %4d   left cut %3d   right cut %3d  array length: %4d  " % (round_ct,  left_cut, right_cut, len(val_array))
                print "left bounds", left_bound, " right bounds", right_bound
                print  "fitted  shape %6.2f  location %6.2f   scale %6.2f    KS pval: %6.2f " % (shape, location, scale, pval)
                print  "fitted  mean %6.2f      stdev %6.2f   skew %6.2f        kurt: %6.2f " % (fit_mean, sqrt(fit_variance), fit_skew, fit_kurtosis)

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
                    if distro=='burr':
                        opt_params = [c, d, location, scale]
                    else:
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
                    if distro=='burr':
                        opt_params = [c, d, location, scale]
                    else:
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

        retvals = optimize_distro_domain (symbol, val_array, init_width, mean, stdev, use_gamma, verbose=False)
        pval = retvals[0]
        if (pval > pval_max):
            pval_max = pval
            opt_results = retvals
        if pval > 0.9:
            break

    return opt_results


#########################################
def process_data_set (cursor, db_name, gene_list, outf):

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
        #scaling[sample_id] = 1.0
        scaling[sample_id] = scaling_factor

    toyset = [ "RPL22", 'WNT11', 'WLS', 'MDM2', 'CDKN2A', 'ACTN', 'ACTB', 'LEP']
    toyset = [ 'GPR161'] # overexpressed in breast cancer
    toyset = [ 'CDK2'] # overexpressed in COAD (should be, but is not
    toyset = [ 'HIF1A'] # cancers in general
    toyset = ['ERBB2'] # brca
    toyset = ['PTEN'] # brca -loss
    toyset = ['CDKN2A'] # melanoma
    #toyset = sample(gene_list, 100)
    toyset = ["AKAP13", "ESR1", "HAND2", "PRKACA", "PRKAR2A", "PRKAR2B", "PRKCA"]


    for symbol in toyset:
        ct += 1
        if not ct%1000: print "%d done in %5.2fs" % (ct, time() - start)
        if use_normal_tissue:
            qry = "select rpkm, experiment_id from rnaseq_rpkm  where symbol='%s' and source_code in (10,11,12,14)" % symbol
        else:
            qry = "select rpkm, experiment_id from rnaseq_rpkm  where symbol='%s' and source_code in (1,3,8,9)" % symbol
        rows = search_db(cursor, qry)
        if not rows:
            print "no match for ", symbol
            no_match.append(symbol)
            continue
        values[symbol] = []
        for row in rows:
            [rpkm, sample_id] = row
            values[symbol].append(rpkm*scaling[sample_id])

        #print "search done in %5.2fs" % (time() - start)
        # TODO looks like some name esolution is in order here in general

        tot    = 0
        normal = 0
        weak   = 0
        weak_with_tail = 0
        not_skewed     = 0
        renormalizable = 0
        wide  = 0
        other = 0

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
            print "weak"
            weak += 1
            description += (0.0, 0.0, 0.0)
            blurb (symbol, description, "weak - presumably not expressed", outf, cursor)
            continue

        if nobs < 20:
            print "too few samples for a fit - bailing out"
            continue

        val_array_orig = val_array[:] # this makes a copy of the original list

        description += (0.0, 0.0, 0.0)
        blurb (symbol, description, "", outf, cursor)

        val_array = sorted(val_array)

        print symbol
        print ("%12s" * 7)% ("distro", "pval", "left_cut", "right_cut", "shape", "location", "scale")

        distro = 'gamma'
        opt_result_gamma = optimization_loop (symbol, val_array, mean, stdev, distro, verbose=False)
        print "%12s" % "gamma",
        if len(opt_result_gamma) <6:
            gamma_ok = False
            print " failed"
        else:
            gamma_ok = True
            print "%12.2f" % opt_result_gamma[0],
            print "%12d%12d" % tuple(opt_result_gamma[1:3]),
            print ("%12.2f" * 3) % tuple (opt_result_gamma[3:])

        distro = 'lognorm'
        opt_result_lognorm = optimization_loop (symbol, val_array, mean, stdev, distro, verbose=False)
        print "%12s" % "lognorm",
        if len(opt_result_lognorm) <6:
            lognorm_ok = False
            print " failed"
        else:
            lognorm_ok = True
            print "%12.2f" % opt_result_lognorm[0],
            print "%12d%12d" % tuple(opt_result_lognorm[1:3]),
            print ("%12.2f" * 3) % tuple (opt_result_lognorm[3:])
        print

        orig_max = val_array[-1]
        if lognorm_ok:
            [left_cut, right_cut] = opt_result_lognorm[1:3]
            [shape, location, scale] = opt_result_lognorm[3:]
            print "lognorm intervals:"
            for alpha in [99, 95, 90]:
                print alpha, stats.lognorm.interval(alpha/100.0, shape, location, scale)
        elif gamma_ok:
            [left_cut, right_cut] = opt_result_gamma[1:3]
            [shape, location, scale] = opt_result_gamma[3:]
            print "gamma intervals:"
            for alpha in [99, 95, 90]:
                print alpha, stats.gamma.interval(alpha/100.0, shape, location, scale)
        else:
            [left_cut, right_cut] = [0, len(val_array)]
        val_left  = val_array[0:left_cut]
        val_right = val_array[right_cut:]
        val_array = val_array[left_cut:right_cut] # note array cut from this point on
        if not val_array: val_array = val_array_orig

        # for concise description of gamma fn using the nomenclature comparapble to scipy.stats, see
        # http://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm
        # the fitting function returns alpha, locus,and beta values
        fig, ax1 = plt.subplots(1, 1)
        title = db_name + ", " + symbol
        if use_normal_tissue:
            title += ", \"normal\""
        else:
            title += ", tumor"
        ax1.set_title (title, fontsize = 24)
        ax1.set_xlabel('RPKM', fontsize = 22)
        ax1.set_ylabel('Number of cases (binned)', fontsize = 22)
        x = np.linspace(0, orig_max+1, 100)

        orig_max += 1 # the numerics can and will get in trouble otherwise
        binsize = (orig_max - 0 )/20
        bins = [i*binsize for i in range(21)]
        histogram  = np.histogram(val_array, bins)[0]
        if histogram[0] > 0.5* len(val_array):
            binsize = (orig_max - 0 )/20.0
            bins = [i*binsize for i in range(21)]
            histogram  = np.histogram(val_array, bins)[0]
        # finding the max of the  gamma distr is possible if shape is > 1
        # otherwise th thing diverges around zero (or even some larger number,
        # see https://en.wikipedia.org/?title=Gamma_distribution#/media/File:Gamma_distribution_pdf.svg)
        # note: wikipedia uses location = 0 without  mentioning the fact

        #############################################
        # gamma:
        if opt_result_gamma:
            model_distro = stats.gamma
            [shape, location, scale] =  opt_result_gamma[3:]
            if shape > 1:
                mode = (shape - 1)*scale + location
                peak_pdf_value = model_distro.pdf([mode], shape,location, scale)[0]

            else:
                # let's use the value in the middle of the first bin
                peak_pdf_value = model_distro.pdf([binsize/2.0], shape,location, scale)[0]
            peakval = np.amax(histogram)/peak_pdf_value
            y = [p*peakval for p in model_distro.pdf(x, shape,location, scale)]
            ax1.plot (x, y, 'r-', lw=5, alpha=0.6, label='gamma')

        #############################################
        # lognormal
        # mode = exp(mu - sigma^2) --- mu is log scale
        if opt_result_lognorm:
            model_distro = stats.lognorm
            [shape, location, scale] = opt_result_lognorm[3:]
            mode = scale*exp(-shape*shape) + location
            peak_pdf_value = model_distro.pdf([mode], shape,location, scale)[0]
            peakval = np.amax(histogram)/peak_pdf_value
            y = [p*peakval for p in model_distro.pdf(x, shape,location, scale)]
            ax1.plot (x, y, 'r-', lw=5, alpha=0.6, label='lognorm', color='b')


        ax1.legend()

        #############################################
        # histograms
        ax1.hist(val_array, bins,  normed=False, histtype='stepfilled', alpha=0.2)
        if val_left:  ax1.hist(val_left, bins,  normed=False, histtype='stepfilled', alpha=0.2, color='g')
        if val_right: ax1.hist(val_right, bins,  normed=False, histtype='stepfilled', alpha=0.2, color='r')

        plt.ylim(0, np.amax(histogram)*1.1)

        filename = "%s_%s.png" % (db_name, symbol)
        if filename:
            plt.savefig(filename)
        else:
            plt.show()





#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BLCA", "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "REA", "UCEC"]
    #db_names  = ["BRCA"]
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
        if to_file:
            if use_normal_tissue:
               outf = open("output/%s.rnaseq.normal.table"%db_name, "w")
            else:
                outf = open("output/%s.rnaseq.table"%db_name, "w")
        else:
            outf = sys.stdout

        if output:
            print >> outf,  "%15s %4s  %6s  %10s     %10s   %6s  %5s     %8s     %5s  %5s "  \
                            % ('symbol', 'pts', 'min', 'max', 'mean', 'stdev', 'skew', 'kurt', 'left', 'right')

        process_data_set (cursor, db_name, gene_list, outf)

        if to_file: outf.close()


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
