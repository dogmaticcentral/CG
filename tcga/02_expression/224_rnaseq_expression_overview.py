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

import MySQLdb
from   tcga_utils.mysql   import  *
from math import log

# this should perhaps be number_of_bins_on_each_side
# i.e. the total number of bins will be  2*number_of_bins + 1 (for the zeroth bin)
number_of_bins = 10
bin_width = 2.0/number_of_bins

#########################################
def  expression_stats(cursor, gene_names, overall_histogram):

    count = {}
    avg   = {}
    avg_sq     = {}
    histogram  = {}
    for symbol in gene_names:

        qry  = "select rpkm, sample_id from rnaseq_rpkm "
        qry += "where symbol = '%s' " % symbol
        qry += "and source_code = 11"
        rows = search_db(cursor, qry)
        if not rows:
            print "no norms for", symbol
            print qry
            exit(1)
            continue
        normal_value = {}
        for row in rows:
            [rpkm, sample_id] = row
            normal_value[sample_id] = rpkm

        qry  = "select rpkm, paired_sample_id from rnaseq_rpkm "
        qry += "where symbol = '%s' " % symbol
        qry += "and not source_code = 11"
        rows = search_db(cursor, qry)
        if not rows: continue

        fold_change = []
        for row in rows:
            [rpkm, paired_sample_id] = row
            if not normal_value.has_key(paired_sample_id):
                print 'Error: no matching normal'
                exit(1)
            norm = normal_value[paired_sample_id]
            # I am trying to avoid ascribing meaning to zeros:
            if norm < 10 and rpkm < 10:
                fc = 0
            else:
                if norm < 0.01: norm = 0.01
                if rpkm < 0.01: rpkm = 0.01
                fc = log(rpkm/norm,2)
            fold_change.append(fc)


        if not count.has_key(symbol):   count[symbol] = 0
        if not avg.has_key(symbol):     avg[symbol] = 0
        if not avg_sq.has_key(symbol):  avg_sq[symbol] = 0
        if not histogram.has_key(symbol):
            histogram[symbol] = {}
            for i in range (-number_of_bins, number_of_bins+1):
                 histogram[symbol][i] = 0
        for fc in fold_change:
            count[symbol] += 1
            avg[symbol]    += fc
            avg_sq[symbol] += fc*fc
            for i in range (-number_of_bins, number_of_bins):
                if fc < i*bin_width + bin_width/2 :
                    histogram[symbol][i] += 1.0
                    break

            if fc >=  (number_of_bins-1)*bin_width + bin_width/2:
                 histogram[symbol][number_of_bins] += 1.0

        for i in range (-number_of_bins, number_of_bins+1):
            overall_histogram[i] += histogram[symbol][i]


    for [symbol, ct] in count.iteritems():
        avg[symbol] /= ct
        avg_sq[symbol] /= ct
        for i in range (-number_of_bins, number_of_bins+1):
            histogram[symbol][i] /= ct


    symbols_sorted = sorted(count.keys(), key= lambda x: x )


    for symbol in symbols_sorted:
        print '\t\t"%s" : {' % symbol
        print '\t\t\t"no_cases" : "%d",' % ct
        print '\t\t\t"avg" : "%.2f",' % avg[symbol]
        print '\t\t\t"avg2" : "%.2f",' % avg_sq[symbol]
        hist = histogram[symbol]
        print '\t\t\t"hist" : ['
        first = True
        for [k,v] in hist.iteritems():
            if first:
                first = False
            else:
                print ","
            print '\t\t\t\t{x:%.2f, y:%d}' % (k*bin_width, int (v*100)),
        print
        print "\t\t\t]"
        if symbol == symbols_sorted[-1]:
            print "\t\t}"
        else:
            print "\t\t},"

    
    return


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()


    #db_names  = ["BLCA","BRCA","HNSC","KIRC","LIHC","LUAD","LUSC","UCEC"]
    db_names  = ["BLCA","BRCA","HNSC","LUAD","LUSC","UCEC"]

    family_name_root = 'RPL'
    switch_to_db (cursor, 'baseline')
    qry  = " select approved_symbol from hgnc_id_translation "
    qry += "where approved_symbol like '%s%%' " % family_name_root
    qry += "and locus_type = 'gene with protein product' and not approved_symbol like '%withdrawn'"
    rows = search_db(cursor, qry)
    if not rows:
        print "no gene names found in baseline.hgnc_id_translation corresponding to root name", family_name_root
        exit(1)
    gene_names = [row[0] for row in rows]

    #gene_names = gene_names [:3]

    overall_histogram = {}
    for i in range (-number_of_bins, number_of_bins+1):
        overall_histogram[i] = 0

    # output - js formatted for use wth expression.js
    print 'var cancers = [',
    for db_name in db_names[:-1]:
        print '"%s",' %db_name,
    print '"%s"];' % db_names[-1]

    print 'var genes = [',
    for gene_name in gene_names[:-1]:
        short = gene_name.replace('RPL', '')
        print '"%s",' % short,
    print '"%s"];' % gene_names[-1]


    print 'var gene_expression = { '
    for db_name in db_names:
        print '\t"%s" : {' % db_name
        switch_to_db (cursor, db_name)
        expression_stats(cursor, gene_names, overall_histogram)
        print "\t},"

    tot = 0
    for i in range (-number_of_bins, number_of_bins+1):
        tot += overall_histogram[i]
    for i in range (-number_of_bins, number_of_bins+1):
        overall_histogram[i] /=  tot

    print '\t"%s" : [' % 'overall'
    first = True
    for [k,v] in overall_histogram.iteritems():
        if first:
            first = False
        else:
            print ","
        print '\t\t{x:%.2f, y:%d}' % (k*bin_width,int (v*100)),
    print
    print "\t]"

    print '}'

#########################################
if __name__ == '__main__':
    main()
