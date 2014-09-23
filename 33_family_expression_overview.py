#!/usr/bin/python

import MySQLdb
from   tcga_utils.mysql   import  *
import commands

#########################################
def  expression_stats(cursor, family_name_root):


    qry = "select symbol, fold_change from gene_expression where symbol like '%s%%'" % family_name_root
    rows = search_db(cursor, qry)
 
    if not rows:
        print 'no returnn for'
        print qry
        return
    count = {}
    avg   = {}
    histogram  = {}
    for row in rows:
        [symbol, fold_change] = row
        if not count.has_key(symbol): count[symbol] = 0
        if not avg.has_key(symbol):     avg[symbol] = 0
        if not histogram.has_key(symbol): 
            histogram[symbol] = {}
            for i in range (-4, 5):
                 histogram[symbol][i] = 0
        count[symbol] += 1
        avg[symbol] += fold_change
        for i in range (-4,4):
            if fold_change < i :
                histogram[symbol][i] += 1
                break
        if fold_change >=4:
             histogram[symbol][4] += 1
 
    for [symbol, ct] in count.iteritems():
        avg[symbol] /= ct

    symbols_sorted = sorted(count.keys(), key= lambda x: x )


    for symbol in symbols_sorted:
        print '\t\t"%s" : {' % symbol
        print '\t\t\t"avg" : "%.2f",' % avg[symbol]
        hist = histogram[symbol]
        print '\t\t\t"hist" : {'
        first = True
        for [k,v] in hist.iteritems():
            if first:
                first = False
            else:
                print ","
            print '\t\t\t\t"%d" : "%d"' % (k,v),
        print
        print "\t\t\t}"
        if symbol == symbols_sorted[-1]:
            print "\t\t}"
        else:
            print "\t\t},"
    
    return


#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    family_name_root = 'RPL'

    db_names  = ["BRCA", "COAD", "KIRC", "KIRP", "LUAD", "LUSC",  "OV", "UCEC"]
    #db_names  =  ["KIRC"]

    print 'var gene_expression = { '
    for db_name in db_names:
        print '\t"%s" : {' % db_name
        switch_to_db (cursor, db_name)
        expression_stats(cursor, family_name_root)
        if db_name == db_names[-1]:
            print "\t}"
        else:
            print "\t},"
        
    print '}'

#########################################
if __name__ == '__main__':
    main()
