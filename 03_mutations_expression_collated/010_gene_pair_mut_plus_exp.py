#!/usr/bin/python

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

import sys, os
import MySQLdb
from   tcga_utils.mysql   import  *
from   tcga_utils.utils   import  *
from scipy import stats

#########################################
def expression_status(cursor, hugo_symbol, sample_barcode_short, scaling):

    [ratio, pval] = [0,"1.0"]
    qry  = "select rpkm, experiment_id from rnaseq_rpkm where symbol='%s' " % hugo_symbol
    qry += "and sample_barcode_short = '%s'" % sample_barcode_short
    rows = search_db(cursor, qry)
    if not rows or 'E' in str(rows[0][0]): # we shouldn't be here
        print "no return for '%s' [?!]" % qry
        exit(1)
    [raw_rpkm, exp_id] = rows[0]
    rpkm = raw_rpkm * scaling[exp_id]

    qry = "select mean, KL_pval, interval_endpoints from rnaseq_distro_description where symbol='%s'" % hugo_symbol
    rows = search_db(cursor, qry)
    if not rows or 'E' in str(rows[0][0]): # we shouldn't be here
        print "no return for '%s' [?!]" % qry
        exit(1)
    [mean, KL_pval, interval_endpoints_str] = rows[0]
    # interval enpoints for 99, 95, 90
    if interval_endpoints_str:
        interval_endpoints = [float(x) for x in interval_endpoints_str.split(';')]
        if  rpkm < interval_endpoints[0] or  rpkm > interval_endpoints[1]:
            pval = "<0.01"
        elif rpkm < interval_endpoints[2] or  rpkm > interval_endpoints[3]:
            pval = "<0.05"
        elif rpkm < interval_endpoints[4] or  rpkm > interval_endpoints[5]:
            pval = "<0.10"
        else:
            pval = ">0.10"
    else:
        pval = 'na'
    #print rpkm, mean, KL_pval, interval_endpoints
    if mean > 0.1:
        ratio = rpkm/mean
    else:
        ratio = rpkm
        pval = "mean=0"

    return [rpkm, ratio,pval]

#########################################
def clasf_short (clasf_long):
    if '_Mut' in clasf_long:
        short = clasf_long.replace('_Mutation','')
    elif 'Frame' in clasf_long:
        short = 'FS_Ins'
    else:
        short = clasf_long
    return short

#########################################
def main():

    if len(sys.argv) < 3:
        print  "usage: %s <gene symbol 1>  <gene symbol 2>" % sys.argv[0]
        exit(1)

    full_name = read_cancer_names ()

    # unbuffered output
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    genes = [x.upper() for x in sys.argv[1:]]

    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA",  "GBM", "HNSC", "KICH" ,"KIRC",
                 "KIRP","LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "REA",
                 "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
    db_names = ["BLCA", "BRCA",  "CESC", "CHOL","COAD" , "ESCA",  "GBM", "HNSC", "KICH" ,"KIRC"]
    #db_names = ["BLCA", "COAD","UCEC"]

    cont_table = {}
    for symbol in genes[1:]:
        cont_table[symbol] = {}
        for ref_status in ['mutated','supressed','overexpressed','normal']:
            cont_table[symbol][ref_status] = {}
            for status in ['mutated','supressed','overexpressed','normal']:
                cont_table[symbol][ref_status][status] = 0
    ############################
    total_over = {}  # refers to the first gene which now I expect to be tp53, that stas we'll collect only for dbs where  our gene of interst is sup of over
    total_sup  = {}
    overexpr    = {}
    suppr       = {}
    overexpr_in = {} # list of tumors where symbol is overxpressed
    suppr_in = {} # list of tumors where symbol is overxpressed

    for symbol in genes[1:]:
        total_over[symbol] = {}
        total_over[symbol]["wt"]  = 0
        total_over[symbol]["mut"] = 0
        total_sup[symbol]  = {}
        total_sup[symbol]["wt"]   = 0
        total_sup[symbol]["mut"]  = 0

        overexpr[symbol] = {}
        overexpr[symbol]["wt"]  = 0
        overexpr[symbol]["mut"] = 0
        suppr[symbol] = {}
        suppr[symbol]["wt"]  = 0
        suppr[symbol]["mut"] = 0
        overexpr_in[symbol] = []
        suppr_in[symbol] = []

    for db_name in db_names:
        print "######################################"
        print db_name, full_name[db_name]
        switch_to_db (cursor, db_name)

        ############################
        total_for_this_db   = {}  # refers to the first gene which now I expect to be tp53
        total_for_this_db["wt"]  = 0
        total_for_this_db["mut"] = 0
        overexpr_for_this_db = {}
        suppr_for_this_db    = {}
        for symbol in genes[1:]:
            overexpr_for_this_db[symbol] = {}
            overexpr_for_this_db[symbol]["wt"]  = 0
            overexpr_for_this_db[symbol]["mut"] = 0
            suppr_for_this_db[symbol] = {}
            suppr_for_this_db[symbol]["wt"]  = 0
            suppr_for_this_db[symbol]["mut"] = 0

        ############################
        table = 'rnaseq_rpkm'
        if not check_table_exists (cursor, db_name, table):
            print table, "not found"
            continue
        print "total number of entries in %s:" % table,
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print  rows[0][0]
        if not rows[0][0]: continue

        table = 'somatic_mutations'
        print "total number of entries in %s:" % table,
        qry = "select count(1) from " + table
        rows = search_db(cursor, qry)
        print  rows[0][0]
        if not rows[0][0]: continue

        ############################
        # for how many patients do we have info about mutations and expression?
        qry  = "select distinct sample_barcode_short from somatic_mutations"
        rows = search_db (cursor, qry)
        patients_w_mutations = set([row[0] for row in rows])

        qry  = "select distinct sample_barcode_short from rnaseq_rpkm"
        rows = search_db (cursor, qry)
        patients_w_expression = set([row[0] for row in rows])

        print "number of samples with mutation   info:", len(patients_w_mutations)
        print "number of samples with expression info:", len(patients_w_expression)
        intersection = patients_w_mutations & patients_w_expression
        print "intersection:", len (intersection)

        if not intersection: continue

        ############################
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


        ############################

        print "%10s %7s  " % ('sample id ', '#muts'),
        print "%10s  %10s  %10s  %7s    "  % ('name', "rpkm", "expr/mean", "pval"),
        print "%10s" % ( 'variant(s)')
        ct = 0
        for sample_barcode_short in intersection:
            ct += 1
            if not ct % 10: print "sample", ct
            # how many mutations in this particular sample?
            qry = "select count(1) from %s where sample_barcode_short = '%s'" % (table, sample_barcode_short)
            rows = search_db(cursor, qry)
            if not rows:
                tot_number_of_mutations_in_sample = 0
            else:
                tot_number_of_mutations_in_sample = rows[0][0]
            outstr = ""
            first = True
            one_off = False
            for symbol in genes:
                if first:
                    outstr +=  "%10s %7s  " % (sample_barcode_short, tot_number_of_mutations_in_sample)
                else:
                    outstr +=  "%10s %7s  " % ( "","")
                #####################################
                [rpkm, expr_ratio, expr_pval] = expression_status(cursor, symbol, sample_barcode_short, scaling)
                if  expr_ratio > 3.0:
                    one_off = True
                outstr +=  "%10s  %10.2f %10.2f  %7s     "  % (symbol, rpkm,  expr_ratio, expr_pval)
                #####################################
                qry = "select  variant_classification, aa_change "
                qry += "from %s " % table
                qry += "where hugo_symbol = '%s' " % symbol
                qry += "and sample_barcode_short= '%s' " % sample_barcode_short
                mutations = search_db (cursor, qry)
                if not mutations:
                    mut_string = "none"
                else:
                    mut_string = "; ".join( [m[0]+","+m[1] for m in mutations])
                outstr +=  mut_string
                outstr += "\n"

                if 'Missense' in mut_string or 'Nonsense' in mut_string or 'Frame' in mut_string:
                    status = 'mutated'
                elif expr_ratio< 0.5:
                    status = 'supressed'
                elif expr_ratio>2.0:
                    status = 'overexpressed'
                else:
                    status = 'normal'

                if first: # here I am taking that the first is TP53 or some such reference protein
                    reference_status = status
                    if reference_status=='mutated':
                        total_for_this_db["mut"] += 1
                    else:
                        total_for_this_db["wt"] += 1
                else:
                    cont_table[symbol][reference_status][status] += 1
                    if (expr_pval=="<0.05" or expr_pval=="<0.01"):
                        if expr_ratio > 1.0:
                            if reference_status=='mutated':
                                overexpr_for_this_db[symbol]["mut"]  += 1
                            else:
                                overexpr_for_this_db[symbol]["wt"]  += 1
                        else:
                            if reference_status=='mutated':
                                suppr_for_this_db[symbol]["mut"] += 1
                            else:
                                suppr_for_this_db[symbol]["wt"]  += 1

                first = False

            if False and one_off:
                print outstr
                print

        for symbol in genes[1:]:
            if overexpr_for_this_db[symbol]["mut"] or overexpr_for_this_db[symbol]["wt"]:
                overexpr_in[symbol].append(db_name)
                for stat in ["mut", "wt"]:
                    total_over[symbol][stat] += total_for_this_db[stat]
                    overexpr[symbol][stat] += overexpr_for_this_db[symbol][stat]

            if suppr_for_this_db[symbol]["mut"] or suppr_for_this_db[symbol]["wt"]:
                suppr_in[symbol].append(db_name)
                for stat in ["mut", "wt"]:
                    total_sup[symbol][stat] += total_for_this_db[stat]
                    suppr[symbol][stat] = suppr_for_this_db[symbol][stat]

    for symbol in genes[1:]:
        # overexpression
        if not overexpr_in[symbol]: continue
        outstr = ""
        outstr += "=======================================================================\n"
        outstr += symbol +  " overexpressed in " + " ".join (overexpr_in[symbol])+ "\n"
        a = overexpr[symbol]["wt"]
        b = total_over[symbol]["wt"] - a
        c = overexpr[symbol]["mut"]
        d = total_over[symbol]["mut"] - c
        [odds,pval1] = stats.fisher_exact([[a, b], [c, d]],"greater")
        if pval1 < 0.05:
            outstr += "%d  %d  %d  %d\n" % (a,b,c,d)
            outstr += "hypothesis: overexpression of "+ symbol +" is equally likely for wt and mut p53\n"
            outstr += "(alternative: wt is more likely to have "+ symbol + " overexpressed)\n"
            outstr +=  'fisher pval:  %5.2e\n' % (pval1)
            outstr += "\n"
        [odds,pval2] = stats.fisher_exact([[c,d], [a, b]],"greater" )
        if pval2 < 0.05:
            outstr += "%d  %d  %d  %d\n" % (c,d,a, b)
            outstr += "hypothesis: overexpression of " + symbol + " is equally likely for wt and mut p53 \n"
            outstr += "(alternative: mut is more likely to have "+ symbol +  "overexpressed)\n"
            outstr +=  'fisher pval:  %5.2e\n' % (pval2)
            outstr += "\n"
            #outstr +=  "fisher.test(matrix(c(%d, %d, %d, %d),nrow=2,ncol=2),alternative=\"greater\")" % (a, b, c, d)
        if pval1 < 0.05 or pval2 < 0.05:
            print outstr
        

        # suppression
        if not suppr_in[symbol]: continue
        outstr = ""
        outstr += "=======================================================================\n"
        outstr += symbol +  " suppressed in " + " ".join (overexpr_in[symbol]) + "\n"
        a = suppr[symbol]["wt"]
        b = total_sup[symbol]["wt"] - a
        c = suppr[symbol]["mut"]
        d = total_sup[symbol]["mut"] - c
        [odds,pval1] = stats.fisher_exact([[a, b], [c, d]],"greater")
        if pval1 < 0.05:
            outstr += "%d  %d  %d  %d\n" % (a,b,c,d)
            outstr += "hypothesis: suppression of " + symbol + " is equally likely for wt and mut p53\n"
            outstr += "(alternative: wt is more likely to have " + symbol +  " suppressed)\n"
            outstr +=  'fisher pval:  %5.2e\n' % (pval1)
            outstr +="\n"
        [odds,pval2] = stats.fisher_exact([[c,d], [a, b]],"greater" )
        if pval2 < 0.05:
            outstr += "%d  %d  %d  %d\n" % (c,d,a, b)
            outstr += "hypothesis: suppression of "+ symbol+ " is equally likely for wt and mut p53\n"
            outstr += "(alternative: mut is more likely to have " +  symbol+ " suppressed)\n"
            outstr +=  'fisher pval:  %5.2e\n' % (pval2)
            outstr += "\n"
            #outstr +=  "fisher.test(matrix(c(%d, %d, %d, %d),nrow=2,ncol=2),alternative=\"greater\")" % (a, b, c, d)
        if pval1 < 0.05 or pval2 < 0.05:
            print outstr

        ############################

    #print "contingecy"
    #for symbol in genes[1:]:
    #    for ref_status in ['mutated','supressed','overexpressed','normal']:
    #        for status in ['mutated','supressed','overexpressed','normal']:
    #            print " p53  %15s     %10s %15s     %4d" % (ref_status, symbol, status, cont_table[symbol][ref_status][status])


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()

