#!/usr/bin/python -u

# to speed up things make index:  (int 01_maf_somatic_mutations_table.sql)
# create index mutation_idx on somatic_mutations (tumor_sample_barcode, chromosome, strand, start_position)
# note: I am assuming then that no tumot will have the exact same mutation in both alleles, for the simple reason that I do nto see
# how would the thwo entried in the database then be distinguished from one another
# (rather if tumor_seq_allele1 == tumor_seq_allele2 != match_norm_seq_allele1  and tumor_seq_allele2 != match_norm_seq_allele2
# then I have somehting like that)


import MySQLdb
from   tcga_utils.mysql   import  *
from math  import sqrt
import os
from  scipy import stats
import  numpy

#########################################
def check_barcode(cursor, barcode):
    switch_to_db(cursor, 'tcga_meta')
    qry = "select * from annotation where item_barcode='%s'" % barcode
    rows = search_db(cursor, qry)
    if rows:
        print "annotation found:"
        print rows
        exit(1)

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
def process_data_set (cursor, data_set, parent_dir, data_files):
    print "processing ", data_set

    values = {}
    symbol_ct = 0
    for data_file in data_files:
        fields = data_file.split( '.')
        barcode = fields[1]
        #print barcode
        check_barcode(cursor, barcode)

        full_path = "/".join([parent_dir, data_file])
        inff = open(full_path)
        for line in inff.readlines( )[1:]:
            if line[0] == '?': continue
            line = line.rstrip()
            fields = line.split ('\t')
            if len(fields) != 4: continue # I don't know what this is in that case
            fields_clean = [x.replace("'", '') for x in fields]
            symbol = fields_clean[0].split('|')[0].upper()
            if not values.has_key(symbol):
                values[symbol] = []
            rpkm = float (fields_clean [-1])
            if rpkm > 2000: continue
            values[symbol].append(rpkm)

    tot = 0
    normal = 0
    weak = 0
    weak_with_tail = 0
    not_skewed = 0
    renormalizable = 0
    wide  = 0
    other = 0
    bins = [0,1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            40, 50, 60, 70, 80, 90 ,100, 200, 300, 400, 500, 1000, 10000]

    for symbol, val_array in values.iteritems():
        if bad (cursor, symbol): continue
        print "working on", symbol
        tot += 1
        description = stats.describe(val_array)
        [min,max]  = description[1]

        if max < 3:
            weak += 1
            print "%10s " % symbol, " weak - presumably not expressed; ", "max = %.2f" % max
            continue


        mean = description[2]
        skewness = description[4]
        std = sqrt(description[3])
        if len(val_array) > 20:
            [teststat, pval] = stats.normaltest(val_array)
        else:
            teststat  = 100

        if teststat < 10:
            normal += 1
            print "%10s " %  symbol, "  well defined normal distro - presumably not modified from healthy; ", "mean = %5.2f,  stdev = %5.2f" % (mean, std)
            continue
        # a coarse way of getting rid of outliers:
        if True:
            done = False
            while not done:
                [hist_numpy, low_range, binsize, extrapoints] = stats.histogram (val_array)
                #hist_numpy  is nd_array some numpy piece of shit
                # I wish people would learn to hide their lame garbage from users
                histogram = hist_numpy.tolist()
                # shave off the last value esp if isolated and small
                #print "   ", symbol,  len(val_array), histogram, histogram[-2] , histogram[-1]
                # histogram values are floats (don't ask me) the second condition means 'nothing in the bin [-2]'
                if len(val_array) and  histogram[-2] < 1 and  histogram[-1] <3:
                    # the default number of bins in this story is 10
                    # see documentation for stats.histogram
                    cutoff_value = low_range + binsize*8
                    # lets make new val_array according to that rule
                    val_array = [x for x in val_array if x < cutoff_value]
                    [hist_numpy, low_range, binsize, extrapoints] = stats.histogram (val_array)

                else:
                    done = True
        val_array = [x for x in val_array if mean-5*std < x < mean+5*std]
        # repeat the stats on the reduced array
        description = stats.describe(val_array)
        [min,max]  = description[1]
        if max < 3:
            weak_with_tail += 1
            print  "%10s " % symbol, " after cutting tails off: weak - presumably not expressed; ", "max = %.2f" % max
            continue

        mean = description[2]
        skewness = description[4]
        std = sqrt(description[3])
        if len(val_array) > 20:
            [teststat, pval] = stats.normaltest(val_array)
        else:
            teststat = 100

        if teststat < 10:
            renormalizable += 1
            print  "%10s " % symbol, " after cutting tails off: well defined normal distro - presumably not modified from healthy; ", "mean = %.2f" % mean
            continue

        if abs(skewness) <3:
            width = max - min
            if width < 20:
                print "%10s " % symbol, " after cutting tails off: weaky skewed - narrow;  ",
                print "mean = %5.2f     min = %5.2f   max = %5.2f" %  (mean, min, max)
                not_skewed += 1
            else:
                print "%10s " % symbol, " after cutting tails off: weaky skewed - wide;    ",
                print "mean = %5.2f     min = %5.2f   max = %5.2f" %  (mean, min, max)
                wide += 1

            continue

        hist_numpy = stats.histogram2 (val_array, bins)
        histogram = hist_numpy.tolist()
        print "%10s " % symbol, " undecided "
        if False:
            i = 0
            for val in histogram[:-1]:
                print " %5d %5d " % (bins[i], bins [i+1]),
                print " %5d " %  val
                i += 1
            print
        other += 1


        if  False and mean > 10 and abs(skewness)>5:
            print symbol
            print "\t min max:  %8.2f, %8.2f" % description[1]
            print "\t mean:     %8.2f" % mean
            print "\t stdev:    %8.2f" % std
            print "\t skewness: %8.2f" % description[4]
            print "\t kurtosis: %8.2f" % description[5]
            print 'normaltest teststat = %6.3f pvalue = %6.4f'  % (teststat, pval)
            print

    print "summary for parent directory", parent_dir
    print "tot %d, weak %d,  weak with tail %d,   normal %d,   not_skewed %d, renormalizable: %d,  wide: %d,  other: %d" % \
          (tot,    weak,     weak_with_tail,      normal,       not_skewed,   renormalizable,  wide, other)




#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BRCA","COAD","GBM","KIRC","KIRP","LAML","LGG","LUAD","LUSC","OV","REA","UCEC"]

    switch_to_db (cursor, 'tcga_meta')

    for db_name in db_names:

        print " ** ", db_name

        table = 'gene_expression'

        db_dir  = '/Users/ivana/databases/TCGA/%s/Expression_Genes' % db_name

        for parent_dir, subdir_list, file_list in os.walk(db_dir):
            data_files = filter (lambda x: x[-3:] == 'txt', file_list)
            if data_files:
                data_set = parent_dir
                process_data_set (cursor, data_set, parent_dir, data_files)

        exit(1)

    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
