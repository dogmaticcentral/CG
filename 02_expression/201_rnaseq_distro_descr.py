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

use_normal_tissue = False
store_in_db       = True

if store_in_db: use_normal_tissue = False # at least until we create a separate table

#########################################
def check_barcode(cursor, barcode):

    # is this primary tumor
    elements = barcode.split('-')
    source  = elements[3][:-1]
    # source can signa additional or metastatic tumors from the same patient
    # to keep our life simple we'll just stick to primary tumors
    # indicated by source code 01, 03, 08, or 09
    if use_normal_tissue:
        if  not source in ['10', '11', '12', '14']:
            #print "non-primary:", barcode
            return False
    else:
        if  not source in ['01', '03', '08', '09']:
            print "non-primary:", barcode
            return False

    # are there any other annotations (usually indicating something is off with the sample
    switch_to_db(cursor, 'tcga_meta')
    qry = "select * from annotation where item_barcode='%s'" % barcode
    rows = search_db(cursor, qry)
    if rows:
        print "annotation found:"
        print rows
        return False

    return True


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
    [nobs, [min,max], mean, variance, skewness, kurtosis, left, right]  = description
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
        update_fields = {'min':min, 'max':max, 'mean':mean, 'stdev':sqrt(variance), 'skewness':skewness, 'kurtosis':kurtosis}

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

#########################################
def process_data_set (cursor, db_name, parent_dirs, files_in_data_dir, outf):

    values = {}
    symbol_ct = 0

    for parent_dir in parent_dirs:
        # the magical formula in check_barcode gives me the barcode from the name of the file
        filtered_files = [x for x in files_in_data_dir[parent_dir] if check_barcode(cursor, x.split('.')[1]) ]
        print "parent dir:", parent_dir, "number of files:", len(filtered_files)
        for data_file in filtered_files:
            fields = data_file.split( '.')
            barcode = fields[1]
            #print barcode
            if not check_barcode(cursor, barcode): continue

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
                #if rpkm > 2000:
                #    print "high val: %10.2f" %  rpkm
                values[symbol].append(rpkm)

    tot = 0
    normal = 0
    weak = 0
    weak_with_tail = 0
    not_skewed = 0
    renormalizable = 0
    wide  = 0
    other = 0

    for symbol, val_array in values.iteritems():
        if bad (cursor, symbol): continue
        if not val_array: continue

        switch_to_db(cursor, db_name)
        tot += 1
        description = stats.describe(val_array)
        [nobs, [min,max], mean, variance, skewness, kurtosis]  = description
        stdev = sqrt(variance)

        if mean < 3:
            weak += 1
            description += (0.0, 0.0)
            blurb (symbol, description, "weak - presumably not expressed", outf, cursor)
            continue

        in_left_tail  = float (len([x for x in val_array if x < mean-2*stdev ]))/len(val_array)
        in_right_tail = float (len([x for x in val_array if x > mean+2*stdev ]))/len(val_array)
        description += (in_left_tail, in_right_tail)
        if len(val_array) > 20:
            [teststat, pval] = stats.normaltest(val_array)
        else:
            teststat  = 100

        if teststat < 10:
            normal += 1
            blurb (symbol, description, "well defined normal distro - presumably not modified from healthy ", outf, cursor)
            continue

        # a coarse way of getting rid of outliers:
        done = False
        val_array_modified = False
        val_array_orig = val_array[:] # this makes a copy of the original list
        while not done:
            [hist_numpy, low_range, binsize, extrapoints] = stats.histogram (val_array)
            #hist_numpy  is nd_array some numpy piece of shit
            # I wish people would learn to hide their lame garbage from users
            histogram = hist_numpy.tolist()
            # shave off the last value esp if isolated and small
            #print "   ", symbol,  len(val_array), histogram, histogram[-2] , histogram[-1]
            # histogram values are floats (don't ask me) the second condition means 'nothing in the bin [-2]'
            start_len = len(val_array)
            if len(val_array) and  histogram[-2] < 1 and  histogram[-1] <3:
                # the default number of bins in this story is 10
                # see documentation for stats.histogram
                cutoff_value = low_range + binsize*8
                # lets make new val_array according to that rule
                val_array = [x for x in val_array if x < cutoff_value]
                [hist_numpy, low_range, binsize, extrapoints] = stats.histogram (val_array)
                val_array_modified = True
                done =  start_len == len(val_array) # if we ar not changing, we are done
                if ( float(start_len-len(val_array))/start_len > 0.1):
                    # we are chopping off more than 10% of the original array
                    # in this case also revert to the original array
                    done = True
                    val_array_modified = False
                    val_array = val_array_orig
            else:
                done = True


        if val_array_modified:
            # repeat the stats on the reduced array
            description = stats.describe(val_array)
            # nos == number of observations, i.e. number of data points
            [nobs, [min,max], mean, variance, skewness, kurtosis]  = description
            stdev = sqrt(variance)
            if mean < 3:
                weak_with_tail += 1
                description += (0.0, 0.0)
                blurb (symbol, description, "tails off:  weak - presumably not expressed", outf, cursor)
                continue
            # note - how many elements stick out from the original array, but using the new stdev
            in_left_tail  = float (len([x for x in val_array_orig if x < mean-2*stdev ]))/len(val_array_orig)
            in_right_tail = float (len([x for x in val_array_orig if x > mean+2*stdev ]))/len(val_array_orig)

            description += (in_left_tail, in_right_tail)

            if len(val_array) > 20:
                [teststat, pval] = stats.normaltest(val_array)
            else:
                teststat = 100

            if teststat < 10:
                renormalizable += 1
                blurb (symbol, description,  "tails off:  well defined normal distro - presumably not modified from healthy; ", outf, cursor)
                continue

            if abs(skewness) <3:
                width = max - min
                if width < 20:
                    blurb (symbol, description, "tails off: weakly skewed - narrow", outf, cursor)
                    not_skewed += 1
                else:
                    blurb (symbol, description, "tails off: weakly skewed - wide", outf, cursor)
                    wide += 1
                continue

        # if we got to here, means we could not classify the distro
        [hist_numpy, low, binsize, extrapoints] = stats.histogram (val_array)
        histogram = hist_numpy.tolist()
        if is_bimodal(histogram):
            comment = "bimodal (?)"
        else:
            comment = "undecided"
        blurb (symbol, description, comment, outf, cursor)
        if False:
            i = 0
            bin_prev = low
            for val in histogram[:-1]:
                print >> outf," %5d %5d " % (int(bin_prev), int(bin_prev+binsize)),
                print >> outf," %5d " %  val
                i += 1
                bin_prev += binsize

            print >> outf,""
        other += 1


    print >> outf, "summary for parent directory", parent_dir
    print >> outf, "tot %d, weak %d,  weak with tail %d,   normal %d,   not_skewed %d, renormalizable: %d,  wide: %d,  other: %d" % \
          (tot,    weak,     weak_with_tail,      normal,       not_skewed,   renormalizable,  wide, other)




#########################################
def main():
    
    db     = connect_to_mysql()
    cursor = db.cursor()

    db_names  = ["BRCA", "COAD","GBM","KIRC","KIRP","LAML","LGG","LUAD","LUSC","OV","REA","UCEC"]
    #db_names  = ["KIRP"]

    switch_to_db (cursor, 'tcga_meta')

    to_file = False

    for db_name in db_names:

        print " ** ", db_name

        if to_file:
            if use_normal_tissue:
               outf = open("output/%s.rnaseq.normal.table"%db_name, "w")
            else:
                outf = open("output/%s.rnaseq.table"%db_name, "w")
        else:
            outf = sys.stdout

        print >> outf,  "%15s %4s  %6s  %10s     %10s   %6s  %5s     %8s     %5s  %5s "  \
                        % ('symbol', 'pts', 'min', 'max', 'mean', 'stdev', 'skew', 'kurt', 'left', 'right')

        table = 'gene_expression'

        db_dir  = '/Users/ivana/databases/TCGA/%s/Expression_Genes' % db_name

        parent_dirs = []
        files_in_data_dir = {}
        for parent_dir, subdir_list, file_list in os.walk(db_dir):
            data_files = filter (lambda x: x[-3:] == 'txt', file_list)
            if data_files:
                parent_dirs.append(parent_dir)
                files_in_data_dir[parent_dir] = data_files
        if not parent_dirs:
            print "no data found for db_name (?)"
        else:
            process_data_set (cursor, db_name, parent_dirs, files_in_data_dir, outf)

        if to_file: outf.close()


    cursor.close()
    db.close()


#########################################
if __name__ == '__main__':
    main()
