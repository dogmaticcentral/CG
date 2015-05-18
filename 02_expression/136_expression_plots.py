#!/usr/bin/python

import json
from pprint import pprint

# gene_expression.json came from 136 script, I just deleted
# var name I needed for incorporation into the webpage

#########################################
def main():
    
    with open('gene_express.json') as data_file:    
        data = json.load(data_file)
    
    data_per_gene = {}
    for tumor, per_gene in data.iteritems():
        for gene, histogram in per_gene.iteritems():
            if not data_per_gene.has_key(gene): data_per_gene[gene] = {}
            data_per_gene[gene][tumor] = histogram

    for gene, per_tumor in data_per_gene.iteritems():
        print ###########################
        print gene
        for tumor, histogram in data_per_gene[gene].iteritems():
            print "%1s" %  tumor
            new_hist = {}
            for item in histogram['hist'].iteritems():
                new_hist[ int(item[0])] = item[1]
            for i in range (-8,9):
                print "%3d   %8.2f " % (i, float(new_hist[i]))



#########################################
if __name__ == '__main__':
    main()
