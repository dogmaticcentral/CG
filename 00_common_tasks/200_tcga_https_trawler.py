#!/usr/bin/python -u
#
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
#
#  requests lib http://docs.python-requests.org/en/latest/
#
# html parser: http://www.crummy.com/software/BeautifulSoup/, 
# see also http://stackoverflow.com/questions/11709079/parsing-html-python
#
#
import requests
#from BeautifulSoup import BeautifulSoup
from bs4 import BeautifulSoup

# Set security warning to always go off by default.
import warnings

not_interesting = ['bio', 'pathology_reports', 'diagnostic_images', 'tissue_images',
                  'mirnaseq', 'miRNASeq', 'mirna',  'transcriptome',
                  'methylation', 'protein_exp',  'mutations_protected', 
                   'cdna', 'cna', 'bisulfiteseq', 'rnaseqv2','totalrnaseqv2',
                   'microsat_i', 'exon', 'tracerel']

target_set = 'mutations'

# this script can be used to locate data files to be downloaded by 201_tcga_http_dnwld.data
############################
# to download rnaseq based expression data use 02_expression/100_rnaseq_http_dwnld.py
############################
# apparently 'transcriptome' is a code word for array based gene expression in TCGA,
# and the method has been somewhat discredited lately
############################


if  target_set == 'mutations': # we are looking for somatic mutations
    not_interesting += ['snp', 'rnaseq']
elif target_set == 'cnv':
    not_interesting += ['mutations', 'rnaseq']
elif target_set == 'expression':
    not_interesting += ['mutations', 'snp']
else:
    print "target set", target_set, " not recognized"



##########################################
def get_page_text (httpaddr):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        http_req  = requests.get(httpaddr)
    return http_req.text

#########################################
def recursive_descent (name_pieces):


    tcga_http = "".join(name_pieces)
    ret_text = get_page_text(tcga_http)

    soup = BeautifulSoup(ret_text)
    #print dir(soup)
    for link in soup.findAll('a'):
        if 'Parent Directory' in link.getText(): continue
        new_name_piece = link.get('href')
        if 'http' in  new_name_piece: continue

        # here we are picking the type of data we are looking for
        # the whole databse is such a mess I do not know how to look for
        # the data type I need in a more robust way
        #print new_name_piece[:-1], "not interesting?",  new_name_piece[:-1] in not_interesting
        if new_name_piece[:-1] in not_interesting: continue
        if target_set == 'mutations':
            if 'Level_1' in new_name_piece: continue
            if 'Level_3' in new_name_piece: continue
        else:
            if 'Level_1' in new_name_piece: continue
            if 'Level_2' in new_name_piece: continue

        if new_name_piece[-1] != '/' :

            if target_set == 'expression':
                if 'trimmed.annotated.gene.quantification.txt' in new_name_piece:
                    print "".join(name_pieces[1:]), ":  ",
                    print new_name_piece


            elif 'tar.gz' in new_name_piece[-6:] and not 'mage-tab' in new_name_piece:
                print "".join(name_pieces[1:]), ":  ",
                print new_name_piece


            continue

        name_pieces.append(new_name_piece)
        recursive_descent (name_pieces)
        name_pieces.pop()

    return

#########################################
def main():
    
    main_url = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/'
    # recursive descent is well, a recursion, which we seed with the core url string
    name_pieces = [main_url]
    recursive_descent (name_pieces)
    

#########################################
if __name__ == '__main__':
    main()
