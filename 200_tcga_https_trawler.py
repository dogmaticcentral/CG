#!/usr/bin/python
#
#  requests lib http://docs.python-requests.org/en/latest/
#
# html parser: http://www.crummy.com/software/BeautifulSoup/, 
# see also http://stackoverflow.com/questions/11709079/parsing-html-python
#
#
import requests
from BeautifulSoup import BeautifulSoup

# Set security warning to always go off by default.
import warnings

def get_page_text (httpaddr):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        http_req  = requests.get(httpaddr)
    return http_req.text


not_interesting = ['bio', 'pathology_reports', 'diagnostic_images', 'tissue_images',
                  'mirnaseq', 'miRNASeq', 'mirna', 
                  'methylation', 'protein_exp', 'mutations', 'mutations_protected', 
                  'rnaseqv2', 'cdna', 'cna', 'bisulfiteseq', 'rnaseq', 'totalrnaseqv2',
                  'transcriptome', 'microsat_i', 'exon', 'tracerel']

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
        if new_name_piece[-1] != '/': 
            if 'tar.gz' in new_name_piece[-6:] and not 'mage-tab' in new_name_piece:
                print "".join(name_pieces[1:]), ":  ", 
                print new_name_piece
            continue

        # here we are picking the type of data we are looking for
        # the whole databse is such a mess I do not know how to look for
        # the data type I need in a more robust way
        if new_name_piece[:-1] in not_interesting: continue
            
        name_pieces.append(new_name_piece)
        recursive_descent (name_pieces)
        name_pieces.pop()

    return

#########################################
def main():
    
    name_pieces = ['https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/']
    recursive_descent (name_pieces)
    

#########################################
if __name__ == '__main__':
    main()
