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
# Set security warning to always go off by default.
import json
import urllib, os
#########################################
def check_and_make(path):
    if not os.path.exists(path):
        print path, "not found, making", path
        os.makedirs(path)
    elif not os.path.isdir(path):
        print ">>> PROBLEM <<<<< ", path,  "is not directory ?"
        exit(1)


########################################
def main():

    main_url = 'https://gdc-api.nci.nih.gov/files'
    filter   = {
        "op" : "=",
        "content":{
            "field": "data_type",
            "value": ["Masked Somatic Mutation"]
        }
    }
    # there are fewere htan 200 cancer types currently, so I should be ok
    params = {'filters':json.dumps(filter), "from":1, "size":200}
    json_hash = requests.get(main_url, params = params).json()
    hits_array  = json_hash["data"]["hits"]
    #data_types = set([ hit["data_type"] for hit in hits_array])
    for hit in hits_array:
        filename = hit["file_name"]
        cancer_type = filename.split('.')[1]
        path = "/mnt/databases/TCGA/%s/Somatic_Mutations" % cancer_type

        print "downloading ", filename, " to ",  path
        check_and_make(path)
        full_url = "https://gdc-api.nci.nih.gov/data/%s" %  hit["file_id"]

        dwnldfile = urllib.URLopener()
        dwnldfile.retrieve(full_url, path + "/" + filename)
# unzip: find /mnt/databases/TCGA -name "*maf.gz" | xargs gunzip

#########################################
if __name__ == '__main__':
    main()

# first attempts ...
# it loks like I am looking for Masked Somatic Mutation cases:
# there the true identity of the carrier is protected by
# 'masking' the original allele by replacing it with the reference genome variant
  #  from=1&size=2 means return hits from 1 to 2
    #  main_url = 'https://gdc-api.nci.nih.gov/files?from=1&size=20000'
    #  ret_text = get_page_text(main_url)
    #  json_hash = json.loads(ret_text)
    #  hits_array  = json_hash["data"]["hits"]
    #data_types = set([ hit["data_type"] for hit in hits_array])
    #for data_type in data_types:
    #    print data_type
    #  for hit in hits_array:
     #     if not "maf" in  hit["file_name"]: continue
    #      if "protected"  in  hit["file_name"]: continue
   #       print hit["data_type"],  hit["file_name"]
    # recursive descent is well, a recursion, which we seed with the core url string
    # name_pieces = [main_url]
    # recursive_descent (name_pieces)

