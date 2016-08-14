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
from mysql import *
import re
import urllib2
from bs4 import BeautifulSoup

###############################################################################################
def segment_from_das(assembly, chrom, start, end, verbose=False):
    # das counts from 1, not from 0
    das_request = "http://genome.ucsc.edu/cgi-bin/das/%s/" % assembly
    das_request += "dna?segment=%s:%s,%s" % (chrom, start, end)
    if verbose:
        print das_request
    response = urllib2.urlopen(das_request)
    html = response.read()
    soup = BeautifulSoup(html, 'html.parser')
    if not soup or not soup.dna or not soup.dna.string: return None
    stripped = soup.dna.string.strip().upper()
    return re.sub('\s','', stripped) # get rid of all space characters

###############################################################################################
def construct_refGene_qry (chrom, fields):
    qry = "select "
    first = True
    for field in fields:
        if first:
            first = False
        else:
            qry += ", "
        qry += "g." + field
    qry += " "
    qry += "from refGene as g, gbCdnaInfo as i, refSeqStatus as s "
    qry += "where g.chrom='%s' " % chrom
    qry += "and i.acc=g.name  and i.type='mRNA' and i.mol='mRNA' "
    qry += "and s.mrnaAcc=g.name and s.status in ('Validated', 'Reviewed')"

    return qry

