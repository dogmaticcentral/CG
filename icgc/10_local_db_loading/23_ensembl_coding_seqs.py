#! /usr/bin/python3
#
# This source code is part of icgc, an ICGC processing pipeline.
#
# Icgc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Icgc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
#
# Contact: ivana.mihalek@gmail.com
#

# from ensembl biomart:  human sequences coding, stable transcript id in the header

from config import Config
from icgc_utils.mysql import *

#########################################
# make one-on-one ENST to ENSG translation table
def make_ensembl_coding_seqs_table(cursor, db_name, table_name):

	if check_table_exists (cursor, db_name, table_name): return
	switch_to_db (cursor, db_name)
	qry = "  CREATE TABLE  %s (" % table_name
	qry += " transcript_id VARCHAR(20) NOT NULL ,"
	qry += " sequence text,"
	qry += "	 PRIMARY KEY (transcript_id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)
	return

def store_seq (cursor, table_name, transcript_id, sequence):
	if transcript_id and len(sequence)>0:
		qry = "insert into %s (transcript_id, sequence) values ('%s','%s')" % (table_name, transcript_id,  "".join(sequence))
		search_db(cursor, qry, verbose=False)


#########################################
#########################################
def main():

	ens_seq_file = "/storage/databases/ensembl-96/human/fasta/ensembl_human_canontrans_coding_seqs.fasta"

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()
	db_name = "icgc"
	table_name = 'ensembl_coding_seqs'

	make_ensembl_coding_seqs_table(cursor, db_name, table_name)

	inf = open(ens_seq_file, "r")
	transcript_id = None
	sequence=[]
	for line in inf:
		if line[0]=='>':
			store_seq (cursor, table_name, transcript_id, sequence)
			transcript_id = line.strip()[1:]
			sequence=[]
		else:
			sequence.append(line.strip())

	store_seq (cursor, table_name, transcript_id, sequence)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()
