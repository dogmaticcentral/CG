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

from icgc_utils.common_queries  import  *
from icgc_utils.icgc import *
from config import Config


tcga_icgc_table_correspondence = {
"ACC_donor" : None,
"ALL_donor" : "ALL_donor",
"BLCA_donor": "BLCA_donor",
"BRCA_donor": "BRCA_donor",
"CESC_donor": "CESC_donor",
"CHOL_donor": None,
"COAD_donor": "COCA_donor",
"DLBC_donor": "DLBC_donor",
"ESCA_donor": "ESAD_donor",
"GBM_donor" : "GBM_donor",
"HNSC_donor": "HNSC_donor",
"KICH_donor": "KICH_donor",
"KIRC_donor": "KIRC_donor",
"KIRP_donor": "KIRP_donor",
"LAML_donor": "AML_donor",
"LGG_donor" : "LGG_donor",
"LIHC_donor": "LICA_donor",
"LUAD_donor": "LUAD_donor",
"LUSC_donor": "LUSC_donor",
"MESO_donor": None,
"OV_donor"  : "OV_donor",
"PAAD_donor": "PACA_donor",
"PCPG_donor": None,
"PRAD_donor": "PRAD_donor",
"READ_donor": "COCA_donor",
"SARC_donor": "SARC_donor",
"SKCM_donor": "MELA_donor",
"STAD_donor": "GACA_donor",
"TGCT_donor": None,
"THCA_donor": "THCA_donor",
"THYM_donor": None,
"UCEC_donor": "UCEC_donor",
"UCS_donor" : "UTCA_donor",
"UVM_donor" : None
}


#########################################
#########################################
def main():

	print("disabled")
	exit()

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name =  "icgc"
	switch_to_db(cursor, db_name)
	for tcga_donor_table, icgc_donor_table in tcga_icgc_table_correspondence.items():
		if icgc_donor_table: continue

		qry = "drop table " + tcga_donor_table
		search_db(cursor,qry, verbose=True)
		make_donors_table(cursor, db_name, tcga_donor_table)

		specimen_table_name =  tcga_donor_table.replace("donor", "specimen")
		qry = "drop table " + specimen_table_name
		search_db(cursor,qry, verbose=True)
		make_specimen_table(cursor, db_name, specimen_table_name)

		somatic_table_name = tcga_donor_table.replace("donor","simple_somatic")
		qry = "drop table " + somatic_table_name
		search_db(cursor,qry, verbose=True)
		make_variants_table(cursor, db_name, somatic_table_name)

	# When reading icgc donor and specimen tables from tsv,
	# we needed id *not* to be autoincremnet. Now we will be doing the
	# the inserts, so turn the autoincrementing back on.
	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_donor'"
	icgc_donor_tables = [field[0] for field in search_db(cursor,qry)]
	for icgc_donor_table in icgc_donor_tables:

		tumor = icgc_donor_table.split("_")[0]
		icgc_specimen_table = "%s_specimen"%tumor
		set_autoincrement(cursor, 'icgc',  icgc_donor_table, 'id')
		set_autoincrement(cursor, 'icgc',  icgc_specimen_table, 'id')

		# cleanup in case we were already fiddling with TCGA merge
		qry = "delete from {} where icgc_donor_id like 'TARGET%' or icgc_donor_id like 'TCGA%'"
		search_db(cursor,qry.format(icgc_donor_table), verbose=True)
		search_db(cursor,qry.format(icgc_specimen_table), verbose=True)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

