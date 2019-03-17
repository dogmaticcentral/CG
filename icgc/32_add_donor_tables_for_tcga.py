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
# icgc_donor_id,submitted_donor_id,donor_sex,donor_diagnosis_icd10
def make_donors_table(cursor, db_name, donor_table):

	switch_to_db (cursor, db_name)
	qry = ""
	qry += "  CREATE TABLE  %s (" % donor_table
	qry += "     id INT NOT NULL, "
	qry += "  	 icgc_donor_id VARCHAR (20) NOT NULL, "
	qry += "     submitted_donor_id VARCHAR (50) NOT NULL, "
	qry += "	 donor_sex VARCHAR (10), "
	qry += "	 donor_diagnosis_icd10  VARCHAR (20), "

	qry += "	 PRIMARY KEY (id) "
	qry += ") ENGINE=MyISAM"

	rows = search_db(cursor, qry)
	print(qry)
	print(rows)


#########################################
#########################################
def main():

	db     = connect_to_mysql(Config.mysql_conf_file)
	cursor = db.cursor()

	db_name =  "icgc"
	for tcga_table, icgc_table in tcga_icgc_table_correspondence.items():
		if not icgc_table and not check_table_exists(cursor, "icgc", tcga_table):
			make_donors_table(cursor, db_name, tcga_table)

	cursor.close()
	db.close()

#########################################
if __name__ == '__main__':
	main()

