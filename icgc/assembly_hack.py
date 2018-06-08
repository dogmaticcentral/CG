#! /usr/bin/python
import subprocess
import time, re

from icgc_utils.common_queries  import  *
from icgc_utils.processes   import  *

tcga_icgc_table_correspondence = {
"ACC_somatic_mutations" :  "ACC_simple_somatic",
"ALL_somatic_mutations" :  "ALL_simple_somatic",
"BLCA_somatic_mutations": "BLCA_simple_somatic",
"BRCA_somatic_mutations": "BRCA_simple_somatic",
"CESC_somatic_mutations": "CESC_simple_somatic",
"CHOL_somatic_mutations": "CHOL_simple_somatic",
"COAD_somatic_mutations": "COCA_simple_somatic",
"DLBC_somatic_mutations": "DLBC_simple_somatic",
"ESCA_somatic_mutations": "ESAD_simple_somatic",
"GBM_somatic_mutations" :  "GBM_simple_somatic",
"HNSC_somatic_mutations": "HNSC_simple_somatic",
"KICH_somatic_mutations": "KICH_simple_somatic",
"KIRC_somatic_mutations": "KIRC_simple_somatic",
"KIRP_somatic_mutations": "KIRP_simple_somatic",
"LAML_somatic_mutations":  "AML_simple_somatic",
"LGG_somatic_mutations" :  "LGG_simple_somatic",
"LIHC_somatic_mutations": "LICA_simple_somatic",
"LUAD_somatic_mutations": "LUAD_simple_somatic",
"LUSC_somatic_mutations": "LUSC_simple_somatic",
"MESO_somatic_mutations": "MESO_simple_somatic",
"OV_somatic_mutations"  :   "OV_simple_somatic",
"PAAD_somatic_mutations": "PACA_simple_somatic",
"PCPG_somatic_mutations": "PCPG_simple_somatic",
"PRAD_somatic_mutations": "PRAD_simple_somatic",
"READ_somatic_mutations": "COCA_simple_somatic",
"SARC_somatic_mutations": "SARC_simple_somatic",
"SKCM_somatic_mutations": "MELA_simple_somatic",
"STAD_somatic_mutations": "GACA_simple_somatic",
"TGCT_somatic_mutations": "TGCT_simple_somatic",
"THCA_somatic_mutations": "THCA_simple_somatic",
"THYM_somatic_mutations": "THYM_simple_somatic",
"UCEC_somatic_mutations": "UCEC_simple_somatic",
"UCS_somatic_mutations" : "UTCA_simple_somatic",
"UVM_somatic_mutations" :  "UVM_simple_somatic"
}

def invert_table_correpondence():

	icgc_to_tcga = {}

	for tcga, icgc in tcga_icgc_table_correspondence.iteritems():
		if not icgc_to_tcga.has_key(icgc): icgc_to_tcga[icgc] = []
		icgc_to_tcga[icgc].append(tcga)

	return icgc_to_tcga

#########################################
def main():


	icgc_to_tcga = invert_table_correpondence()

	# divide by cancer types, because I have duplicates within each cancer type
	# that I'll resolve as I go, but I do not want the threads competing)
	db     = connect_to_mysql()
	cursor = db.cursor()

	qry  = "select table_name from information_schema.tables "
	qry += "where table_schema='icgc' and table_name like '%_simple_somatic'"
	#qry += "where table_schema='icgc' and table_name like 'mutations_chrom_%'"
	tables = [field[0] for field in search_db(cursor,qry)]

	switch_to_db(cursor,'icgc')

	for table in tables:
		print table
		qry  = "select distinct(icgc_mutation_id) from %s " % table
		ret = search_db(cursor,qry)
		if not ret: continue
		for row in ret:
			icgc_mutation_id = row[0]
			assembly = None
			qry  = "select icgc_donor_id, submitted_sample_id, chromosome from %s " % table
			qry += "where icgc_mutation_id='%s' limit 1 " % icgc_mutation_id
			icgc_donor_id, submitted_sample_id, chromosome = search_db(cursor,qry)[0]
			if icgc_donor_id[2] == 'T': # this came from tcga - the assembly field may be empty
				switch_to_db(cursor,'tcga')
				for tcga_table in icgc_to_tcga[table]:
					meta_table = tcga_table.split("_")[0] + "_mutations_meta"
					qry2  = "select  distinct(m.assembly) from %s t, %s m " %(tcga_table, meta_table)
					qry2 += "where t.tumor_sample_barcode='%s' " % submitted_sample_id
					qry2 += "and t.meta_info_id=m.id"
					ret = search_db(cursor,qry2)
					if ret:
						assembly = ret[0][0]
						break
				switch_to_db(cursor,'icgc')
			else:
				qry2  = "select assembly from %s_temp " % table
				qry2 += "where icgc_mutation_id='%s' limit 1 " % icgc_mutation_id

				ret = search_db(cursor,qry2)
				if ret:
					assembly = ret[0][0]
			if assembly:
				qry2  = "update mutations_chrom_%s " % chromosome
				qry2 += "set assembly='%s' " % assembly
				qry2 += "where icgc_mutation_id='%s' " % icgc_mutation_id
				search_db(cursor,qry2)
	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()
