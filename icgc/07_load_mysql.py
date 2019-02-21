#!/usr/bin/python


from icgc_utils.mysql   import  *


#########################################
#########################################
def main():
	print("disabled")
	return

	homedir = "/data/icgc"
	cancer_types = []
	for name in os.listdir(homedir):
		if os.path.isdir("/".join([homedir,name])): cancer_types.append(name)

	db     = connect_to_mysql()
	cursor = db.cursor()

	db_name = "icgc"
	switch_to_db(cursor, db_name)
	for ct in cancer_types:
		for filetype in ['simple_somatic_temp']:
		#for filetype in ["donor", "specimen"]:
			table = "{}_{}".format(ct, filetype)
			qry = "load data local infile 'tsvs/%s.tsv' into table %s" % (table,table)
			search_db(cursor,qry,verbose=True)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()

