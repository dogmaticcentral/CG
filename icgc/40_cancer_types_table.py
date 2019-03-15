#!/usr/bin/python

import os

def main():

	country_codes  = ['US', 'CA', 'UK', 'FR', 'CN', 'KR', 'SG','EU','JP', 'AU','BR','IN','SA','DE','ES']

	base_dir = "/storage/databases/icgc"
	meta = base_dir+"/meta.tsv"
	data_dir =  base_dir+"/v27"
	for dependency in [base_dir, data_dir, meta]:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()

	cancer_name = {}
	inf = open(meta,"r")
	for line in  inf:
		fields = line.strip().split("-")
		cancer_type = fields[0]
		description = "-".join(fields[1:-1])
		# get rid of the country code
		description = (" ".join(description.split("\t")[1:])).strip()
		if cancer_type not in cancer_name: cancer_name[cancer_type] =set ([])
		cancer_name[cancer_type].add(description)
	inf.close()

	cancer_types = [nm for nm in os.listdir(data_dir) if nm[-4:] not in ['json']]

	for cancer_type in sorted(cancer_types):
		subtypes = [nm for nm in os.listdir("{}/{}".format(data_dir,cancer_type)) if nm not in country_codes]
		description = set([])
		for ctype in [cancer_type]+subtypes:
			description |= cancer_name.get(ctype)
		print("%s\t%s\t%s" % (cancer_type, ", ".join(subtypes), "; ".join(description)))


	return

#########################################
if __name__ == '__main__':
	main()
