#!/usr/bin/python

import os

def main():

	country_codes  = ['US', 'CA', 'UK', 'FR', 'CN', 'KR', 'SG','EU','JP', 'AU','BR','IN','SA','DE','ES']

	data_dir = "/storage/databases/icgc"
	meta = data_dir+"/meta.txt"
	for dependency in [data_dir, meta]:
		if not os.path.exists(dependency):
			print dependency, "not found"
			exit()

	cancer_name = {}
	inf = open(meta,"r")
	for line in  inf:
		fields = line.strip().split("-")
		cancer_type = fields[0]
		description = "-".join(fields[1:-1])
		# get rid of the country code
		description = (" ".join(description.split(" ")[1:])).strip()
		if not cancer_name.has_key(cancer_type): cancer_name[cancer_type] =set ([])
		cancer_name[cancer_type].add(description)
	inf.close()

	cancer_types = [nm for nm in os.listdir(data_dir) if nm[-3:]!='txt']

	for cancer_type in sorted(cancer_types):

		subtypes = [nm for nm in os.listdir("{}/{}".format(data_dir,cancer_type)) if nm not in country_codes]
		description = set([])
		for ctype in [cancer_type]+subtypes:
			description |= cancer_name.get(ctype)
		print "%s\t%s\t%s" % (cancer_type, ", ".join(subtypes), "; ".join(description))


	return

#########################################
if __name__ == '__main__':
	main()
