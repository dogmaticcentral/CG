#!/usr/bin/python -u

import matplotlib.pyplot as plt


#########################################
def main():

	extras = ['RPL5', 'RPL11', 'RPL22',
	          'TP53','TTN', 'PTEN', 'SMAD4','SETD2',
	          'RB1', 'VHL', 'PIK3CA', 'IDH1',
	          'WDR41',
	          'RABAC1','BEX5', 'DEFB116', 'NODAL','DDX3X','RAD9B', 'LUC7L2', 'FZD10',
	           'ATM' ,  'BRCA2', 'NF1', 'CASP8', 'TSGA10',
	          'BBX','B2M',
	          'SOS1','FBXW7','GNAS',
	          'PIK3CA', 'ZC3H13', 'APC',
	          'BRCA1','MTOR', 'KRAS' , 'MDM2',
	          'TBL1XR1', 'UBC', 'HSPA1L', 'SP8',
	          'VHL', 'BRAF', 'MUC2', 'MUC4', 'NPM1']

	inf = open("silent_ratio.tsv","r")
	x = []
	y = []
	coords = {}
	for line in inf:
		[gene, ratio, totct] = line.split('\t')
		ratio = float(ratio)
		totct = int(totct)
		x.append(ratio)
		y.append(totct)
		for extra in extras:
			if gene == extra: coords[extra] = [ratio,totct]


	# s is marker size in points**2
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(0.025,1)
	plt.ylim(10,11000)
	plt.scatter(x,y,s=4)
	for extra in extras:
		plt.plot(coords[extra][0], coords[extra][1],'ro',mfc='none')
		plt.annotate(extra, coords[extra], color="r")


	plt.show()

#########################################
if __name__ == '__main__':
	main()
