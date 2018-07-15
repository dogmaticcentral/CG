#!/usr/bin/python -u

import matplotlib.pyplot as plt


#########################################
def main():

	# extras = ['RPL5', 'RPL11', 'RPL22',
	#           'TP53','TTN', 'PTEN', 'SMAD4','SETD2',
	#           'RB1', 'VHL', 'PIK3CA', 'IDH1',
	#           'WDR41',
	#           'RABAC1','BEX5', 'DEFB116', 'NODAL','DDX3X','RAD9B', 'LUC7L2', 'FZD10',
	#            'ATM' ,  'BRCA2', 'NF1', 'CASP8', 'TSGA10',
	#           'BBX','B2M',
	#           'SOS1','FBXW7','GNAS',
	#           'PIK3CA', 'ZC3H13', 'APC',
	#           'BRCA1','MTOR', 'KRAS' , 'MDM2',
	#           'TBL1XR1', 'UBC', 'HSPA1L', 'SP8',
	#           'VHL', 'BRAF', 'MUC2', 'MUC4', 'NPM1']

	zoom = False

	if zoom:
		extras = ['RPL5', 'RPL11', 'RPL22',
		          'TP53','TTN', 'PTEN', 'SMAD4','SETD2',
		          'RB1', 'VHL',
		          'WDR41',
		          'RABAC1','BEX5', 'DEFB116', 'RAD9B', 'LUC7L2', 'FZD10',
		          'NPM1' ,  'BRCA2',  'CASP8',  'B2M',
		          'SOS1',
		          'PIK3CA', 'APC',
		          'BRCA1','MTOR', 'KRAS' , 'MDM2',
		          'UBC', 'HSPA1L', 'SP8',
		          'VHL', 'BRAF', 'MUC2', 'MUC4']
	else:
		extras = ['BEX5','PTEN','TP53','TTN', 'HSPA1L','MUC16']


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
	if zoom:
		plt.xlim(0.06,0.3)
		plt.ylim(50,1000)
		plt.xlabel('SILENT MUTATION FRACTION', fontsize=16)
		plt.ylabel('TOTAL NUMBER OF MUTATIONS', fontsize=16)
		plt.scatter(x,y,s=0.5, c='#17aad3')
		for extra in extras:
			plt.plot(coords[extra][0], coords[extra][1],'ro',mfc='none')
			annotcoords = [coords[extra][0] + 0.003, coords[extra][1]+5]
			plt.annotate(extra, annotcoords, color="r", fontsize=12) #, fontweight='bold')

	else: # full pic
		plt.xlim(0.025,1)
		plt.ylim(10,11000)
		plt.scatter(x,y,s=1,c='#17aad3')
		for extra in extras:
			plt.plot(coords[extra][0], coords[extra][1],'ro',mfc='none')
			plt.annotate(extra, coords[extra], color="r", fontsize=15)#, fontweight='bold')


	plt.savefig("zoom.svg", format="svg")
	plt.show()

#########################################
if __name__ == '__main__':
	main()


'''
% HSPA1L
% 
https://www.nature.com/articles/onc2017263?WT.feed_name=subjects_cancer

% BEX5 > poorly described, brain expressed
% PTEN well-stuided tumor suppressor 
#https://www.annualreviews.org/doi/abs/10.1146/annurev.pathol.4.110807.092311

'''