import os, subprocess
from icgc_utils.common_queries import *
# what is the difference between GRCh37 and hg19 ?
# https://www.biostars.org/p/123767/
# The genomic content for the two is identical, except for the mitochondrial contig
# annovar humandd firectory:
# annovar_downdb.log	       hg18_ensGeneMrna.fa  hg19_ensGeneMrna.fa	         hg19_MT_ensGeneMrna.fa  hg19_refGeneVersion.txt
# genometrax-sample-files-gff  hg18_ensGene.txt	    hg19_ensGene.txt	         hg19_MT_ensGene.txt     hg19_refGeneWithVerMrna.fa
# GRCh37_MT_ensGeneMrna.fa     hg18_refGeneMrna.fa  hg19_example_db_generic.txt  hg19_refGeneMrna.fa     hg19_refGeneWithVer.txt
# GRCh37_MT_ensGene.txt	       hg18_refGene.txt	    hg19_example_db_gff3.txt     hg19_refGene.txt


##################################
def assembly_name_translate(assembly, mitochondrial):
	if assembly=='hg18':
		return 'hg18'
	if assembly=='hg19':
		if mitochondrial:
			return 'hg19_MT'
		else:
			return 'hg19'
	if assembly=='GRCh37':
		if mitochondrial:
			return 'GRCh37_MT'
		else:
			return 'hg19'
	print("unrecognized assembly:", assembly)
	exit()


##################################
def run_annovar(avinput, assembly, table_name, mitochondrial=False):

	print ("running annovar ...")
	avoutname = "%s.%s_multianno.txt" % (table_name, assembly)
	# danger zone: this file batter does not exist if it is outdated
	if os.path.exists(avoutname) and os.path.getsize(avoutname)!=0:
		print("\t %s found"%avoutname)
		return avoutname

	translated_assembly_name = assembly_name_translate(assembly, mitochondrial)
	cmd  = "/home/ivana/third/annovar/table_annovar.pl %s " % avinput
	cmd += "/home/ivana/third/annovar/humandb/ -buildver %s -out %s " % (translated_assembly_name, table_name)
	cmd += " -protocol ensGene  -operation g  -nastring ."
	subprocess.call(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	# clean the junk
	cmd = "rm %s.ensGene.variant_function " % table_name
	cmd +="%s.ensGene.exonic_variant_function %s.ensGene.log" % (table_name, table_name)
	subprocess.call(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

	return avoutname


##################################
def parse_annovar_header_fields(annovar_input_line):
	header_fields = None
	if annovar_input_line[:3]=="Chr":  # header
		fields = annovar_input_line.rstrip().split('\t')
		header_fields = [f.split(".")[0].lower() for f in fields] # get rid of ".ensGene"
		# the header_fields now should be chr, start, end, ref, alt, func, gene, genedetail, exonicfunc, aachange
	return header_fields


#########################################
# annovar has a different idea about what  gene-relative location means
def gene_location_annovar_to_icgc(annovar_description):
	if annovar_description in ['upstream','downstream','intergenic']:
		return annovar_description
	if annovar_description in ['intronic', 'exonic','splicing', 'UTR5', 'UTR3']:
		return 'intragenic'
	if 'ncRNA' in annovar_description: # ncRNA_exonic, ncRNA_intronic, ncRNA_splicing
		return 'intragenic'
	return 'unk'


def annovar2gene_relative_string(annovar_named_field):
	gene_relative_string = None
	# icgc fields in locations_chrom_*: position, gene_relative, transcript relative
	# gene relative is in two annovar fields: func and gene
	gene_relative = {}
	for gene_id, funct in list(zip(annovar_named_field['gene'].split(";"),annovar_named_field['func'].split(";") )):
		if not gene_id in gene_relative: gene_relative[gene_id] = set()
		gene_relative[gene_id].add(gene_location_annovar_to_icgc(funct))
	for gene_id,annotation in gene_relative.items():
		if len(annotation)!=1:
			print("inconsistent or nonexistent annotation for", gene_id)
		gene_relative_string = ";".join(["%s:%s"%(gene_id,list(annotation)[0])])
	return gene_relative_string

#########################################
def parse_annovar_location_fields(cursor, annovar_named_field):
	# note: just write your own annotator - this is is  garbage

	gene_relative_string = annovar2gene_relative_string(annovar_named_field)
	# what to do with introns? can introns be  assigned to a transcript?
	# tcga thinks they can, while annovar things they are only assignable to a gene as a whole (I'd agree)
	# I am disregarding intronic mutations as zero impact
	# as an ad hoc measure, I will assign the annotation to the canonical transcript
	tr_relative = []
	# transcript relative contains a bit more info than I need - perhaps if I was doing it
	# again I would stick to annovar annotation, but not now
	# I would in any case change the names in the annovar header - they do not reflect the content
	# for these imbeciles the separator in the first field is ';' and in the second it is ','
	for trrel in annovar_named_field['genedetail'].split(";") + annovar_named_field['aachange'].split(","):
		fields = trrel.split(":") # the content of each field is not fixed
		[enst,cdna_change, aa_change, annotation] = [None]*4
		for field in fields:
			if field[:4] == "ENSG":
				# just drop, we already have that info (plus it is recoverable from ENST)
				pass
			elif field[:4] == "ENST":
				enst = field
			elif field[:2] == "c.":
				cdna_change = field
			elif field[:2] == "p.":
				aa_change = field
			elif field[:4] == "exon":
				annotation = 'exon'
			elif field[:4] == "UTR3":
				# there is no UTR annotation here - it is in he gene function or some such
				# so what - UTR3 is always the same ...perhaps
				annotation = 'UTR3'
			elif field[:4] == "UTR5":
				annotation = 'UTR5'

		if not enst:
			if aa_change:
				print("aa change without specified ENST (this should not have happened)")
				print(fields)
				exit()
			else:
				continue

		if cdna_change and ('-' in cdna_change or '+' in cdna_change):
			if 'UTR3' in annovar_named_field['func']:
				annotation = 'UTR3'
			elif 'UTR5' in annovar_named_field['func']:
				annotation = 'UTR5'
			elif not annotation:
				annotation='intron'
			elif annotation=='exon':
				annotation = 'splice_region'

		tr_relative.append("{}:{}".format(enst, annotation))

	if len(tr_relative)==0 and 'intron' in annovar_named_field['func']:
		# I really need to get rid of annovar here
		genes = annovar_named_field['gene']
		if genes:
			for gene in genes.split(";"):
				if gene  and 'ENSG' in gene:
					canonical = gene_stable_id_2_canonical_transcript_id(cursor, gene)
					if canonical: tr_relative.append("{}:{}".format(canonical, 'intron'))


	tr_relative_string  = ";".join(list(tr_relative))
	start = annovar_named_field['start']
	end = annovar_named_field['end']
	return [start, end, gene_relative_string, tr_relative_string]


########################################
def parse_annovar_line(cursor, header_fields, line):
	fields = line.rstrip().split('\t')
	annovar_named_field = dict(list(zip(header_fields,fields)))
	return parse_annovar_location_fields(cursor, annovar_named_field)