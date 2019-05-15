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
def run_annovar(avinput, assembly, out_name_root, mitochondrial=False):

	avoutname = ".".join(avinput.split(".")[:-1]+["avout"])
	# danger zone: this file better does not exist if it is outdated
	if os.path.exists(avoutname) and os.path.getsize(avoutname)!=0:
		print("\t\t %s found"%avoutname)
		return avoutname
	cmd = "rm -f *multianno.txt"
	subprocess.call(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

	# the actual run
	translated_assembly_name = assembly_name_translate(assembly, mitochondrial)
	cmd  = "/home/ivana/third/annovar/table_annovar.pl %s " % avinput
	cmd += "/home/ivana/third/annovar/humandb/ -buildver %s -out %s " % (translated_assembly_name, out_name_root)
	cmd += " -protocol ensGene  -operation g  -nastring ."
	subprocess.call(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

	# move the output to the name we want (annovar insists on chopping down the name to the first few tokens)
	cmd = "mv *multianno.txt " + avoutname
	subprocess.call(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

	# clean the junk
	cmd = "rm %s.ensGene.variant_function " % out_name_root
	cmd += "%s.ensGene.exonic_variant_function %s.ensGene.log" % (out_name_root, out_name_root)
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
					if canonical:
						tr_relative.append("{}:{}".format(canonical, 'intron'))
					#else:
					#	gene_stable_id_2_canonical_transcript_id(cursor, gene, verbose=True)


	tr_relative_string  = ";".join(list(tr_relative))
	start = annovar_named_field['start']
	end = annovar_named_field['end']
	return [start, end, gene_relative_string, tr_relative_string]


########################################
def get_icgc_mutation_type(start, end, ref, alt):
	frame_consequences = set([])
	if ref=='-':
		if len(alt)%3==0:
			frame_consequences.add('inframe')
		else:
			frame_consequences.add('frameshift')
		return 'insertion', frame_consequences
	if alt=='-':
		if len(ref)%3==0:
			frame_consequences.add('inframe')
		else:
			frame_consequences.add('frameshift')
		return 'deletion', frame_consequences
	if start==end:
		return 'single', frame_consequences
	return 'multiple', frame_consequences

#########################################
def parse_annovar_mutation_fields(cursor, annovar_named_field):

	# I should probably get rid of mutation_type columns - it is completely derivable from other columns
	# note though that it might complement consequence notation which only says inframe or frameshift
	tr_relative = []
	aa_mutation = []
	consequences = set([])
	mutation_type, frame_consequences = get_icgc_mutation_type(annovar_named_field['start'], annovar_named_field['end'],
	                                                           annovar_named_field['ref'], annovar_named_field['alt'])
	# icgc fields in locations_chrom_*: position, gene_relative, transcript relative
	# gene relative is in two annovar fields: func and gene
	gene_relative = dict(list(zip(annovar_named_field['gene'].split(";"),annovar_named_field['func'].split(";") )))
	gene_relative_string = ";".join(["%s:%s"%(k,gene_location_annovar_to_icgc(v))
	                                 for k,v in gene_relative.items()])
	# what to do with introns? can introns be  assigned to a transcript?
	# tcga thinks they can, while annovar things they are only assignable to a gene as a whole (I'd agree)
	# I am disregarding intronic mutations as zero impact
	# as an ad hoc measure, I will assign the annotation to the canonical transcript
	if 'intronic' in annovar_named_field['func']:
		consequences.add('intronic')
		for gene_stable_id, annot in gene_relative.items():
			if annot != 'intronic': continue
			canonical_transcript_id = gene_stable_id_2_canonical_transcript_id(cursor, gene_stable_id)
			tr_relative.append("{}:{}".format(canonical_transcript_id, 'intronic'))
	# transcript relative contains a bit more info than I need - perhaps if I was doing it
	# again I would stick to annovar annotation, but not now
	# I would in any case change the names in the annovar header - they do not reflect the content
	for trrel in annovar_named_field['aachange'].split(","):
		fields = trrel.split(":") # the content of each field is not fixed
		[enst,cdna_change, aa_change, annotation] = [None]*4
		for field in fields:
			if field[:4] == "ENSG":
				# just drop, we already have that info
				pass
			elif field[:4] == "ENST":
				enst = field
			elif field[:2] == "c.":
				cdna_change = field
				pass
			elif field[:2] == "p.":
				aa_change = field
			elif field[:4] == "exon":
				annotation = 'exon'
			elif field[:4] == "UTR3":
				annotation = 'UTR3'
			elif field[:4] == "UTR5":
				annotation = 'UTR5'

		if not enst:
			if aa_change:
				print("aa change without specified ENST", fields)
				exit()
			else:
				continue

		if annotation=='exon' and cdna_change and ('-' in cdna_change or '+' in cdna_change):
			annotation = 'splice_region'
			consequences.add('splice_region')
		tr_relative.append("{}:{}".format(enst, annotation))

		if aa_change:
			aa_mutation.append("{}:{}".format(enst, aa_change[2:]))
			if aa_change[-1] in ['X','x','*']:
				consequences.add('stop_gained')
			elif aa_change[-2:].lower()=='fs':
				consequences.add('frameshift')
			elif aa_change[-2:].lower() in ['ins','del']:
				pass
			elif aa_change[2]==aa_change[-1]:
				consequences.add('synonymous')
			else:
				consequences.add('missense')

	tr_relative_string  = ";".join(list(tr_relative))
	if 'exon' in tr_relative_string: consequences |= frame_consequences
	consequences_string = ";".join(list(consequences))
	aa_mutation_string = ";".join(list(aa_mutation))

	return [mutation_type, consequences_string, aa_mutation_string, gene_relative_string, tr_relative_string]

########################################
def parse_annovar_line(cursor, header_fields, line):
	fields = line.rstrip().split('\t')
	named_field = dict(list(zip(header_fields,fields)))
	dictkey = "_".join([named_field['chr'], named_field['start'], named_field['end'], named_field['ref'], named_field['alt']])
	return dictkey, parse_annovar_location_fields(cursor, named_field), parse_annovar_mutation_fields(cursor, named_field)

########################################
def parse_annovar(cursor, avout):
	annotation = {}
	inf = open(avout, "r")
	header_fields = None
	for line in inf:
		if not header_fields:
			header_fields = parse_annovar_header_fields(line)
			if not header_fields: # header must be the first line
				print("header not found in", avout)
				exit()
		else:
			dictkey, loc, mut = parse_annovar_line(cursor, header_fields, line)
			annotation[dictkey]= {'loc':loc, 'mut': mut}

	return annotation

