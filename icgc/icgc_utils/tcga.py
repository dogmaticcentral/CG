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

import os
from icgc_utils.utils import *
from icgc_utils.mysql import *

#########################################
def get_position_translation(cursor, tcga_somatic_table, ref_assembly):
	standard_chroms = [str(i) for i in range(1,23)] + ['X','Y']
	meta_table_name = tcga_somatic_table.split("_")[0] + "_mutations_meta"

	# get the info that annovar needs
	qry  = "select s.tumor_sample_barcode, s.chromosome, s.start_position, s.end_position, "
	qry += "s.reference_allele, s.tumor_seq_allele1, s.tumor_seq_allele2, m.assembly  "
	qry += "from tcga.%s s, tcga.%s m " % (tcga_somatic_table, meta_table_name)
	qry += "where s.meta_info_id=m.id"
	rows = search_db(cursor, qry)

	position_translation = {}
	positions = {}
	for row in rows:
		(tumor_sample_barcode, chromosome, start_position, end_position, reference_allele,
			tumor_seq_allele1, tumor_seq_allele2, assembly) = \
			[str(entry, 'utf-8') if type(entry)==bytes else str(entry) for entry in row]
		if not chromosome in standard_chroms: continue
		if not assembly in positions: positions[assembly] = {}
		if not chromosome in positions[assembly]: positions[assembly][chromosome] = set()
		positions[assembly][chromosome] |= {start_position, end_position} # set literal

	for assembly, pos_per_assembly in positions.items():
		if not assembly in position_translation: position_translation[assembly] = {}
		for chromosome, pos in pos_per_assembly.items():
			# some positions may remain untranslated - translate_positions() will issue warning
			position_translation[assembly][chromosome] = translate_positions(pos, chromosome, assembly, ref_assembly)

	return position_translation

#########################################
def tcga_sample2tcga_donor(s):
	return "-".join(s.split("-")[:3])

#########################################
def tcga_sample2tcga_specimen(sample_barcode):
	pieces = sample_barcode.split("-")
	sample_id = "-".join(pieces[:4])
	# the last character here is a "vial"
	# is "vial' the same as "sample" in ICGC parlance?
	return sample_id

#########################################
def specimen_type_from_TCGA_barcode(sample_barcode):
	# we want to translate this to something similar to what ICGC is using
	# roughly: Normal, Primary, Metastatic, Recurrent, Cell_line
	tcga_sample_code = sample_barcode.split("-")[3][:2]
	if tcga_sample_code in ['01','03','05', '09']:
		return 'Primary'

	elif tcga_sample_code in ['02','04','40']:
		return 'Recurrent'

	elif tcga_sample_code in ['06','07']:
		return 'Metastatic'

	elif tcga_sample_code in ['10','11','12','13','14']:
		return 'Normal'

	elif tcga_sample_code in ['08','50']:
		return 'Cell_line'

	return "Other"



