
from icgc_utils.mysql   import  *

########################################
def get_donors(cursor, table):
	qry  = "select distinct(icgc_donor_id) from %s " % table
	return [ret[0] for ret in search_db(cursor,qry)]

def get_mutations(cursor, table):
	qry = "select  distinct(icgc_mutation_id)  from %s " % table
	return [ret[0] for ret in search_db(cursor,qry)]

def get_specimens_from_donor(cursor, table, icgc_donor_id):
	qry = "select  distinct(icgc_specimen_id)  from %s " % table
	qry += "where icgc_donor_id = '%s'" % icgc_donor_id
	return [r[0] for r in search_db(cursor,qry)]

def get_specimen_type(cursor, table, spec_ids):
	specimen_type = {}
	for spec_id in spec_ids:
		qry = " select specimen_type from %s " % table.replace("simple_somatic","specimen")
		qry += "where icgc_specimen_id = '%s'" % spec_id
		specimen_type[spec_id] = search_db(cursor,qry)[0][0]
	return specimen_type

def get_mutations_from_donor(cursor, table, icgc_donor_id):
	qry = "select  distinct(icgc_mutation_id)  from %s " % table
	qry += "where icgc_donor_id = '%s'" % icgc_donor_id
	return [r[0] for r in search_db(cursor,qry)]

def mutation_provenance(cursor, table, icgc_donor_id, icgc_mutation_id):
	qry = "select  distinct(icgc_specimen_id)  from %s " % table
	qry += "where icgc_donor_id='%s' and icgc_mutation_id='%s'" % (icgc_donor_id, icgc_mutation_id)
	return [r[0] for r in search_db(cursor,qry)]

