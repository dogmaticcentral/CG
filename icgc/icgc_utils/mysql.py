import MySQLdb, sys, os
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
#
import traceback
#
# This source code is part of tcga, a TCGA processing pipeline, written by Ivana Mihalek.
# Copyright (C) 2014-2016 Ivana Mihalek.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see<http://www.gnu.org/licenses/>.
# 
# Contact: ivana.mihalek@gmail.com
#

########
def set_autoincrement(cursor, db_name, table_name, column):
	# what should be the next value for column
	qry = "select max({}) from {}.{}".format(column, db_name, table_name)
	new_column_value = 1
	ret =  search_db (cursor, qry, verbose=False)
	if ret and ret[0][0]: new_column_value += ret[0][0]

	qry = "alter table {}.{} set auto_increment={}".format(db_name, table_name, new_column_value)
	search_db (cursor, qry, verbose=False)

	qry = "alter table {}.{} modify column {} int not null auto_increment".format(db_name, table_name, column)
	rows = search_db (cursor, qry, verbose=False)
	if (rows):
		print(rows)
		exit()


########
def check_null (variable):

	if variable is None:
		return None
	if (type(variable) is str and variable=="None"):
		return None
	return variable


########
def switch_to_db (cursor, db_name):
	qry = "use %s" % db_name
	rows = search_db (cursor, qry, verbose=False)
	if (rows):
		print(rows)
		return False
	return True


########
def get_table_size(cursor, db_name, tables):
	table_size = {}
	for table in tables:
		qry = "select count(1) from %s.%s" % (db_name,table)
		ret = search_db (cursor, qry, verbose=True)
		if ret and ret[0] and type(ret[0][0])==int:
			table_size[table] = ret[0][0]
		else:
			table_size[table] = 0

	return table_size


########
def val2mysqlval(value):
	if  value is None:
		return  "null "
	elif type(value) is str:
		return "\'%s\'" % value
	return "{}".format(value)


########
def store_without_checking(cursor, table, fields, verbose=False, database=None, ignore=False):

	qry = "insert "
	# if we use INSERT IGNORE, the duplication attempt is ignored
	if ignore: qry += "ignore "
	qry += "%s.%s "%(database,table) if database else "%s "%table

	qry += "("
	qry += ",".join(fields.keys())
	qry += ")"

	qry += " values "
	qry += "("
	qry += ",".join([val2mysqlval(v) for v in fields.values()])
	qry += ")"

	rows  = search_db (cursor, qry, verbose)
	if verbose: print("qry:",qry,"\n", "rows:", rows)

	if rows:
		rows   = search_db (cursor, qry, verbose=True)
		print(rows)
		return -1

	rows = search_db (cursor, "select last_insert_id()" )
	try:
		row_id = int(rows[0][0])
	except:
		row_id = -1
	return row_id


########
def store_or_update (cursor, table, fixed_fields, update_fields, verbose=False, primary_key='id'):

	conditions = " and ".join(["{}={}".format(k,val2mysqlval(v)) for k,v in fixed_fields.items()])

	# check if the row exists
	qry = "select %s from %s  where %s "  % (primary_key, table, conditions)
	rows   = search_db (cursor, qry, verbose)
	exists = rows and (type(rows[0][0]) is int)

	row_id = -1
	if exists: row_id = rows[0][0]
	if verbose: print("\n".join(["", qry, "exists? {}".format(exists), str(row_id)]))
	if exists and not update_fields: return row_id

	if exists: # if it exists, update
		if verbose: print("exists; updating")
		qry  = "update %s set " % table
		qry += ",".join(["{}={}".format(k,val2mysqlval(v)) for k,v in update_fields.items()])
		qry += " where %s " % conditions

	else: # if not, make a new one
		if verbose: print("does not exist; making new one")
		qry  = "insert into %s " % table
		keys = list(fixed_fields.keys())
		vals = list(fixed_fields.values())
		if update_fields:
			keys += list(update_fields.keys())
			vals += list(update_fields.values())
		qry += "(" + ",".join(keys) + ")"
		qry += " values "
		qry += "(" + ",".join([val2mysqlval(v) for v in vals]) + ")"

	rows   = search_db (cursor, qry, verbose)

	if verbose: print("qry:",qry,"\n", "rows:", rows)
	# if there is a return, it is an error msg
	if rows:
		rows   = search_db (cursor, qry, verbose=True)
		print(rows[0])
		return -1

	if row_id==-1:
		rows = search_db (cursor, "select last_insert_id()" )
		try:
			row_id = int(rows[0][0])
		except:
			row_id = -1
	return row_id


#########################################
def create_index (cursor, db_name, index_name, table, columns, verbose=False):

	if  not switch_to_db (cursor, db_name):
		return False

	# check whether this index exists already
	qry = "show index from %s where key_name like '%s'" % ( table, index_name)
	rows = error_intolerant_search(cursor, qry)
	if (rows):return True

	# columns is a list of columns that we want to have indexed
	qry = "create index %s on %s (%s)" % (index_name, table, ",".join(columns))
	rows = error_intolerant_search(cursor, qry)
	if (rows): return False
	return True


#########################################
def get_column_names (cursor, db_name, table_name):

	qry  = "select c.column_name from information_schema.columns c "
	qry += "where c.table_schema='%s' and c.table_name='%s'" % (db_name, table_name)

	rows = search_db (cursor, qry, verbose=False)
	if (rows):
		if ( 'Error' in rows[0]):
			rows = search_db (cursor, qry, verbose=True)
			return False
		else:
			return [row[0] for row in rows]
	else:
		return False


#########################################
def column_exists (cursor, db_name, table_name, column_name):

	if  not switch_to_db (cursor, db_name):
		return False

	qry = "show columns from "+ table_name + " like '%s'" % column_name
	rows = search_db (cursor, qry, verbose=False)
	if (rows):
		if ( 'Error' in rows[0]):
			return False
		else:
			return True
	else:
		return False

def add_column(cursor, db_name, table_name, column_name, col_type, default=None, after_col=None):
	if not column_exists (cursor, db_name, table_name, column_name):
		qry = "alter table  %s.%s add  %s %s  " %(db_name, table_name, column_name, col_type)
		if default: qry += "default %s " % default
		if after_col: qry += "after %s" % after_col
		error_intolerant_search(cursor,qry)
	return

#########################################
def add_boolean_column(cursor, db_name, table_name, column_name):
	return add_column(cursor, db_name, table_name, column_name, 'boolean', '0')

#########################################
def add_float_column(cursor, db_name, table_name, column_name):
	return add_column(cursor, db_name, table_name, column_name, 'float', '0.0')


#########################################
def entry_exists(cursor, db_name, table_name, column_name, column_value):
	qry = "select 1 from %s.%s where %s=%s" % (db_name, table_name, column_name, column_value)
	rows = search_db (cursor, qry, verbose=False)
	if rows:
		return True
	else:
		return False


#########################################
def check_table_exists (cursor, db_name, table_name):

	if  not switch_to_db (cursor, db_name):
		return False

	qry = "show tables like '%s'" % table_name
	rows = search_db (cursor, qry, verbose=False)
	if (rows):
		if ( 'Error' in rows[0]):
			return False
		else:
			return True
	else:
		return False

############
def check_and_drop(cursor, db_name, table):
	search_db(cursor, "drop table if exists %s.%s"% (db_name, table))
	return

#########################################
def table_create_time (cursor, db_name, table_name):

	qry  = "select create_time from information_schema.tables where "
	qry += "table_schema   = '%s' " % db_name
	qry += "and table_name = '%s' " % table_name

	rows = search_db (cursor, qry, verbose=False)
	if (not rows or  'Error' in rows[0]):
		search_db (cursor, qry, verbose=True)
		return ""
	else:
		return rows[0][0]

import warnings
#######
def search_db(cursor, qry, verbose=False):
	warnings.filterwarnings('ignore', category=MySQLdb.Warning)
	try:
		cursor.execute(qry)
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.execute() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Error"], e.args]
	except MySQLdb.Warning as e: # this does not work for me - therefore filterwarnings
		if verbose:
			print("Warning running cursor.execute() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Warning"], e.args]

	try:
		rows = cursor.fetchall()
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.fetchall() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Error"], e.args]

	if len(rows) == 0:
		if verbose:
			print("No return for query:\n%s" % qry)
		return False

	# since python3 fetchall returns bytes inst of str in some  random fashion
	# not clear what's going on
	# here is a rather useless issue page on github
	# https://github.com/PyMySQL/mysqlclient-python/issues/145#issuecomment-283936456
	rows_clean = []
	for row in rows:
		rows_clean.append([r.decode('utf-8') if type(r)==bytes else r for r in row])
	return rows_clean

#########################################
def error_intolerant_search(cursor, qry):
	ret =  search_db(cursor, qry)
	if not ret: return ret
	if type(ret[0][0])==str and 'error' in ret[0][0].lower():
		search_db(cursor, qry, verbose=True)
		exit()
	return ret

#########################################
def hard_landing_search(cursor, qry):
	ret =  search_db(cursor, qry)
	if not ret or (type(ret[0][0])==str and 'error' in ret[0][0].lower()):
		search_db(cursor, qry, verbose=True)
		exit()
	return ret



########
def connect_to_mysql (conf_file):
	try:
		mysql_conn_handle = MySQLdb.connect(read_default_file=conf_file)
	except  MySQLdb.Error as e:
		print("Error connecting to mysql (%s) " % (e.args[1]))
		sys.exit(1)
	return mysql_conn_handle


########
def connect_to_db (db_name, user=None, passwd=None):

	try:
		if not user is None:
			db = MySQLdb.connect(user=user, passwd=passwd, db=db_name)
		else:
			db = MySQLdb.connect(user="root", db=db_name)
	except  MySQLdb.Error as e:
		print("Error connecting to %s: %d %s" % (db_name, e.args[0], e.args[1]))
		exit(1)

	return db

