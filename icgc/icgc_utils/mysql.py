import MySQLdb, sys, os

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
        print rows
        return False
    return True

########
def store_without_checking(cursor, table, fields, verbose=False):
    qry = "insert into %s " % table
    qry += "("
    first = True
    for field in fields.keys(): # again will have to check for the type here
        if (not first):
            qry += ", "
        qry += field
        first = False
    qry += ")"

    qry += " values "
    qry += "("
    first = True
    for value in fields.values(): # again will have to check for the type here
        if (not first):
            qry += ", "
        if  value is None:
            qry += " null "
        elif type(value) is int:
            qry += " %d" % value
        elif type(value) is float:
            qry += " %f" % value
        else:
            qry += " \'%s\'" % value
        first = False
    qry += ")"


    rows   = search_db (cursor, qry)
    if verbose:
        print
        print " ** ", qry
        print " ** ", rows

    if (rows):
        rows   = search_db (cursor, qry, verbose=True)
        print rows
        return False

    return True

########
def store_or_update (cursor, table, fixed_fields, update_fields, verbose=False):

    conditions = ""
    first = True
    for [field, value] in fixed_fields.iteritems():
        if (not first):
            conditions += " and "
        if ( type (value) is int):
            conditions += " %s= %d "  % (field, value)
        elif value is None:
            conditions += " %s= null" % field
        else:
            conditions += " %s='%s' " % (field, value)
        first = False

    # check if the row exists
    qry = "select exists (select 1 from %s  where %s) "  % (table, conditions)
    rows   = search_db (cursor, qry)
    exists = rows and (type(rows[0][0]) is long) and (rows[0][0]==1)

    if verbose:
        print
        print qry
        print rows
        print "exists ?", exists


    if exists and not update_fields: return True


    if exists: # if it exists, update
        if verbose: print "exists; update"
        qry = "update %s set " % table
        first = True
        for field, value in update_fields.iteritems():
            if (not first):
                qry += ", "
            qry += " %s = " % field
            if  value is None:
                qry += " null "
            elif type(value) is int:
                qry += " %d" % value
            else:
                qry += " \'%s\'" % value

            first = False
        qry += " where %s " % conditions

    else: # if not, make a new one
        if verbose: print "does not exist; make new one"

        qry = "insert into %s " % table
        qry += "("
        first = True
        for field in fixed_fields.keys()+update_fields.keys(): # again will have to check for the type here
            if (not first):
                qry += ", "
            qry += field
            first = False
        qry += ")"
       
        qry += " values "
        qry += "("
        first = True
        for value in fixed_fields.values()+update_fields.values(): # again will have to check for the type here
            if (not first):
                qry += ", "
            if  value is None:
                qry += " null "
            elif type(value) is int:
                qry += " %d" % value
            elif type(value) is float:
                qry += " %f" % value
            else:
                # strip its own quotes
                value =value.replace("'", "").replace('"', '')
                qry += " \"%s\"" % value
            first = False
        qry += ")"
        
    rows   = search_db (cursor, qry)

    if verbose:
        print
        print " ** ", qry
        print " ** ", rows

    if (rows):
        rows   = search_db (cursor, qry, verbose=True)
        return False
       
    
    return True

#########################################
def create_index (cursor, db_name, index_name, table, columns):

    if  not switch_to_db (cursor, db_name):
        return False
    
    # check whether this index exists already
    qry = "show index from %s where key_name like '%s'" % ( table, index_name) 
    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        print rows
        return True
   
    # columns is a list of columns that we want to have indexed
    qry = "create index %s  on %s " % (index_name, table)
    qry += " ("
    first = True
    for column in columns:
        if (not first):
            qry += ", "
        qry += column
        first = False
    qry += ")"


    rows = search_db (cursor, qry, verbose=True)
    if (rows):
        print rows
        return False
   
    return True

#########################################
def get_column_names (cursor, db_name, table_name):

    if  not switch_to_db (cursor, db_name):
        return False

    qry = "show columns from "+ table_name
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
def check_column_exists (cursor, db_name, table_name, column_name):
    
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

#######
def search_db (cursor, qry, verbose=False):

    try:
        cursor.execute(qry)
    except MySQLdb.Error, e:
        if verbose:
            print "Error running cursor.execute() for  qry: %s: %s " % (qry, e.args[1])
        return  [["ERROR: "+e.args[1]]]


    try:
        rows = cursor.fetchall()
    except MySQLdb.Error, e:
        if verbose:
            print "Error running cursor.fetchall() for  qry: %s: %s " % (qry, e.args[1])
        return  ["ERROR: "+e.args[1]]
    
    if (len(rows) == 0):
        if verbose:
            print "No return for query %s"  % qry
        return False

    return rows

########
def connect_to_mysql (user=None, passwd=None, host=None, port=None, conf_file=None):

    if conf_file:
        if not  os.path.isfile(conf_file):
            print conf_file, "not found or is not a file"
            return None
        try: # user is spelled out in the conf file
            mysql_conn_handle = MySQLdb.connect(read_default_file=conf_file)
        except  MySQLdb.Error, e:
            print "Error connecting to mysql (%s) " % (e.args[1])
            sys.exit(1) 
        return mysql_conn_handle
    
    else:
       
        try:
            if (user is None):
                mysql_conn_handle = MySQLdb.connect(user="root")
            elif (host is None):
                mysql_conn_handle = MySQLdb.connect(user=user, passwd=passwd)
            else:
                mysql_conn_handle = MySQLdb.connect(user=user, passwd=passwd, host=host, port=port)
            
        except  MySQLdb.Error, e:
            print "Error connecting to mysql as %s" % (e.args[1])
            sys.exit(1)

        return mysql_conn_handle

########
def connect_to_db (db_name, user=None, passwd=None):

    try:
        if not user is None:
            db = MySQLdb.connect(user=user, passwd=passwd, db=db_name)
        else:
            db = MySQLdb.connect(user="root", db=db_name)
    except  MySQLdb.Error, e:
        print "Error connecting to %s: %d %s" % (db_name, e.args[0], e.args[1])
        exit(1)

    return db

