#!/usr/bin/env python
import os
import sys
print(sys.path)
import glob
from collections import namedtuple
import argparse
import imp




from IPython import embed
from sqlalchemy import create_engine
import sqlalchemy.exc as exc
import sqlalchemy.orm.exc as oexc

from obspy import read


import pisces as ps

from pisces.util import db_connect
import pisces.schema.css3 as kba
import pisces.tables.css3 as kb
import pisces.io.mseed as ms



# user supplies their own class, inherited from kbcore, or just uses .tables
# the prototype tables have a .from_sac or .from_mseed classmethod.

# for readability, use these named tuples for core tables, like:
# tab = CORETABLES[7]
# tab.name is 'site', tab.prototype is the abstract Site class,
# and tab.table is an actual Site table
CoreTable = namedtuple('CoreTable', ['name', 'prototype', 'table'])
CORETABLES = [CoreTable('affiliation', kba.Affiliation, kb.Affiliation),
			  CoreTable('arrival', kba.Arrival, kb.Arrival),
			  CoreTable('assoc', kba.Assoc, kb.Assoc),
			  CoreTable('event', kba.Event, kb.Event),
			  CoreTable('instrument', kba.Instrument, kb.Instrument),
			  CoreTable('lastid', kba.Lastid, kb.Lastid),
			  CoreTable('origin', kba.Origin, kb.Origin),
			  CoreTable('site', kba.Site, kb.Site),
			  CoreTable('sitechan', kba.Sitechan, kb.Sitechan),
			  CoreTable('wfdisc', kba.Wfdisc, kb.Wfdisc)]

# HELPER FUNCTIONS
def expand_glob(option, opt_str, value, parser):
	"""Returns an iglob iterator for file iteration. Good for large file lists."""
	setattr(parser.values, option.dest, glob.iglob(value))

def get_parser():
	"""
	This is where the command-line options are defined, to be parsed from argv

	Returns
	-------
	argparse.ArgumentParser instance

	Examples
	-------
	Test the parser with this syntax:

	>>> from sac2db import get_parser
	>>> parser = get_parser()
	>>> options = parser.parse_args(['--origin', 'origin', '--affiliation',
									 'myaffiliation', 'sqlite:///mydb.sqlite',
									 '*.sac'])
	>>> print options
	Namespace(affiliation='my.affiliation', arrival=None,
	assoc=None, event=None, files=['*.sac'], instrument=None, lastid=None,
	origin='origin', absolute_path=False, site=None, sitechan=None,
	url='sqlite:///mydb.sqlite', wfdisc=None)

	"""
	parser = argparse.ArgumentParser(prog='sac2db',
			formatter_class=argparse.RawDescriptionHelpFormatter,
			description="""
			Write data from SAC files into a database.

			Standard lowercase core table names are used by default, but can be
			prefixed (see --prefix). Individual table names can also be specified.

			Plugins
			-------
			Plugin functions can be added for custom processing.  Plugins must be
			functions with a signature like "myplugin(**kwargs)", which
			accept keyword arguments that include the canonical table names as
			keywords and lists of corresponding table instances as values,
			modifies them, and returns them. The --plugin flag can be specified
			multiple times, and each plugin is executed in the order in which
			it was specified.  Specified plugins are executed for each SAC
			file, after a collection of the core table instances are created.

			e.g.
			A plugin like this, which stuffs 'mysta' into all site.refsta
			attributes...

			# in plugins/sac.py
			def stuff_refsta(**kwargs):
				for site in kwargs.get('site', []):
					site.refsta = 'mysta'

				return kwargs

			...is called like this:
			sac2db.py --plugin plugins/sac:stuff_refsta sqlite:///mydb.sqlite *.sac

			More complex plugin function configuration may be done with the
			standard library ConfigParser module and text config files.


			Examples
			--------
			# use standard table names to local test.sqlite file
			sac2db.py sqlite:///test.sqlite datadir/*.sac

			# prefix all tables in an oracle account with my_, prompt for password
			sac2db.py oracle://user@server.domain.com:port/dbname --prefix my_ datadir/*.sac

			# if there are too many SAC files for the shell to handle, use a list:
			find datadir -name "*.sac" -print > saclist.txt
			sac2db.py sqlite:///test.sqlite saclist.txt

			""",
			version='0.2')
	# ----------------------- Add core table arguments ------------------------
	#The following loop adds the core table owner/name options.
	table_group = parser.add_argument_group('table name overriding',
			"Optionally, override individual standard table names using 'owner.tablename' naming. No owner for sqlite.")
	for coretable in CORETABLES:
		table_group.add_argument('--' + coretable.name,
							default=None,
							metavar='owner.tablename',
							dest=coretable.name)
	# -------------------------------------------------------------------------
	parser.add_argument('--plugin',
			default=None,
			metavar='path/to/module_file:plugin_function',
			dest='plugins',
			action='append',
			help="For each SAC file, import and execute an additional plugin \
				  function.  The function must accept keyword arguments that \
				  include canonical table name keywords and lists of table \
				  instances as values. ")

	parser.add_argument('--prefix',
			default='',
			metavar='owner.prefix',
			dest='prefix',
			help="Target tables using 'account.prefix naming.'\
				  e.g. myaccount.test_ will target tables \
				  like myaccount.test_origin, myaccount.test_sitechan.")

	parser.add_argument('--absolute_paths',
			default=False,
			help="Write database 'dir' directory entries as absolute paths, not relative.",
			action='store_true',
			dest='absolute_paths')



	parser.add_argument('url',
			help="SQLAlchemy-style database connection string, such as \
			sqlite:///mylocaldb.sqlite or oracle://myuser@myserver.lanl.gov:8000/mydb")
	
	parser.add_argument('siteTABLE',
		help="get site table")

	parser.add_argument('files',
			nargs='+',
			help="file names, including shell name expansions.")

	return parser


def get_plugins(options):
	"""Returns a list of imported plugin function objects."""
	plugin_functions = []
	if options.plugins:
		for plugin_string in options.plugins:
			try:
				pth, fn = plugin_string.split(':')
			except ValueError:
				msg = "Must specify plugin like: path/to/modulename:plugin_function"
				raise ValueError(msg)
			pth = pth.split(os.path.sep)
			modname = pth.pop(-1)
			f, pathname, descr = imp.find_module(modname, pth)
			mod = imp.load_module(modname, f, pathname, descr)
			plugin_functions.append(getattr(mod, fn))

	return plugin_functions


def get_session(options):
	# accept command line arguments, return a database-bound session.
	#session = url_connect(options.url)
	session =db_connect(options.url)
	return session


def get_files(options):
	"""
	Return a sequence of WFDISC file names from either a list of file names
	(trivial) or a text file list (presumable because there are too many files
	to use normal shell expansion).

	"""
	if len(options.files) == 1:
		#make a generator of non-blank lines
		try:
			listfile = open(options.files[0], 'r')
			files = (line.strip() for line in listfile if line.strip())
		except IOError:
			msg = "{0} does not exist.".format(options.files[0])
			raise IOError(msg)
	else:
		files = options.files
	return files

def get_wfdisc(options):
	"""
	Return a sequence of WFDISC file names from either a list of file names
	(trivial) or a text file list (presumable because there are too many files
	to use normal shell expansion).

	"""
	if len(options.files) == 1 and not isSAC(options.files[0]):
		#make a generator of non-blank lines
		try:
			listfile = open(options.files[0], 'r')
			files = (line.strip() for line in listfile if line.strip())
		except IOError:
			msg = "{0} does not exist.".format(options.files[0])
			raise IOError(msg)
	else:
		files = options.files

	return files


def get_or_create_tables(options, session, create=True):
	"""
	Load or create canonical ORM KB Core table classes.

	Parameters
	----------
	options : argparse.ArgumentParser
	session : sqlalchemy.orm.Session

	Returns
	-------
	tables : dict
		Mapping between canonical table names and SQLA ORM classes.
		e.g. {'origin': MyOrigin, ...}

	"""
	# The Plan:
	# 1. For each core table, build or get the table name
	# 2. If it's a vanilla table name, just use a pre-packaged table class
	# 3. If not, try to autoload it.
	# 4. If it doesn't exist, make it from a prototype and create it in the database.

	# TODO: check options for which tables to produce.

	dbout = options.prefix

	tables = {}
	for coretable in CORETABLES:
		# build the table name
		if getattr(options, coretable.name, None):
			fulltablename = getattr(options, coretable.name)
		else:
			fulltablename = dbout + coretable.name

		# fulltablename is either an arbitrary string or dbout + core name, but not None

		# put table classes into the tables dictionary
		if fulltablename == coretable.name:
			# it's a vanilla table name. just use a pre-packaged table class instead of making one.
			tables[coretable.name] = coretable.table
		else:
			tables[coretable.name] = ps.make_table(fulltablename, coretable.prototype)

		tables[coretable.name].__table__.create(session.bind, checkfirst=True)

	session.commit()

	return tables


def dicts2rows(dicts, classes):
	for table, dcts in list(dicts.items()):
		cls = classes[table]
		dicts[table] = [cls(**dct) for dct in dcts]

	return dicts


def make_atomic(last, **rows):
	"""
	Unify related table instances/row, including: ids, dir, and dfile
	"""
	# last is an AttributeDict of {'keyvalue': lastid instance, ...}
	# rows is a dictionary of {'canonical tablename': [list of instances], ...}
	# of _related_ instances from a single SAC header?
	# TODO: check existance of rows before changing their ids.

	#print rows
	# the order matters here

	# for SAC, only 1
	for event in rows.get('event', []):
		# skips if no 'event' key and if 'event' value is empty []
		# XXX: check for existence first
		event.evid = next(last.evid)

	# for SAC, only 1
	for origin in rows.get('origin', []):
		# XXX: check for existance first
		origin.orid = next(last.orid)
		origin.evid = event.evid

	# for SAC, only 1
	for affil in rows.get('affiliation', []):
		pass

	# for SAC, only 1
	for sitechan in rows.get('sitechan', []):
		# XXX: check for existance first
		sitechan.chanid = next(last.chanid)

	# arrivals correspond to assocs
	for (arrival, assoc) in zip(rows.get('arrival', []), rows.get('assoc', [])):
		arrival.arid = next(last.arid)
		arrival.chanid = sitechan.chanid

		assoc.arid = arrival.arid
		assoc.orid = origin.orid

	# for SAC, only 1
	for wfdisc in rows.get('wfdisc', []):
		wfdisc.wfid = next(last.wfid)


def apply_plugins(plugins, **rows):
	for plugin in plugins:
		rows = plugin(**rows)

	return rows


def ms2db(files, url, tables, plugins=None, abs_paths=False):
	pass

# TODO: make this main also accept a get_iterable and get_row_dicts functions,
#   so it can be renamed to iter2db and re-used in a sac2db.py and miniseed2db.py
import time



def main(argv=None):
	"""
	Command-line arguments are created and parsed, fed to functions.

	"""
	parser = get_parser()
	
	options = parser.parse_args(argv)

	# options_dict = vars(options)

	session = get_session(options)

	files = get_files(options)
	tables = get_or_create_tables(options, session, create=True)
	#lastids = ['arid', 'chanid', 'evid', 'orid', 'wfid']
	#last = get_lastids(session, tables['lastid'], lastids, create=True)

	# for item in iterable:
	
	last_time = time.time()
	countI=0
	staTOT=[]
	

	if not (options.siteTABLE==False):
		with open(options.siteTABLE) as si:
			for line in si:
				isi=tables['site'].from_string(line,default_on_error=['lddate'])
				session.add(isi)
		session.commit()
	#continue
	
	
	k=0
	embed()
	for wfile in files:		
		with open(wfile) as wf:
			for line in wf:
				iwf=tables['wfdisc'].from_string(line)
				iwf.dir='2012/'+iwf.dir
				iwf.wfid=k
				k=k+1
				session.add(iwf)
		session.commit()
	

		
		
		
		'''		
		#print 'last:', last_time-time.time()
		countI=countI+1
		start_time = time.time()
		print msfile

		# row_dicts = get_row_dicts(item)

		tr = read(msfile, format='MSEED', headonly=True)[0]
		#print time.time()-start_time
		# rows needs to be a dict of lists, for make_atomic
		# row_dicts = get_row_dicts(tr.stats.sac) # put in the whole trace, to determine byte order?
		# sac
		

		dicts = ms.mshdr2tables(tr.stats, tables=tables.keys())
		#print time.time()-start_time
		# row_instances = dicts_to_instances(row_dicts, tables)

		rows = dicts2rows(dicts, tables)
		#print time.time()-start_time
		
		# manage dir, dfile, datatype
		# XXX: hack.  replace with an updated obspy.io.sac later.

		bo = tr.data.dtype.byteorder
		if bo == '<':
			datatype = 'f4'
		elif bo == '>':
			datatype = 't4'
		elif bo == '=':
			if sys.byteorder == 'little':
				datatype = 'f4'
			else:
				datatype = 't4'

		for wf in rows['wfdisc']:
			wf.datatype = 'ms'
			wf.dfile = os.path.basename(msfile)
			if options.absolute_paths:
				idir = os.path.dirname(os.path.realpath(msfile))
			else:
				idir = os.path.dirname(msfile)
			wf.dir = idir

		# manage the ids
		#print time.time()-start_time
		make_atomic(last, **rows)

		#plugins = get_plugins(options)

		#rows = apply_plugins(plugins, **rows)

		# add rows to the database
		# XXX: not done very elegantly.  some problem rows are simply skipped.
		#countI=0
		#staTOT=[]
		for table, instances in rows.items():

			#print table
			if not table=='site':
				continue
			
			if instances[0].sta in staTOT:
				continue
				
			   
			if instances:
				
				# could be empty []
				try:
					
					session.add_all(instances)
					session.commit()
					staTOT.append(instances[0].sta)
					
						
				except exc.IntegrityError as e:
					# duplicate or nonexistant primary keys
					session.rollback()
					print("rollback {}".format(table))
				except exc.OperationalError as e:
					# no such table, or database is locked
					session.rollback()
					print("rollback {}".format(table))
					
					
		#countI=0
		for table, instances in rows.items():
			#print table
			if not table=='wfdisc':
				continue
				
			
			if instances:
				try:
					session.add_all(instances)
					#countI=countI+1
					if countI>100:
						print 'commiting'
						session.commit()
						countI=0
						
				except exc.IntegrityError as e:
					# duplicate or nonexistant primary keys
					session.rollback()
					print("rollback {}".format(table))
				except exc.OperationalError as e:
					# no such table, or database is locked
					session.rollback()
					print("rollback {}".format(table))
					

		#print time.time()-start_time
		'''

		#last_time = time.time()


if __name__ == '__main__':
	main(sys.argv[1:])
