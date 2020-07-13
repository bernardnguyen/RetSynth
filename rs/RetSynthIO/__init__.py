'''
RetSynth I/O Interface
'''
import os,sys,re

import datetime
from Database import query as Q

import rs
# BP_FILE = sys.path[0] + '/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC.sp'

def generate_unique_ID():
	now = datetime.datetime.now()
	return 'RS_%02d%02d%d_%02d%02d%02d' % (now.month, now.day, now.year, now.hour, now.minute, now.second)

def get_new_output_folder(PATH, count):
	'''Check if folder to store images is already present if so new temp folder is made'''
	count+=1
	try:
		os.mkdir(PATH+str(count))
		return PATH+str(count)
	except OSError:
		PATH_NEW = get_new_output_folder(PATH, count)
		return PATH_NEW
		
def generate_output_directory(run_id, directory=''):
	HOMEDIR = os.path.expanduser('~')
	
	if directory:
		output_dir = get_new_output_folder(directory+'/output_retsynth', 0)
		return output_dir
	else:
		if os.path.isdir(HOMEDIR + '/Desktop'):
			try:
				os.mkdir(HOMEDIR + '/Desktop/RetSynth')			
			except:
				pass
		else:
			try:
				os.mkdir(HOMEDIR + '/RetSynth')				
			except:
				pass

		if os.path.isdir(HOMEDIR + '/Desktop/RetSynth'):			
			OUTPUT_DIR = HOMEDIR + '/Desktop/RetSynth/'
		else:
			OUTPUT_DIR = HOMEDIR + '/RetSynth/'
			
		os.mkdir(OUTPUT_DIR + run_id)
		return OUTPUT_DIR + run_id + '/'

def generate_target_file(run_id, start, target, directory):
	f = open(directory + '/targetfile.txt','w')	
	f.write('#compoundid\torganismid\n')
	f.write('%s\t%s' % (target,start))
	f.close()

def parameterize_retsynth(run_id, request_form):
	sys.argv = ['rs.py', '-t', request_form.get("output_path") + '/targetfile.txt']
	sys.argv += ['-op', request_form.get("output_path")]
	# sys.argv += ['-rank_bp','-rank_bp_file',BP_FILE]
	sys.argv += ['--figures']
	sys.argv += ['-v']
	print('')	

	# DATABASE PARAMETERS
	if request_form.get("database_file"):
		sys.argv += ['-db',request_form.get("database_file")]
		print("PARAMETERS:\tUsing database %s" % request_form.get("database_file").split('/')[-1])

		corresponding_constraints = request_form.get("database_file")[:-3]
		# EVALRXNS
		if request_form.get("evalrxns") != 'all':
			sys.argv += ['-evalrxns', request_form.get("evalrxns")]
			print('PARAMETERS:\tOnly evaluate %s reactions' % request_form.get("evalrxns"))

			corresponding_constraints += "_" + request_form.get("evalrxns")
		corresponding_constraints += ".constraints"

		if os.path.isfile(corresponding_constraints):
			sys.argv += ['-dbc',corresponding_constraints]
		else:
			sys.argv += ['-gdbc',corresponding_constraints]		
			print("PARAMETERS:\t\tNo corresponding constraints file found - generating a new one")	
	elif request_form.get("gdbname"):
		OUTPUT_DIR = '/'.join(request_form.get("output_path").split('/')[:-2])
		sys.argv += ['-gdb','%s/%s.db' % (OUTPUT_DIR,request_form.get("gdbname"))]
		sys.argv += ['-gdbc','%s/%s.constraints' % (OUTPUT_DIR,request_form.get("gdbname"))]
		print("PARAMETERS:\tGenerating database and constraints - %s" % request_form.get("gdbname"))

		# InChI DB
		if int(request_form.get("inchi_check")) != 0:
			sys.argv += ['--inchidb']
			print('PARAMETERS:\tUsing InChI IDs')

		# PATRIC DB
		if int(request_form.get("patric_check")) != 0:
			sys.argv += ['--patric_models']
			print('PARAMETERS:\tUsing PATRIC data')

			if request_form.get("patric_un") != '' and request_form.get("patric_pw") != '':
				sys.argv += ['-p_un',request_form.get("patric_un"),'-p_pw',request_form.get("patric_pw")]

			if request_form.get("patric_media") != 'Complete':
				sys.argv += ['-patricrxntype',request_form.get("patric_media")]
				print('PARAMETERS:\t\tPATRIC media - %s' % request_form.get("patric_media"))

			if request_form.get("patric_file"):
				sys.argv += ['--patricfile',request_form.get("patric_file")]
				print('PARAMETERS:\t\tPATRIC file - %s' % request_form.get("patric_file"))

		# METACYC DB
		if int(request_form.get("metacyc_check")) != 0:
			sys.argv += ['--metacyc']
			print('PARAMETERS:\tUsing MetaCyc data')

			if request_form.get("metacyc_addition"):
				sys.argv += ['--metacyc_addition',request_form.get("metacyc_addition")]
				print('PARAMETERS:\t\tUsing MetaCyc XML file - %s' % request_form.get("metacyc_addition"))


		# KEGG DB
		if int(request_form.get("kegg_check")) != 0:
			sys.argv += ['--kegg']
			print('PARAMETERS:\tUsing KEGG data')

			if request_form.get("kegg_org_type") != 'all':
				sys.argv += ['-keggorganismtype',request_form.get("kegg_org_type")]
				print('PARAMETERS:\t\tInclude KEGG organism of type - %s' % request_form.get("kegg_org_type"))

			if int(request_form.get("kegg_nunorgs")) != 0:
				sys.argv += ['-keggnunorganisms', str(request_form.get("kegg_nunorgs"))]
				print('PARAMETERS:\t\t%d Organisms' % request_form.get("kegg_nunorgs"))

			if int(request_form.get("kegg_nunorg_paths")) != 0:
				sys.argv += ['-keggnunorganismpaths',str(request_form.get("kegg_nunorg_paths"))]
				print('PARAMETERS:\t\t%d Organism paths' % request_form.get("kegg_nunorg_paths"))

	# K PARAMETER
	if int(request_form.get("kparam")) != 0:
		sys.argv += ['-k', str(request_form.get("kparam"))]
		print('PARAMETERS:\tK Parameter: %s' % request_form.get("kparam"))

	# XLSX OUTPUT
	if request_form.get("oxf_check"):
		sys.argv += ['-oxf']
		print('PARAMETERS:\tOutput as .xlsx')

	# TANIMOTO THRESHOLD ANALYSIS
	if request_form.get("tanimoto_check"):
		tan_thresh = request_form.get("tanimoto_thresh")
		sys.argv += ['-run_tan_thresh']	
		sys.argv += ['-tan_thresh', str(tan_thresh)]
		print('PARAMETERS:\tTanimoto Structural Similarity Analysis: Threshold = %.2f' % tan_thresh)
	elif request_form.get("tanimoto_thresh") != 0.9:
		tan_thresh = request_form.get("tanimoto_thresh")
		sys.argv += ['-tan_thresh', str(tan_thresh)]
		print('PARAMETERS:\tTanimoto Threshold = %.2f' % tan_thresh)

	# FBA
	if request_form.get("fba_check"):
		sys.argv += ['-fba']
		print('PARAMETERS:\tFlux Balance Analysis')
	
	# WITH KNOCKOUTS
	if request_form.get("ko_check"):
		sys.argv += ['-ko']
		print('PARAMETERS:\t\tWith Knockouts')

	# GENE COMPATIBILITY
	if request_form.get("gene_compatibility_check"):
		sys.argv += ['-gc']
		print('PARAMETERS:\tOptimize gene compatibility')
		
		caithresh = request_form.get("cai_optimal_threshold")
		if caithresh != 0.5 :
			sys.argv += ['-cai_thresh', str(caithresh)]
			print('PARAMETERS:\t\tCodon Adaptation Index (CAI) Optimal Threshold = %.2f' % caithresh)


	# IUPAC NAME
	if request_form.get("iupac_name_check"):
		sys.argv += ['--use_iupac_names']
		print('PARAMETERS:\tUsing IUPAC names for reaction solvents and catalysts')

	# CHEMICAL RXN INFO
	if request_form.get("chemical_rxn_info_check"):
		sys.argv += ['--show_rxn_info']
		print('PARAMETERS:\tDisplay Chemical Reaction Information')


def find_matches(search_term, database):
	DB = Q.Connector(database)
	return DB.get_all_cpd_with_search(search_term)


def run_retsynth(request_form, stdout_q=None):
	run_id = generate_unique_ID()
	startingOrganism = request_form.get("start").split(' -- ')[0]
	targetCompound = request_form.get("target").split(' -- ')[1]

	if request_form.get("output_path"):
		directory = request_form.get("output_path")
		output_directory = generate_output_directory(run_id,directory=directory)
		generate_target_file(run_id, startingOrganism, targetCompound, output_directory)
		request_form['output_path']=output_directory
	else:
		output_directory = generate_output_directory(run_id)		
		generate_target_file(run_id, startingOrganism, targetCompound, output_directory)
		request_form["output_path"] = output_directory

	if stdout_q:
		sys.stdout = stdout_q
		logfile = request_form['output_path'] + '/logfile.txt'
		stdout_q.setLogFile(logfile)
		sys.stderr = open(logfile,'a')

	print('STATUS:\tRunning RetSynth (ID: %s)...\n' % run_id)
	print('PARAMETERS:\tStart From Model:\t%s' % startingOrganism)
	print('PARAMETERS:\tTarget Compound:\t%s' % targetCompound)

	parameterize_retsynth(run_id, request_form)
	
	try:
		rs.main()
		print('\nSTATUS:\tRetSynth Complete')
		print('INFO:\tResults available at %s\n' % str(request_form["output_path"]))
	except:
		print('\nSTATUS:\tRetSynth has failed...')
		print('INFO:\tError details may be available at %s' % logfile)
		pass

	return run_id