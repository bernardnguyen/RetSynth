"""
RetroSynth GUI
"""
import os, sys, tarfile


# if application is frozen, need local executable path
if getattr(sys, 'frozen', False):
	application_path = os.path.dirname(sys.executable)
	os.environ["PATH"] += str(os.pathsep + application_path)
	DEFAULT_DB = application_path + '/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC_cas_SPRESI_reduced1x.db'
else:
	DEFAULT_DB = sys.path[0] + '/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC_cas_SPRESI_reduced1x.db'


from rs_gui.retsynthwindows import MainWindow
# FIRST TIME RUNNING RETSYNTH AFTER WRAPPING
if not os.path.isfile(DEFAULT_DB):
	# Make temporary directories for gene compatibility module
	from GeneCompatibility.gc import make_tmp_folders
	make_tmp_folders()
	
	# Unzip database
	DB_DIR = '/'.join(DEFAULT_DB.split('/')[:-1])
	with tarfile.open(DEFAULT_DB[:-2] + 'tar.gz','r:gz') as tar:
		tar.extractall(DB_DIR)


from Database import query as Q
DB = Q.Connector(DEFAULT_DB)
models = ['%s -- %s' % (m[0],m[1]) for m in DB.get_all_fba_models()]



def main():
	MainWindow(DB,DEFAULT_DB,models)


if __name__ == '__main__':
	main()