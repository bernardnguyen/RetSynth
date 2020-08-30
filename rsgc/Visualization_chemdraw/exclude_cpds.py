import os
import re
PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = os.sep.join(PATH.split(os.sep)[:-1]) + os.sep + 'Database' + os.sep + 'data'

def retrieve_promiscuous_mets(DB):
	promiscuous = {}
	with open(PPATH+os.sep+'promiscuous_cpds_inchi.txt') as f:
		line = f.readline()[:-1]
		while line:
			cpd_name = DB.get_compound_name(line)
			
			if cpd_name == 'None':
				print('\nINFO:\tPromiscuous compound %s not found' % line)
			else:
				promiscuous[cpd_name] = True
			
			line = f.readline()[:-1]
	return promiscuous