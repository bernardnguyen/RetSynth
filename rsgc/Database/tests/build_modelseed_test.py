from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that add reactions and compounds to database from model seed/patric database'

import os
import re
import sqlite3
import unittest
from copy import deepcopy
from rsgc.Database import query as Q
from rsgc.Database import initialize_database as init_db
from rsgc.Database import build_modelseed as bms
from rsgc.Database import mackinac
username = raw_input("Please input patric username...")
password = raw_input("Please input patric password...")

PATH = os.path.dirname(os.path.abspath(__file__))

init_db.Createdb(PATH+'/testPATRIC.db', False)

# init_db.Createdb(PATH+'/testinchiPATRIC.db', False)
token = mackinac.get_token(username=username, password=password)
models_in_ws = mackinac.list_patric_models()
models_in_ws_id = [m['id'] for m in models_in_ws]

print ('Get 551115.6')
if 'testmodel1' not in models_in_ws_id:
    patric_model1 = mackinac.create_patric_model(str(551115.6), 'testmodel1', media_reference='/{}@patricbrc.org/home/media/Complete'.format(username)) #Carbon-D-Glucose'
print ('Get 329726.14')
if 'testmodel2' not in models_in_ws_id:
    patric_model2 = mackinac.create_patric_model(str(329726.14), 'testmodel2', media_reference='/{}@patricbrc.org/home/media/Complete'.format(username)) #Carbon-D-Glucose'

print ('Get cobra model')
cobra_model1 = mackinac.create_cobra_model_from_patric_model('testmodel1')
cobra_model2 = mackinac.create_cobra_model_from_patric_model('testmodel2')

cobra_model1_met = len(cobra_model1.metabolites)
cobra_model1_rxn = len(cobra_model1.reactions)

cobra_model2_met = len(cobra_model2.metabolites)
cobra_model2_rxn = len(cobra_model2.reactions)

cobra1_set_mets = set()
cobra2_set_mets = set()
cobra1_set_rxns = set()
cobra2_set_rxns = set()
for m in cobra_model1.metabolites:
    cobra1_set_mets.add(m.id)
for m in cobra_model2.metabolites:
    cobra2_set_mets.add(m.id)
for m in cobra_model1.reactions:
    cobra1_set_rxns.add(m.id)
for m in cobra_model2.reactions:
    cobra2_set_rxns.add(m.id)

combined_mets = cobra1_set_mets.union(cobra2_set_mets)
combined_rxns = cobra1_set_rxns.union(cobra2_set_rxns)

class BuildKEGG(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")
    
    def test_patric(self):
        bms.BuildModelSeed(username=username, password=password, rxntype='bio', inchidb=False,
                           DBpath=PATH+'/testPATRIC.db', output_folder=PATH, newdb=True,  media='Complete',
                           sbml_output=False, processors=1, patricfile=PATH+'/data8/PATRIC_genome_complete_07152018.csv')
        DB = Q.Connector(PATH+'/testPATRIC.db')
        m1_cpds = DB.get_compounds_in_model(str(551115.6))
        m2_cpds = DB.get_compounds_in_model(str(329726.14))
        self.assertEqual(cobra_model1_met, len(m1_cpds))
        self.assertEqual(cobra_model2_met, len(m2_cpds))
        m1_rxns = DB.get_reactions_in_model(str(551115.6))
        m2_rxns = DB.get_reactions_in_model(str(329726.14))
        self.assertEqual(cobra_model1_rxn, len(m1_rxns))
        self.assertEqual(cobra_model2_rxn, len(m2_rxns))

        all_cpds = DB.get_all_compounds()
        all_rxns = DB.get_all_reactions()
        self.assertEqual(len(all_cpds), len(combined_mets))
        self.assertEqual(len(all_rxns)-1, len(combined_rxns)) # -1 to deal with biomass reactions
if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/testPATRIC.db')