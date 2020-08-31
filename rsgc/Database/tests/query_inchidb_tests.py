from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that query the database'

import sys
import glob
import os
import re
import sqlite3
import unittest
import cobra
from Database import initialize_database as init_db
from Database import build_kbase_db as bkdb
from Database import query as Q

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)

sbml_files = glob.glob(os.path.join(PATH+'/data', '*'))
all_mets = []
all_reactions = []
all_rev = {}

'''Import model information directly from sbml files'''
print ('STATUS OF TESTS: LOADING TEST MODELS FROM SBML FILES')
for file_name in sbml_files:
    model = cobra.io.read_sbml_model(file_name)
    for r in model.reactions:
        i = r.id
        revers = r.reversibility
        all_reactions.append(i)
        if i not in all_rev:
            all_rev[i] = []
            all_rev[i].append(revers)
        elif i in all_rev and revers not in all_rev[i]:
            all_rev[i].append(revers)
    for m in model.metabolites:
        i = m.id
        all_mets.append(i)
all_mets = list(set(all_mets))
all_reactions = list(set(all_reactions))
for k, v in all_rev.iteritems():
    if len(v) > 2:
        print (' ERROR IN TEST CODE')
        sys.exit()
print ('STATUS OF TESTS: FINISHED LOADING TEST MODELS FROM SBML FILES')
init_db.Createdb(PATH+'/testinchi.db', True)
bkdb.BuildKbase(PATH+'/data', PPATH+'/data/KbasetoKEGGCPD.txt', PPATH+'/data/KbasetoKEGGRXN.txt',
                True, PATH+'/testinchi.db', 'bio')
DB = Q.Connector(PATH+'/testinchi.db')

class Generate_databaseTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_get_all_models(self):
        print ("Testing the get_all_models() function")
        hits = DB.get_all_models()
        self.assertEqual(len(hits), 3)

    def test_get_compound_name(self):
        print ("Testing the get_compound_name() function")
        hit = DB.get_compound_name('InChI=1S/C6H10O6/c7-2-5(10)3(8)1-4(9)6(11)12/h3,5,7-8,10H,1-2H2,(H,11,12)/t3-,5+/m0/s1_c0')
        self.assertEqual('2_keto_3_deoxygluconate_c0', hit)

    def test_get_reactions(self):
        print ("Testing the get_reactions() function")
        hits = DB.get_reactions('InChI=1S/C6H10O6/c7-2-5(10)3(8)1-4(9)6(11)12/h3,5,7-8,10H,1-2H2,(H,11,12)/t3-,5+/m0/s1_c0', 1)
        self.assertIn('rxn13_c0', hits)
        hits = DB.get_reactions('InChI=1S/C6H10O6/c7-2-5(10)3(8)1-4(9)6(11)12/h3,5,7-8,10H,1-2H2,(H,11,12)/t3-,5+/m0/s1_c0', 0)
        self.assertIn('rxn7_c0', hits)

    def test_get_reaction_name(self):
        print ("Testing the get_reaction_name() function")
        hit = DB.get_reaction_name("rxn7_c0")
        self.assertEqual('rxn7_c0', hit)

    def test_get_reaction_species(self):
        print ("Testing the get_reaction_species() function")
        hit = DB.get_reaction_species('rxn7_c0')
        self.assertIn('t2', hit)

    def test_get_reactants(self):
        print ("Testing the get_reactants() function")
        hits = DB.get_reactants('rxn6_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('cpdE_c0', hits)

        hits = DB.get_reactants('rxn3_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('cpdB_c0', hits)

    def test_get_products(self):
        print ("Testing the get_products() function")
        hits = DB.get_products('rxn6_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', hits)

        hits = DB.get_products('rxn3_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('InChI=1S/2ClH.2H2N.Pt/h2*1H;2*1H2;/q;;2*-1;+4/p-2_c0', hits)

    def test_get_compounds_in_model(self):
        print ("Testing the get_compounds_in_model() function")
        hits = DB.get_compounds_in_model('t2')
        model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
        self.assertEqual(len(hits), len(model.metabolites))
        self.assertIn("cpdB_c0", hits)
        self.assertIn("InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0", hits)

    def test_get_reactions_in_model(self):
        print ("Testing the get_reactions_in_model() function")
        hits = DB.get_reactions_in_model('t2')
        model = cobra.io.read_sbml_model(PATH+'/data/test2.xml')
        self.assertEqual(len(hits), len(model.reactions))
        self.assertIn("rxn7_c0", hits)

    def test_is_reversible(self):
        print ("Testing the is_reversible() function")
        hits = DB.is_reversible('t1', "rxn3_c0")
        self.assertEqual(hits, 'true')

        hits = DB.is_reversible('t1', "rxn2_c0")
        self.assertEqual(hits, 'false')

    def test_get_all_reactions(self):
        print ("Testing the get_all_reactions() function")
        hits = DB.get_all_reactions()
        self.assertEqual(len(all_reactions), len(hits))

    def test_get_stoichiometry(self):
        print ("Testing the get_stoichiometry() function")
        hit = DB.get_stoichiometry('rxn11_c0', 'cpdI_c0', 0)
        self.assertEqual(hit[0], 1)
        hit = DB.get_stoichiometry('rxn11_c0', 'cpdH_c0', 1)
        self.assertEqual(hit[0], 1)

        hit = DB.get_stoichiometry('rxn1_c0', 'InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0', 0)
        self.assertEqual(hit[0], 1)
        hit = DB.get_stoichiometry('rxn1_c0', 'cpdT_c0', 1)
        self.assertEqual(hit[0], 1)

    def test_is_reversible_all(self):
        print ("Testing the is_reversible_all() function")
        if r in all_rev:
            if len(all_rev[r]) == 2:
                hit = DB.is_reversible_all(r)
                self.assertEqual(hit, 'true')
            else:
                hit = DB.is_reversible_all(r)
                self.assertEqual(hit, all_rev[r])
        self.assertEqual(DB.is_reversible_all('rxn2_c0'), 'true')

    def test_get_compound_compartment(self):
        print ("Testing the get_compound_compartment() function")
        hit = DB.get_compound_compartment('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
        self.assertEqual(hit, 'c0')
        hit = DB.get_compound_compartment('cpdI_e0')
        self.assertEqual(hit, 'e0')

    def test_get_genes(self):
        print ("Testing the get_genes() function")
        hits = DB.get_genes('rxn11_c0', 't1')
        self.assertEqual('(g11)', hits)

        hits = DB.get_genes('rxn12_c0', 't1')
        self.assertEqual('(g2)', hits)

        hits = DB.get_genes('EX_I_e0', 't1')
        self.assertEqual(hits, 'None')

    def test_get_proteins(self):
        print ("Testing the get_proteins() function")
        hits = DB.get_proteins('rxn11_c0', 't1')
        self.assertEqual('(p11)', hits)

        hits = DB.get_proteins('rxn12_c0', 't1')
        self.assertEqual('None', hits)

        hits = DB.get_proteins('EX_I_e0', 't1')
        self.assertEqual(hits, 'None')

    def test_get_organism_name(self):
        print ("Testing the get_organism_name() function")
        hit = DB.get_organism_name('t3')
        self.assertEqual(hit, 'test3.xml')

    def test_get_uniq_metabolic_clusters(self):
        print ("Testing the get_uniq_metabolic_clusters_function()")
        hit = DB.get_uniq_metabolic_clusters()
        self.assertEqual(len(hit), 2)

    def test_get_models_from_cluster(self):
        print ("Testing the get_models_from_cluster() function")
        cnx = sqlite3.connect(PATH+'/testinchi.db')
        Q = cnx.execute('SELECT DISTINCT cluster_num  FROM cluster')
        hits = list(set(Q.fetchall()))
        self.assertEqual(len(hits), 2)

    def test_get_reactions_based_on_type(self):
        hits = DB.get_reactions_based_on_type('bio')
        self.assertEqual(len(hits), 12)

    def test_get_compartment(self):
        hits = DB.get_compartment('cytosol')
        self.assertIn('c0', hits)

    # def test_get_all_cpd_fp(self):
    #     print ("Testing the get_cpd_chemicalformula() function")
    #     hit = DB.get_all_cpd_fp()
    #     cpd_hit = DB.get_all_compounds()
    #     self.assertEqual(len(hit), len(cpd_hit))

    # def test_get_cpd_fp(self):
    #     print ("Testing the get_cpd_fp() function")
    #     hit = DB.get_cpd_fp('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
    #     self.assertEqual(hit, '15,12,0,32,-128,8,3,91,65,2,34,18,10,-75,2,-42,26,-124,114,-127,26,18,2,65,-54,-48,2,4,-128,16,17,-126,33,-20,4,32,16,10,-20,1,-104,0,77,0,1,-95,67,65,-112,88,7,-8,-16,-80,-128,117,40,64,4,-128,80,44,16,8,80,36,-8,42,0,1,-62,-90,36,66,-30,24,-31,56,-120,-83,34,32,64,-62,34,8,35,4,27,-117,24,2,34,77,81,-126,-92,9,65,1,-112,65,6,0,4,21,16,20,-92,28,-120,113,2,2,84,99,12,0,-124,-60,32,-29,80,116,88,32,-124,-62,80,84,81,-56,5,32,-30,-20,-74,64,48,-96,4,25,76,32,-124,8,-108,34,115,22,-128,48,88,96,37,72,-127,113,-54,1,32,12,37,72,6,36,2,3,54,65,35,64,-103,-84,68,-128,22,8,-112,-128,56,80,10,2,64,8,-96,-124,74,18,-120,0,64,-126,-112,21,49,70,32,-31,24,67,45,0,64,-128,27,5,64,10,0,0,-124,32,-103,48,0,-16,32,-58,0,72,-95,60,2,98,-122,0,0,16,70,0,1,32,0,0,68,80,96,-124,74,41,0,-54,1,72,42,4,2,4,80,-128,-60,0,5,-64,36,65,8,0,4,0,6,0,-104,-128,72,-43,-31,102,-54,-49,60,11,36,9,-87,-77,-62,60,2,62,-115,8,-24,14,117,40,122,-63,-31,-50,32,-78,90,35,-28,54,105,35,117,18,9,96,66,-112,-47,-127,-121,-87,-34,35,20,-9,-24,91,-124,-13,-42,-68,2,-108,-2,-122,4,-13,108,57,-111,-37,-115,-42,8,113,-64,-125,-15,34,122,-23,-40,-112,12,-70,-32,-81,105,-124,84,89,22,-57,98,30,-88,-78,3,-120,23,10,48,-46,25,-23,-24,-47,120,24,124,-79,-122,62,64,64,120,-80,39,-96,65,119,-128,22,-24,46,6,88,3,41,114,-90,-126,12,-40,26,95,-126,-112,-35,31,80,19,32,-25,66,8,109,-52,71,-86,-55,33,-36,-25,5,82,65,-35,12,88,-27,2,-115,-59,-21,-96,-8,-55,2,64,-114,111,43,-118,13,-122,-124,27,-47,1,78,115,-98,74,115,11,19,83,20,67,19,22,-96,-116,81,89,50,75,32,64,8,11,2,33,96,-84,0,54')

    def test_get_cpd_cf(self):
        print ("Testing the get_cpd_chemicalformula() function")
        hit = DB.get_cpd_chemicalformula('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
        self.assertEqual(hit, 'C20H15ClF2N2O2')    

    def test_get_all_cpd_chemicalformulas(self):
        print ("Testing the get_cpd_chemicalformula() function")
        hit = DB.get_all_cpd_chemicalformulas()
        cpd_hit = DB.get_all_compounds()
        self.assertEqual(len(hit), len(cpd_hit))

    def test_get_all_cpd_with_chemicalformula(self):
        print ("Testing the get_cpd_chemicalformula() function")
        hit = DB.get_all_cpd_with_chemicalformula('None')
        self.assertEqual(len(hit), 8) 

    def test_get_cpd_casnumber(self):
        print ("Testing the get_cpd_casnumber() function")
        hit = DB.get_cpd_casnumber('InChI=1S/C20H15ClF2N2O2/c21-14-6-4-13(5-7-14)12-27-15-8-9-24-19(10-15)25-20(26)11-16-17(22)2-1-3-18(16)23/h1-10H,11-12H2,(H,24,25,26)_c0')
        self.assertEqual(hit, 'None') 

if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/testinchi.db')
