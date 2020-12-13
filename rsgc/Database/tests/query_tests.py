from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'leanne382@gmail.com'
__description__ = 'Runs tests on codes that query the database'

import sys
import glob
import os
import re
import unittest
import cobra
import sqlite3
from shutil import copyfile
from rsgc.Database import query as Q
from sys import platform
if platform == "cygwin":
    header_path = 'C:\cygwin64'
else:
    header_path = ''
PATH = os.path.dirname(os.path.abspath(__file__))
copyfile(PATH+'/datam/testPATRICinchi.db', PATH+'/datam/testPATRICinchi_query.db')

'''GENERATE TEST DATABASE'''
sbml_files = glob.glob(os.path.join(PATH+'/testssbml_models/', '*'))
all_mets = []
all_reactions = []
all_rev = {}
'''Import model information directly from sbml files'''
print ('STATUS OF TESTS: LOADING TEST MODELS FROM SBML FILES')
for file_name in sbml_files:
    model = cobra.io.read_sbml_model(header_path+file_name)
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
for k, v in all_rev.items():
    if len(v) > 2:
        print (' ERROR IN TEST CODE')
        sys.exit()
print ('STATUS OF TESTS: FINISHED LOADING TEST MODELS FROM SBML FILES')
DB = Q.Connector(PATH+'/datam/testPATRICinchi_query.db')


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
        self.assertEqual(len(hits), 2)

    def test_get_compound_name(self):
        print ("Testing the get_compound_name() function")
        hit = DB.get_compound_name('C00001_c0')
        self.assertEqual('H2O', hit)

    def test_get_reactions(self):
        print ("Testing the get_reactions() function")
        hits = DB.get_reactions('C14787_c0', 0)
        self.assertEqual(len(hits), 2)
        self.assertIn('R07003_c0', hits)

    def test_get_reaction_name(self):
        print ("Testing the get_reaction_name() function")
        hit = DB.get_reaction_name("R07004_c0")
        self.assertEqual("(1S)-hydroxy-(2S)-glutathionyl-1,2-dihydronaphthalene glutathione-lyase (epoxide-forming)", hit)

    def test_get_reaction_species(self):
        print ("Testing the get_reaction_species() function")
        hit = DB.get_reaction_species('R07003_c0')
        self.assertIn(str(329726.14), hit)

    def test_get_reactants(self):
        print ("Testing the get_reactants() function")
        hits = DB.get_reactants('R07003_c0')
        self.assertEqual(len(hits), 2)
        self.assertIn('C14787_c0', hits)

    def test_get_products(self):
        print ("Testing the get_products() function")
        hits = DB.get_products('R07003_c0')
        self.assertEqual(len(hits), 1)
        self.assertIn('C14793_c0', hits)

    def test_get_compounds_in_model(self):
        print ("Testing the get_compounds_in_model() function")
        hits = DB.get_compounds_in_model(str(329726.14))
        model = cobra.io.read_sbml_model(header_path+PATH+'/testssbml_models/Acaryochloris_marina_MBIC11017_Complete_Complete')
        self.assertEqual(len(hits), len(model.metabolites))

    def test_get_reactions_in_model(self):
        print ("Testing the get_reactions_in_model() function")
        hits = DB.get_reactions_in_model(str(329726.14))
        model = cobra.io.read_sbml_model(header_path+PATH+'/testssbml_models/Acaryochloris_marina_MBIC11017_Complete_Complete')
        self.assertEqual(len(hits), len(model.reactions))

    def test_is_reversible(self):
        print ("Testing the is_reversible() function")
        hits = DB.is_reversible(str(551115.6), "R05592_c0")
        self.assertEqual(hits, '1')

        hits = DB.is_reversible(str(329726.14), 'R07003_c0')
        self.assertEqual(hits, '1')

    def test_get_all_compounds(self):
        print ("Testing the get_all_compounds() function")
        hits = DB.get_all_compounds()
        self.assertEqual(len(hits), len(all_mets))

    def test_get_all_reactions(self):
        print ("Testing the get_all_reactions() function")
        hits = DB.get_all_reactions()
        self.assertEqual(len(all_reactions), len(hits)-1)

    def test_get_stoichiometry(self):
        print ("Testing the get_stoichiometry() function")
        hit = DB.get_stoichiometry('R08191_c0', 'C00080_c0', 1)
        self.assertEqual(hit[0], 2)
        hit = DB.get_stoichiometry('R08191_c0', 'C16170_c0', 0)
        self.assertEqual(hit[0], 1)

    def test_is_reversible_all(self):
        print ("Testing the is_reversible_all() function")
        if r in all_rev:
            if len(all_rev[r]) == 2:
                hit = DB.is_reversible_all(r)
                self.assertEqual(hit, True)
            else:
                hit = DB.is_reversible_all(r)
                self.assertEqual(hit, all_rev[r])

    def test_get_compound_compartment(self):
        print ("Testing the get_compound_compartment() function")
        hit = DB.get_compound_compartment('C16170_c0')
        self.assertEqual(hit, 'c0')
        hit = DB.get_compound_compartment('C00080_e0')
        self.assertEqual(hit, 'e0')

    def test_get_genes(self):
        print ("Testing the get_genes() function")
        hits = DB.get_genes('R02568_c0', str(329726.14))
        self.assertEqual('329726.14.peg.1346', hits)

        hits = DB.get_genes('R02568_c0', str(551115.6))
        self.assertEqual('551115.6.peg.6020', hits)

    def test_get_proteins(self):
        print ("Testing the get_proteins() function")
        hits = DB.get_proteins('R02568_c0', str(329726.14))
        self.assertEqual('(EC 4.1.2.13)', hits)


    def test_get_organism_name(self):
        print ("Testing the get_organism_name() function")
        hit = DB.get_organism_name(str(329726.14))
        self.assertEqual(hit, 'Acaryochloris_marina_MBIC11017_Complete')

    def test_get_uniq_metabolic_clusters(self):
        print ("Testing the get_uniq_metabolic_clusters_function()")
        hit = DB.get_uniq_metabolic_clusters()
        self.assertEqual(len(hit), 2)

    def test_get_reactions_based_on_type(self):
        hits = DB.get_reactions_based_on_type('bio')
        hits1 = DB.get_all_reactions()
        self.assertEqual(len(hits), len(hits1))

    def test_get_compartment(self):
        print ("Testing getting compartment ID from name")
        hits = DB.get_compartment('cytosol')
        self.assertIn('c0', hits)

    def test_get_cpd_cf(self):
        print ("Testing the get_cpd_chemicalformula() function")
        hit = DB.get_cpd_chemicalformula('C00080_c0')
        self.assertEqual(hit, 'H')
        hit = DB.get_cpd_chemicalformula('C00437_c0')
        self.assertEqual(hit, 'C7H14N2O3')

    def test_get_cpd_casnumber(self):
        print ("Testing the get_cpd_casnumber() function")
        hit = DB.get_cpd_casnumber('C00437_c0')
        self.assertEqual(hit, '6205-08-9') 

    def test_get_cpd_inchi(self):
        print ("Testing the inchi functions() function")
        hit = DB.get_inchi_from_compoundID('C00437_c0')
        self.assertEqual(hit, 'InChI=1S/C7H14N2O3/c1-5(10)9-6(7(11)12)3-2-4-8/h6H,2-4,8H2,1H3,(H,9,10)(H,11,12)/t6-/m0/s1') 
        hit = DB.get_compound_ID_from_inchi('InChI=1S/C7H14N2O3/c1-5(10)9-6(7(11)12)3-2-4-8/h6H,2-4,8H2,1H3,(H,9,10)(H,11,12)/t6-/m0/s1')
        self.assertEqual(hit, 'C00437_c0') 

if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/datam/testPATRICinchi_query.db')


