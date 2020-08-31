from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes add reactions and compounds to existing database'

import os
import re
import sqlite3
import unittest
from rsgc.Database import query as Q
from rsgc.Database import initialize_database as init_db
from rsgc.Database.build_metacyc_db import MetaCyc
from rsgc.Database import build_metacyc_db as bmcdb
from rsgc.Database import build_kbase_db as bkdb

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)

if os.path.isfile(PATH+'/kbasetest.db') is True:
    os.remove(PATH+'/kbasetest.db')
init_db.Createdb(PATH+'/kbasetest.db', True)
bkdb.BuildKbase(PATH+'/data', PPATH+'/data/KbasetoKEGGCPD.txt', PPATH+'/data/KbasetoKEGGRXN.txt',
                True, PATH+'/kbasetest.db', 'bio')
DB = Q.Connector(PATH+'/kbasetest.db')
file_name = open(PATH+'/metacyc_data/MetaCyc.aliases')
line = file_name.readline()
BIOCYC_translator = {}
BIOCYC_translator['rxn'] = {}
BIOCYC_translator['compound'] = {}
if line.startswith('#'):
    pass
for count, line in enumerate(file_name):
    larray = line.strip('\n').split('\t')
    larray[0] = re.sub('\.\w+$', '', larray[0])
    if larray[1] != '':
        if larray[0] not in BIOCYC_translator['rxn'].keys():
            BIOCYC_translator['rxn'][larray[0]] = []
            BIOCYC_translator['rxn'][larray[0]].append(larray[1])
        else:
            BIOCYC_translator['rxn'][larray[0]].append(larray[1])
    elif larray[2] != '':
        if larray[0] not in BIOCYC_translator['compound']:
            BIOCYC_translator['compound'][larray[0]] = []
            BIOCYC_translator['compound'][larray[0]].append(larray[2])
        else:
            BIOCYC_translator['compound'][larray[0]].append(larray[2])

class Translate_metacycTests(unittest.TestCase):
    def test_metacyc(self):
        """Tests that metacyc can be added to a kbase database correctly"""
        print ("Testing that metacyc can be added to a kbase database correctly")
        MC = MetaCyc(DB, True, sqlite3.connect(PATH+'/kbasetest.db'), True)
        MC.read_metacyc_file(BIOCYC_translator, PATH+'/metacyc_data/metabolic-reactions_inchi_test.xml')
        self.assertNotIn(('R03563_c0', 'META', '()'), MC.genelist)
        self.assertIn(('R03563_c0', 'META', 'false'), MC.modelreactions)
        self.assertEqual(MC.all_rxnreversibility['R03563_c0'], 'false')
        self.assertEqual(5, len(MC.all_reaction_compound['R03563_c0']))

        self.assertIn(('R03563_c0', 'InChI=1S/H2O/h1H2_c0', False, '1'),
                      MC.all_reaction_compound['R03563_c0'])
        self.assertIn(('R03563_c0', 'InChI=1S/C17H23NO3/c1-18-13-7-8-14(18)10-15(9-13)21-17(20)16(11-19)12-5-3-2-4-6-12/h2-6,13-16,19H,7-11H2,1H3/p+1/t13-,14+,15+,16?_c0', False, '1'),
                      MC.all_reaction_compound['R03563_c0'])
        self.assertIn(('R03563_c0', 'InChI=1S/p+1_c0', True, '1'),
                      MC.all_reaction_compound['R03563_c0'])
        self.assertIn(('R03563_c0', 'InChI=1S/C9H10O3/c10-6-8(9(11)12)7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/p-1_c0', True, '1'),
                      MC.all_reaction_compound['R03563_c0'])
        self.assertIn(('R03563_c0', 'InChI=1S/C8H15NO/c1-9-6-2-3-7(9)5-8(10)4-6/h6-8,10H,2-5H2,1H3/p+1/t6-,7+,8+_c0', True, '1'),
                      MC.all_reaction_compound['R03563_c0'])

        self.assertIn(('R03563_c0', 'META', '3.1.1.10'), MC.proteinlist)
        self.assertIn('InChI=1S/C9H10O3/c10-6-8(9(11)12)7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/p-1_c0', MC.all_cpds)
        print (MC.all_cpds['InChI=1S/C9H10O3/c10-6-8(9(11)12)7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/p-1_c0'])
        self.assertEqual(['tropate', 'c0', 'C01456', 'C9H9O3', '529-64-6'],MC.all_cpds['InChI=1S/C9H10O3/c10-6-8(9(11)12)7-4-2-1-3-5-7/h1-5,8,10H,6H2,(H,11,12)/p-1_c0'])

    def test_BIOCYC_translator(self):
        """Test that translator dictionary is built correctly"""
        print ("Testing that translator dictionary is built correctly")
        T = bmcdb.Translate(PATH+'/kbasetest.db', PATH+'/metacyc_data/metabolic-reactions_inchi_test.xml',
                            True, 'bio', False, add=False)
        self.assertEqual(len(T.BIOCYC_translator['rxn']), len(BIOCYC_translator['rxn']), 'bio')
        self.assertEqual(len(T.BIOCYC_translator['compound']), len(BIOCYC_translator['compound']))
        self.assertEqual(len(BIOCYC_translator['compound'])+len(BIOCYC_translator['rxn']), 26233)
        self.assertEqual(T.BIOCYC_translator['compound']['1-RADYL-2-ACYL-SN-GLYCERO-3-PHOSPHOLIPID'],
                         ['cpd21786'])
        self.assertIn('rxn16863', T.BIOCYC_translator['rxn']['1.1.1.190-RXN'])
        self.assertIn('rxn01938', T.BIOCYC_translator['rxn']['1.1.1.190-RXN'])
        self.assertEqual(len(T.BIOCYC_translator['rxn']['1.1.1.190-RXN']), 2)

if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/kbasetest.db')
