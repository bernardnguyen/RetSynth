from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'leanne382@gmail.com'
__description__ = 'Runs tests on codes add reactions and compounds to existing database'

import os
import re
import sqlite3
import unittest
from shutil import copyfile
from rsgc.Database import query as Q
from rsgc.Database.build_metacyc_db import MetaCyc
from rsgc.Database import build_metacyc_db as bmcdb
path2metacyc = input("Please input path to full metacyc ... ")

PATH = os.path.dirname(os.path.abspath(__file__))
copyfile(PATH+'/datam/testPATRIC.db', PATH+'/datam/testPATRIC_mc.db')


DB = Q.Connector(PATH+'/datam/testPATRIC_mc.db')
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
file_name.close()

class Translate_metacycTests(unittest.TestCase):
    def test_metacyc(self):
        """Tests that metacyc can be added to a patric database correctly"""
        print ("Testing that metacyc can be added to a patric database correctly")
        MC = MetaCyc(DB, False, sqlite3.connect(PATH+'/datam/testPATRIC_mc.db'), False)
        MC.read_metacyc_file(BIOCYC_translator, path2metacyc)
        self.assertIn(('R04139_c0', 'META', '(G-44401)'), MC.genelist)
        self.assertIn(('R04139_c0', 'META', True), MC.modelreactions)
        self.assertIn(('R00850_c0', 'META', True), MC.modelreactions)
        self.assertIn(('R00850_c0_v1', 'META', True), MC.modelreactions)
        self.assertIn(('R00850_c0_v2', 'META', True), MC.modelreactions)
        self.assertIn(('R00850_c0_v3', 'META', True), MC.modelreactions)
        self.assertEqual(MC.all_rxnreversibility['R04139_c0'], True)
        self.assertEqual(MC.all_rxnreversibility['R00850_c0'], True)
        self.assertEqual(MC.all_rxnreversibility['R00850_c0_v1'], True)
        self.assertEqual(MC.all_rxnreversibility['R00850_c0_v2'], True)
        self.assertEqual(MC.all_rxnreversibility['R00850_c0_v3'], True)
        self.assertEqual(5, len(MC.all_reaction_compound['R04139_c0']))

        self.assertIn(('R04139_c0', 'C00003_c0', False, '1'),
                      MC.all_reaction_compound['R04139_c0'])
        self.assertIn(('R04139_c0', 'C00004_c0', True, '1'),
                      MC.all_reaction_compound['R04139_c0'])
        self.assertIn(('R04139_c0', 'C00003_c0', False, '1'),
                      MC.all_reaction_compound['R04139_c0'])
        # self.assertIn(('TRANS__45__RXN0__45__443_t0', 'C00638_p0', False, '1'),
        #               MC.all_reaction_compound['TRANS__45__RXN0__45__443_t0'])
        # self.assertIn(('TRANS__45__RXN0__45__443_t0', 'C00638_c0', True, '1'),
        #               MC.all_reaction_compound['TRANS__45__RXN0__45__443_t0'])

        self.assertIn(('R06683_c0', 'C12424_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00005_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00080_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00007_c0', False, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C12425_c0', True, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C00006_c0', True, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'C12425_c0', True, '1'),
                      MC.all_reaction_compound['R06683_c0'])
        self.assertIn(('R06683_c0', 'META', '(G-17686)'),
                      MC.genelist)
        self.assertIn(('R01173_c0', 'META', '(CA_P0035) or (CA_P0162)'), MC.genelist)
        self.assertIn(('R06683_c0', 'META', '1.14.13.180'),
                      MC.proteinlist)
        self.assertIn(('R00850_c0', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn(('R00850_c0_v1', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn(('R00850_c0_v2', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn(('R00850_c0_v3', 'META', '2.7.1.142'), MC.proteinlist)
        self.assertIn('C01456_c0', MC.all_cpds)
        self.assertEqual(['tropate', 'c0', 'C01456', None, '529-64-6', None], MC.all_cpds['C01456_c0'])

    def test_BIOCYC_translator(self):
        """Test that translator dictionary is built correctly"""
        print ("Testing that translator dictionary is built correctly")
        T = bmcdb.Translate(PATH+'/datam/testPATRIC_mc.db',path2metacyc,
                            False, 'bio', False, add=False)
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
    os.remove(PATH+'/datam/testPATRIC_mc.db')
