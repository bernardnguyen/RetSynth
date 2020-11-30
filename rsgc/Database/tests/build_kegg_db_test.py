from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'leanne382@gmail.com'
__description__ = 'Runs tests on codes that add reactions and compounds to database from kegg'

import os
import re
import sqlite3
import unittest
from copy import deepcopy
from shutil import copyfile
from rsgc.Database import query as Q
from rsgc.Database import build_KEGG_db as bkeggdb
PATH = os.path.dirname(os.path.abspath(__file__))
copyfile(PATH+'/datam/testPATRIC.db', PATH+'/datam/testPATRIC_kegg.db')
copyfile(PATH+'/datam/testPATRICinchi.db', PATH+'/datam/testPATRICinchi_kegg.db')

DB = Q.Connector(PATH+'/datam/testPATRIC_kegg.db')
rxnsadd1 = deepcopy(DB.get_all_reactions())

DBinchi = Q.Connector(PATH+'/datam/testPATRICinchi_kegg.db')
rxns1inchi = deepcopy(DBinchi.get_all_reactions())

class BuildKEGG(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_KEGG_with_only_KEGG_orgs(self):
        """test KEGG database construction"""
        K = bkeggdb.CompileKEGGIntoDB(PATH+'/datam/testPATRIC_kegg.db', 'bacteria', False, 1, 2, 2, 'bio', True)
        DBadd = Q.Connector(PATH+'/datam/testPATRIC_kegg.db')
        newrxns = set()
        with open('newrxns.txt') as fin:
            for line in fin:
                line = line.strip()
                newrxns.add(line)

        inter = set(newrxns) - set(rxnsadd1)

        kegg_rxnsadd = DBadd.get_all_reactions()
        self.assertEqual(len(rxnsadd1)+len(inter), len(kegg_rxnsadd)-2) # minus 2 biomass reactions

        reactants = DBadd.get_reactants('R02187_c0')
        products = DBadd.get_products('R02187_c0')
        self.assertIn('C00404_c0', reactants)
        self.assertIn('C00221_c0', reactants)
        self.assertIn('C00404_c0', products)
        self.assertIn('C01172_c0', products)

        cf = DBadd.get_cpd_chemicalformula('C00404_c0')
        self.assertEqual(cf, None)
        cf = DBadd.get_cpd_chemicalformula('C00221_c0')
        self.assertEqual(cf, None)

        cas = DBadd.get_cpd_casnumber('C00404_c0')
        self.assertEqual(cas, None)
        cas = DBadd.get_cpd_casnumber('C00221_c0')
        self.assertEqual(cas, None)       

        reactants = DBadd.get_reactants('R07618_c0')
        products = DBadd.get_products('R07618_c0')
        self.assertIn('C00579_c0', reactants) #different than below because Kbase kegg translation fie 
        self.assertIn('C00003_c0', reactants)
        self.assertIn('C15972_c0', products)
        self.assertIn('C00004_c0', products)
        self.assertIn('C00080_c0', products)

        allcpds = DBadd.get_all_compounds()
        self.assertEqual(len(allcpds), len(set(allcpds)))

    def test_KEGG_with_only_KEGG_orgsinchi(self):
        """test KEGG database construction"""
        K = bkeggdb.CompileKEGGIntoDB(PATH+'/datam/testPATRICinchi_kegg.db', 'bacteria', True, 4, 2, 2, 'bio', True)
        DBaddinchi = Q.Connector(PATH+'/datam/testPATRICinchi_kegg.db')
        newrxns = set()
        with open('newrxns.txt') as fin:
            for line in fin:
                line = line.strip()
                newrxns.add(line)

        inter = set(newrxns) - set(rxns1inchi)
        kbaserxnsadd = DBaddinchi.get_all_reactions()
        self.assertEqual(len(rxns1inchi)+len(inter), len(kbaserxnsadd)-2) # minus 2 biomass reactions

        reactants = DBaddinchi.get_reactants('R02187_c0')
        products = DBaddinchi.get_products('R02187_c0')
        self.assertIn('C00404_c0', reactants)
        self.assertIn('C00221_c0', reactants)
        self.assertIn('C00404_c0', products)
        self.assertIn('C01172_c0', products)
        inchis = DBaddinchi.get_all_compounds_inchi()
        self.assertIn('InChI=1S/H5O10P3/c1-11(2,3)9-13(7,8)10-12(4,5)6/h(H,7,8)(H2,1,2,3)(H2,4,5,6)', inchis)
        self.assertIn('InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1', inchis)
        self.assertIn('InChI=1S/H5O10P3/c1-11(2,3)9-13(7,8)10-12(4,5)6/h(H,7,8)(H2,1,2,3)(H2,4,5,6)', inchis)
        self.assertIn('InChI=1S/C6H13O9P/c7-3-2(1-14-16(11,12)13)15-6(10)5(9)4(3)8/h2-10H,1H2,(H2,11,12,13)/t2-,3-,4+,5-,6-/m1/s1', inchis)
        self.assertIn('InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p+1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1', inchis)
        self.assertIn('InChI=1S/C21H29N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,25)/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1', inchis)

        cas = DBaddinchi.get_cpd_casnumber('C00404_c0')
        self.assertEqual(cas, None) 
        cas = DBaddinchi.get_cpd_casnumber('C00311_c0')
        self.assertEqual(cas, '320-77-4')
   
        cf = DBaddinchi.get_cpd_chemicalformula('C00404_c0')
        self.assertEqual(cf, 'H5O10P3')
        cf = DBaddinchi.get_cpd_chemicalformula('C01172_c0')
        self.assertEqual(cf, 'C6H13O9P')
        
        cf = DBaddinchi.get_cpd_chemicalformula('C15973_c0')
        self.assertEqual(cf, None)

        reactants = DBaddinchi.get_reactants('R07618_c0')
        products = DBaddinchi.get_products('R07618_c0')
        self.assertIn('C00579_c0', reactants)
        self.assertIn('C00003_c0', reactants)
        self.assertIn('C15972_c0', products)
        self.assertIn('C00080_c0', products)
        self.assertIn('C00004_c0', products)

        allcpds = DBaddinchi.get_all_compounds()
        self.assertEqual(len(allcpds), len(set(allcpds)))

if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/datam/testPATRIC_kegg.db')
    os.remove(PATH+'/datam/testPATRICinchi_kegg.db')
