from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Runs tests on codes that add reactions and compounds to database from MINE'

import os
import re
import sqlite3
import unittest
from copy import deepcopy
from rsgc.Database import query as Q
from rsgc.Database import initialize_database as init_db
from rsgc.Database import build_kbase_db as bkdb
from rsgc.Database import build_MINE_db as bminedb

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)

init_db.Createdb(PATH+'/kbasetestaddMINE.db', False)
bkdb.BuildKbase(PATH+'/datam', PPATH+'/data/KbasetoKEGGCPD.txt', PPATH+'/data/KbasetoKEGGRXN.txt',
                False, PATH+'/kbasetestaddMINE.db', 'bio')

init_db.Createdb(PATH+'/kbasetestaddinchiMINE.db', True)
bkdb.BuildKbase(PATH+'/data', PPATH+'/data/KbasetoKEGGCPD.txt', PPATH+'/data/KbasetoKEGGRXN.txt',
                True, PATH+'/kbasetestaddinchiMINE.db', 'bio')

class BuildMINEtests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_MINE_no_inchi(self):
        ''' test addition of reaction info from MINE raw files'''
        bminedb.BuildMINEdb(PATH+'/data5', PATH+'/kbasetestaddMINE.db', False, 'bio')
        DB = Q.Connector(PATH+'/kbasetestaddMINE.db')
        cnx = sqlite3.connect(PATH+'/kbasetestaddMINE.db')
        QC = cnx.execute("""SELECT * FROM model WHERE ID = ?""", ("MINE",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('MINE', results)

        QC = cnx.execute("""SELECT * FROM cluster WHERE ID = ?""", ("MINE",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('18', results)

        MINEreactions = DB.get_reactions_in_model('MINE')
        self.assertEqual(len(MINEreactions), 6)
       	self.assertIn('rxn1_m_c0', MINEreactions)

        MINEcpds = DB.get_compounds_in_model('MINE')
        self.assertIn('C9ff38ba752ec99fdf7db17911973b1447ed37e2e_c0', MINEcpds)
        self.assertEqual(len(MINEcpds), 12)

        allcpds = DB.get_all_compounds()
        self.assertIn('C9ff38ba752ec99fdf7db17911973b1447ed37e2e_c0', allcpds)

        allrxns = DB.get_all_reactions()
        self.assertIn('rxn1_m_c0', allrxns)

        cf = DB.get_cpd_chemicalformula('C9ff38ba752ec99fdf7db17911973b1447ed37e2e_c0')
        self.assertEqual(cf, 'None')

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ?""", ("C9ff38ba752ec99fdf7db17911973b1447ed37e2e_c0",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

    def test_MINE_inchi(self):
        ''' test addition of reaction info from MINE raw files'''
        bminedb.BuildMINEdb(PATH+'/data5', PATH+'/kbasetestaddinchiMINE.db', True, 'bio')
        DBinchi = Q.Connector(PATH+'/kbasetestaddinchiMINE.db')
        cnx = sqlite3.connect(PATH+'/kbasetestaddinchiMINE.db')

        QC = cnx.execute("""SELECT * FROM model WHERE ID = ?""", ("MINE",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('MINE', results)

        QC = cnx.execute("""SELECT * FROM cluster WHERE ID = ?""", ("MINE",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('3', results)

        MINEreactions = DBinchi.get_reactions_in_model('MINE')
        self.assertEqual(len(MINEreactions), 6)
       	self.assertIn('rxn1_m_c0', MINEreactions)

        MINEcpds = DBinchi.get_compounds_in_model('MINE')
        self.assertIn('InChI=1S/C5H9NO6S/c7-5(8)4-1-3(2-6-4)12-13(9,10)11/h3-4,6H,1-2H2,(H,7,8)(H,9,10,11)_c0', MINEcpds)
        self.assertEqual(len(MINEcpds), 12)

        allcpds = DBinchi.get_all_compounds()
        self.assertIn('InChI=1S/C5H9NO6S/c7-5(8)4-1-3(2-6-4)12-13(9,10)11/h3-4,6H,1-2H2,(H,7,8)(H,9,10,11)_c0', allcpds)

        cf = DBinchi.get_cpd_chemicalformula('InChI=1S/C5H9NO6S/c7-5(8)4-1-3(2-6-4)12-13(9,10)11/h3-4,6H,1-2H2,(H,7,8)(H,9,10,11)_c0')
        self.assertEqual(cf, 'C5H9NO6S')

        allrxns = DBinchi.get_all_reactions()
        self.assertIn('rxn1_m_c0', allrxns)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ?""", ("InChI=1S/C5H9NO6S/c7-5(8)4-1-3(2-6-4)12-13(9,10)11/h3-4,6H,1-2H2,(H,7,8)(H,9,10,11)_c0",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)


if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/kbasetestaddMINE.db')
    os.remove(PATH+'/kbasetestaddinchiMINE.db')
