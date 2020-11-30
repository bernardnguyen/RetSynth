from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'leanne382@gmail.com'
__description__ = 'Runs tests on codes that add reactions and compounds to database from MINE'

import os
import re
import sqlite3
import unittest
from copy import deepcopy
from shutil import copyfile
from rsgc.Database import query as Q
from rsgc.Database import build_ATLAS_db as batlasdb

PATH = os.path.dirname(os.path.abspath(__file__))
copyfile(PATH+'/datam/testPATRIC.db', PATH+'/datam/testPATRIC_atlas.db')
copyfile(PATH+'/datam/testPATRICinchi.db', PATH+'/datam/testPATRICinchi_atlas.db')

class BuildATLAStests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_ATLAS_no_inchi(self):
        ''' test addition of reaction info from MINE raw files'''
        batlasdb.build_atlas(PATH+'/data_atlas', PATH+'/datam/testPATRIC_atlas.db', False, 1, 'bio')
        DB = Q.Connector(PATH+'/datam/testPATRIC_atlas.db')
        cnx = sqlite3.connect(PATH+'/datam/testPATRIC_atlas.db')
        QC = cnx.execute("""SELECT * FROM model WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('ATLAS', results)

        QC = cnx.execute("""SELECT * FROM cluster WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('3', results)

        MINEreactions = DB.get_reactions_in_model('ATLAS')
        # print (len(set(MINEreactions)))
        self.assertEqual(len(MINEreactions), 39)
       	self.assertIn('rat000006_c0', MINEreactions)
        self.assertIn('R00045_c0', MINEreactions)
        
        allrxns = DB.get_all_reactions()
        self.assertIn('rat000006_c0', allrxns)
        self.assertIn('R00045_c0', allrxns)
        MINEcpds = DB.get_compounds_in_model('ATLAS')
        
        self.assertIn('C00007_c0', MINEcpds)

        allcpds = DB.get_all_compounds()
        self.assertIn('C00007_c0', allcpds)
        self.assertIn('C04547_c0', allcpds)

        cf = DB.get_cpd_chemicalformula('C00007_c0')
        self.assertEqual(cf, None)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00007_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C04547_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00755_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12361_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00003_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12352_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00004_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00080_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

    def test_ATLAS_inchi(self):
        batlasdb.build_atlas(PATH+'/data_atlas', PATH+'/datam/testPATRICinchi_atlas.db', True, 1, 'bio')
        DB = Q.Connector(PATH+'/datam/testPATRICinchi_atlas.db')
        cnx = sqlite3.connect(PATH+'/datam/testPATRICinchi_atlas.db')
        QC = cnx.execute("""SELECT * FROM model WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('ATLAS', results)

        QC = cnx.execute("""SELECT * FROM cluster WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('3', results)

        MINEreactions = DB.get_reactions_in_model('ATLAS')
        self.assertEqual(len(MINEreactions), 39)
        self.assertIn('rat000006_c0', MINEreactions)
        self.assertIn('R00045_c0', MINEreactions)
        
        allrxns = DB.get_all_reactions()
        self.assertIn('rat000006_c0', allrxns)
        self.assertIn('R00045_c0', allrxns)
        MINEcpds = DB.get_compounds_in_model('ATLAS')
        
        self.assertIn('C00007_c0', MINEcpds)

        allcpds = DB.get_all_compounds()
        self.assertIn('C00007_c0', allcpds)
        self.assertIn('C04547_c0', allcpds)

        inchis = DB.get_all_compounds_inchi()
        self.assertIn('InChI=1S/O2/c1-2', inchis)
        self.assertIn('InChI=1S/C16H16O4/c1-19-15-9-11(5-7-13(15)17)3-4-12-6-8-14(18)16(10-12)20-2/h3-10,17-18H,1-2H3/b4-3+', inchis)

        cf = DB.get_cpd_chemicalformula('C04547_c0')
        self.assertEqual(cf, 'C16H16O4')

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00007_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C04547_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00755_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12361_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12361_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12352_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00003_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C00004_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)


if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/datam/testPATRIC_atlas.db')
    os.remove(PATH+'/datam/testPATRICinchi_atlas.db')
