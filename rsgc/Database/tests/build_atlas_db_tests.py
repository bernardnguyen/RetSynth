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
from rsgc.Database import build_ATLAS_db as batlasdb

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('/tests', '', PATH)
init_db.Createdb(PATH+'/kbasetestadd.db', False)
bkdb.BuildKbase(PATH+'/datam', PPATH+'/data/KbasetoKEGGCPD.txt', PPATH+'/data/KbasetoKEGGRXN.txt',
                False, PATH+'/kbasetestadd.db', 'bio')

init_db.Createdb(PATH+'/kbasetestaddinchi.db', True)
bkdb.BuildKbase(PATH+'/data', PPATH+'/data/KbasetoKEGGCPD.txt', PPATH+'/data/KbasetoKEGGRXN.txt',
                True, PATH+'/kbasetestaddinchi.db', 'bio')

class BuildATLAStests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def test_ATLAS_no_inchi(self):
        ''' test addition of reaction info from MINE raw files'''
        batlasdb.build_atlas(PATH+'/data4', PATH+'/kbasetestadd.db', False, 1, 'bio')
        DB = Q.Connector(PATH+'/kbasetestadd.db')
        cnx = sqlite3.connect(PATH+'/kbasetestadd.db')
        QC = cnx.execute("""SELECT * FROM model WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('ATLAS', results)

        QC = cnx.execute("""SELECT * FROM cluster WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('18', results)

        MINEreactions = DB.get_reactions_in_model('ATLAS')
        self.assertEqual(len(MINEreactions), 36)
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
        self.assertEqual(cf, 'None')

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
        batlasdb.build_atlas(PATH+'/data4', PATH+'/kbasetestaddinchi.db', True, 1, 'bio')
        DB = Q.Connector(PATH+'/kbasetestaddinchi.db')
        cnx = sqlite3.connect(PATH+'/kbasetestaddinchi.db')
        QC = cnx.execute("""SELECT * FROM model WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('ATLAS', results)

        QC = cnx.execute("""SELECT * FROM cluster WHERE ID = ?""", ("ATLAS",))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertIn('3', results)

        MINEreactions = DB.get_reactions_in_model('ATLAS')
        self.assertEqual(len(MINEreactions), 46)
        self.assertIn('rat000006_c0', MINEreactions)
        self.assertIn('R00045_c0', MINEreactions)
        
        allrxns = DB.get_all_reactions()
        self.assertIn('rat000006_c0', allrxns)
        self.assertIn('R00045_c0', allrxns)
        MINEcpds = DB.get_compounds_in_model('ATLAS')
        
        self.assertIn('InChI=1S/O2/c1-2_c0', MINEcpds)

        allcpds = DB.get_all_compounds()
        self.assertIn('InChI=1S/O2/c1-2_c0', allcpds)
        self.assertIn('InChI=1S/C16H16O4/c1-19-15-9-11(5-7-13(15)17)3-4-12-6-8-14(18)16(10-12)20-2/h3-10,17-18H,1-2H3/b4-3+_c0', allcpds)

        cf = DB.get_cpd_chemicalformula('InChI=1S/C16H16O4/c1-19-15-9-11(5-7-13(15)17)3-4-12-6-8-14(18)16(10-12)20-2/h3-10,17-18H,1-2H3/b4-3+_c0')
        self.assertEqual(cf, 'C16H16O4')

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("InChI=1S/O2/c1-2_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("InChI=1S/C16H16O4/c1-19-15-9-11(5-7-13(15)17)3-4-12-6-8-14(18)16(10-12)20-2/h3-10,17-18H,1-2H3/b4-3+_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("InChI=1S/C8H8O3/c1-11-8-4-6(5-9)2-3-7(8)10/h2-5,10H,1H3_c0","R00043_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12361_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("InChI=1S/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p+1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("C12352_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("InChI=1S/C21H29N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1,3-4,7-8,10-11,13-16,20-21,29-32H,2,5-6H2,(H2,23,33)(H,34,35)(H,36,37)(H2,22,24,25)/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)

        QC = cnx.execute("""SELECT * FROM reaction_compound WHERE cpd_ID = ? and reaction_ID = ?""", ("InChI=1S/p+1_c0","rat000008_c0"))
        hits = QC.fetchall()
        results = [i[0] for i in hits]
        self.assertEqual(len(results), 1)


if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/kbasetestadd.db')
    os.remove(PATH+'/kbasetestaddinchi.db')
