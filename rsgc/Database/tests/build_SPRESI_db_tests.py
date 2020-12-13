from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests RDF reader'
import sqlite3
import os
import re
import unittest
from shutil import copyfile
from rsgc.Database import query as Q
from rsgc.Database import build_SPRESI_db as bspresidb

PATH = os.path.dirname(os.path.abspath(__file__))

copyfile(PATH+'/datam/testPATRIC.db', PATH+'/datam/testPATRIC_spresi.db')
copyfile(PATH+'/datam/testPATRICinchi.db', PATH+'/datam/testPATRICinchi_spresi.db')

DB = Q.Connector(PATH+'/datam/testPATRIC_spresi.db')
compartmentID = DB.get_compartment('cytosol')
compartmentID = compartmentID[0]

bspresidb.RDF_Reader(PATH+'/data_spresi/', PATH+'/datam/testPATRIC_spresi.db', 'chem', compartmentID, 1, temp_option=False, pressure_option=False,
                          yield_option=False, time_option=False, catalyst_option=False, solvent_option=False)

DBinchi = Q.Connector(PATH+'/datam/testPATRICinchi_spresi.db')
bspresidb.RDF_Reader(PATH+'/data_spresi/', PATH+'/datam/testPATRICinchi_spresi.db', 'chem', compartmentID, 1, temp_option=False, pressure_option=False,
                          yield_option=False, time_option=False, catalyst_option=False, solvent_option=False)



class RDFileReaderTests(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
        print ("Clearing out test suite")

    def tests_RDFileReader(self):
        '''Tests RDFReader'''
        print ("Tests RDFReader")
        cnx = sqlite3.connect(PATH+'/datam/testPATRIC_spresi.db')

        '''Reaction check'''
        Q = cnx.execute("SELECT reaction_ID FROM model_reaction WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 6)
        self.assertIn(('rxn1_s',), result)
        self.assertIn(('rxn2_s',), result)
        self.assertIn(('rxn3_s',), result)
        self.assertIn(('rxn4_s',), result)
        self.assertIn(('rxn5_s',), result)
        self.assertIn(('rxn6_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn1_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn1_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn2_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn2_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn3_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn3_s',), result)

        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn5_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn5_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn6_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn6_s',), result)

        '''Reversibility check'''
        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn1_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn1_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn2_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn2_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn3_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn3_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn4_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn5_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn5_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn6_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn6_s', False), result)

        '''Compound check'''
        Q = cnx.execute("SELECT cpd_ID FROM model_compound WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 14)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1C=CC(N)=CC=1_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1C=CC(N)=CC=1_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1=CC2=C(N3CN(C4=C(C3)C=C(C=C4)C)C2)C=C1_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1=CC2=C(N3CN(C4=C(C3)C=C(C=C4)C)C2)C=C1_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1(COB(OC1)C1C=CC=CC=1)C_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1(COB(OC1)C1C=CC=CC=1)C_c0',),
                         result)
        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1(COB(OC1)C1C=CC=CC=1)C_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1(COB(OC1)C1C=CC=CC=1)C_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('COC1=CC2C(C=C1)=CC=CC=2_c0',))
        result = Q.fetchone()
        self.assertEqual(('COC1=CC2C(C=C1)=CC=CC=2_c0',),
                         result)
 
        Q = cnx.execute("SELECT chemicalformula FROM compound WHERE ID = ?", ('COC1=CC2C(C=C1)=CC=CC=2_c0',))
        result = Q.fetchone()
        self.assertEqual(('C11H10O',), result)

        Q = cnx.execute("SELECT chemicalformula FROM compound WHERE ID = ?", ('C1C=CC(N)=CC=1_c0',))
        result = Q.fetchone()
        self.assertEqual(('C6H7N',), result)

        '''reaction_compound check'''
        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn1_s',))
        result = Q.fetchall()
        self.assertIn(('rxn1_s', 'CC1C=CC(N)=CC=1_c0', 0, 1, 0), result)
        self.assertIn(('rxn1_s', 'CC1=CC2=C(N3CN(C4=C(C3)C=C(C=C4)C)C2)C=C1_c0', 1, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn3_s',))
        result = Q.fetchall()
        self.assertIn(('rxn3_s', 'CC1(COB(OC1)C1C=CC=CC=1)C_c0', 0, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn4_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4_s', 'CC1(C(ON=C1C1C=CC=CC=1)=O)C_c0', 0, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn6_s',))
        result = Q.fetchall()
        self.assertIn(('rxn6_s', 'CCOC(=O)C1=C(Cl)C(C(OCC)=O)=C2C1=CC=CC=C2C1C2C(C=CC=1)=CC=CC=2_c0', 0, 1, 0), result)

        '''Catalysts check'''
        Q = cnx.execute("SELECT * FROM reaction_catalysts WHERE reaction_ID = ?", ('rxn1_s',))
        result = Q.fetchall()
        self.assertIn(('rxn1_s', 'InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)', 'None'), result)

        Q = cnx.execute("SELECT * FROM reaction_catalysts")
        result = Q.fetchall()
        self.assertEqual(len(result), 7)

        Q = cnx.execute("SELECT * FROM reaction_solvents")
        result = Q.fetchall()
        self.assertEqual(len(result), 3)

        '''Solvents check'''
        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?", ('rxn3_s',))
        result = Q.fetchall()
        self.assertIn(('rxn3_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn3_s',))
        result = Q.fetchall()
        self.assertIn(('rxn3_s', '120 degree', 'None', '12 h', '93.0-93.0',
                       'JOURNAL ARTICLE Tobisu Mamoru; Shimasaki Toshiaki; Chatani Naoto; Angew. Chem., Int. Ed. +Engl., 2008, Vol. 47, P. 4866-4869'), result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn4_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4_s', '25 degree,80 degree,None, 600 degree', 'None,None,None,None',
                       '1 h, 30 min,None,None', '44.0-44.0',
                       'JOURNAL ARTICLE Begue Didier; Dargelos Alain; Berstermann Hans M.; Netsch Klaus P.; Bedna+rek Pawel; Wentrup Curt; J. Org. Chem., 2014, Vol. 79, P. 1247-1253'), result)
    
        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn5_s',))
        result = Q.fetchall()
        self.assertIn(('rxn5_s', 'None', 'None', 'None', 'None',
                       'PATENT CATALYST AND PROCESS FOR THE SELECTIVE DIMERIZATION OF PROPYLENE TO METHY+L-1-PENTENE, STEVENS JAMES C.; FORDYCE WILLIAM A., 1991, Patent number-4855523, US, Patent Class-4 C 07 C 2/10, Patent Owner-THE DOW CHEMICAL CO.'), result)

    def tests_RDFileReader_inchi(self):
        '''Tests RDFReader with inchi database'''
        print ("Tests RDFReader with inchi database")
        cnx = sqlite3.connect(PATH+'/datam/testPATRICinchi_spresi.db')

        '''Reaction check'''
        Q = cnx.execute("SELECT reaction_ID FROM model_reaction WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 6)
        self.assertIn(('rxn1_s',), result)
        self.assertIn(('rxn2_s',), result)
        self.assertIn(('rxn3_s',), result)
        self.assertIn(('rxn4_s',), result)
        self.assertIn(('rxn5_s',), result)
        self.assertIn(('rxn6_s',), result)

        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn1_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn1_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn2_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn2_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn3_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn3_s',), result)

        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn4_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn4_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn5_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn5_s',), result)
        Q = cnx.execute("SELECT ID FROM reaction WHERE ID = ?", ('rxn6_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn6_s',), result)

        '''Reversibility check'''
        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn1_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn1_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn2_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn2_s', False), result)

        Q = cnx.execute("SELECT * FROM reaction_reversibility WHERE reaction_ID= ?",
                        ('rxn3_s',))
        result = Q.fetchone()
        self.assertEqual(('rxn3_s', False), result)

        '''Compound check'''
        Q = cnx.execute("SELECT cpd_ID FROM model_compound WHERE model_ID = ?", ('SR1',))
        result = Q.fetchall()
        self.assertEqual(len(result), 14)

        Q = cnx.execute("SELECT inchistring FROM compound")
        inchis = Q.fetchall()

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1C=CC(N)=CC=1_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1C=CC(N)=CC=1_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1=CC2=C(N3CN(C4=C(C3)C=C(C=C4)C)C2)C=C1_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1=CC2=C(N3CN(C4=C(C3)C=C(C=C4)C)C2)C=C1_c0',), result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1(COB(OC1)C1C=CC=CC=1)C_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1(COB(OC1)C1C=CC=CC=1)C_c0',),
                         result)
        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('CC1(COB(OC1)C1C=CC=CC=1)C_c0',))
        result = Q.fetchone()
        self.assertEqual(('CC1(COB(OC1)C1C=CC=CC=1)C_c0',),
                         result)

        Q = cnx.execute("SELECT ID FROM compound WHERE ID = ?",
                        ('COC1=CC2C(C=C1)=CC=CC=2_c0',))
        result = Q.fetchone()
        self.assertEqual(('COC1=CC2C(C=C1)=CC=CC=2_c0',),
                         result)
 
        Q = cnx.execute("SELECT chemicalformula FROM compound WHERE ID = ?", ('COC1=CC2C(C=C1)=CC=CC=2_c0',))
        result = Q.fetchone()
        self.assertEqual(('C11H10O',), result)

        Q = cnx.execute("SELECT chemicalformula FROM compound WHERE ID = ?", ('C1C=CC(N)=CC=1_c0',))
        result = Q.fetchone()
        self.assertEqual(('C6H7N',), result)


        self.assertIn(('InChI=1S/C7H9N/c1-6-2-4-7(8)5-3-6/h2-5H,8H2,1H3',), inchis)
        self.assertIn(('InChI=1S/C17H18N2/c1-12-3-5-16-14(7-12)9-18-11-19(16)10-15-8-13(2)4-6-17(15)18/h3-8H,9-11H2,1-2H3',),inchis)

        self.assertIn(('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3',),inchis)
        self.assertIn(('InChI=1S/C11H15BO2/c1-11(2)8-13-12(14-9-11)10-6-4-3-5-7-10/h3-7H,8-9H2,1-2H3',), inchis)
        
        '''reaction_compound check'''
        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn1_s',))
        result = Q.fetchall()
        self.assertIn(('rxn1_s', 'CC1C=CC(N)=CC=1_c0', 0, 1, 0), result)
        self.assertIn(('rxn1_s', 'CC1=CC2=C(N3CN(C4=C(C3)C=C(C=C4)C)C2)C=C1_c0', 1, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn3_s',))
        result = Q.fetchall()
        self.assertIn(('rxn3_s', 'CC1(COB(OC1)C1C=CC=CC=1)C_c0', 0, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn4_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4_s', 'CC1(C(ON=C1C1C=CC=CC=1)=O)C_c0', 0, 1, 0), result)

        Q = cnx.execute("SELECT * FROM reaction_compound WHERE reaction_ID = ?", ('rxn6_s',))
        result = Q.fetchall()
        self.assertIn(('rxn6_s', 'CCOC(=O)C1=C(Cl)C(C(OCC)=O)=C2C1=CC=CC=C2C1C2C(C=CC=1)=CC=CC=2_c0', 0, 1, 0), result)

        '''Catalysts check'''
        Q = cnx.execute("SELECT * FROM reaction_catalysts WHERE reaction_ID = ?",
                        ('rxn1_s',))
        result = Q.fetchall()
        self.assertIn(('rxn1_s', 'InChI=1S/C2HF3O2/c3-2(4,5)1(6)7/h(H,6,7)', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_catalysts")
        result = Q.fetchall()
        self.assertEqual(len(result), 7)

        Q = cnx.execute("SELECT * FROM reaction_solvents")
        result = Q.fetchall()
        self.assertEqual(len(result), 3)

        '''Solvents check'''
        Q = cnx.execute("SELECT * FROM reaction_solvents WHERE reaction_ID = ?", ('rxn3_s',))
        result = Q.fetchall()
        self.assertIn(('rxn3_s', 'InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3', 'None'),
                      result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn3_s',))
        result = Q.fetchall()
        self.assertIn(('rxn3_s', '120 degree', 'None', '12 h', '93.0-93.0',
                       'JOURNAL ARTICLE Tobisu Mamoru; Shimasaki Toshiaki; Chatani Naoto; Angew. Chem., Int. Ed. +Engl., 2008, Vol. 47, P. 4866-4869'), result)

        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn4_s',))
        result = Q.fetchall()
        self.assertIn(('rxn4_s', '25 degree,80 degree,None, 600 degree', 'None,None,None,None',
                       '1 h, 30 min,None,None', '44.0-44.0',
                       'JOURNAL ARTICLE Begue Didier; Dargelos Alain; Berstermann Hans M.; Netsch Klaus P.; Bedna+rek Pawel; Wentrup Curt; J. Org. Chem., 2014, Vol. 79, P. 1247-1253'), result)
        
        Q = cnx.execute("SELECT * FROM reaction_spresi_info WHERE reaction_ID = ?",
                        ('rxn5_s',))
        result = Q.fetchall()
        self.assertIn(('rxn5_s', 'None', 'None', 'None', 'None',
                       'PATENT CATALYST AND PROCESS FOR THE SELECTIVE DIMERIZATION OF PROPYLENE TO METHY+L-1-PENTENE, STEVENS JAMES C.; FORDYCE WILLIAM A., 1991, Patent number-4855523, US, Patent Class-4 C 07 C 2/10, Patent Owner-THE DOW CHEMICAL CO.'), result)

if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/datam/testPATRIC_spresi.db')
    os.remove(PATH+'/datam/testPATRICinchi_spresi.db')
