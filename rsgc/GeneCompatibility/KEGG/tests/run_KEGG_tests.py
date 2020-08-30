from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests GenBank module'

import os
import re
import unittest
from KEGG import kegg_functions as kf

PATH = os.path.dirname(os.path.abspath(__file__))
PPATH = re.sub('tests', '', PATH)

class RunKEGGFunctions(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""
    
    def test_extract_KEGG_orgIDs(self):
        """test function extract_KEGG_orgIDs"""

        orgs_gbs, keggorganisms = kf.extract_KEGG_orgIDs({}, taxid_info=False, kegg_id=True)

        self.assertIn('T00700', keggorganisms.keys())
        self.assertIn('T00007', keggorganisms.keys())

        orgs_gbs, keggorganisms = kf.extract_KEGG_orgIDs({}, taxid_info=True, kegg_id=True, kegg_org_types='bacteria')

        self.assertEqual(keggorganisms['T00700']['sbc'][0], 'GCA_000020185.1')
        self.assertEqual(keggorganisms['T00007']['eco'][0], 'GCF_000005845.2')
        self.assertIn('sbc', orgs_gbs)
        self.assertIn('eco', orgs_gbs)

    def test_extract_KEGG_data(self):
        """test function extract_KEGG_data"""

        darray = kf.extract_KEGG_data('get/gn:T00007')
        self.assertTrue(darray)
        darray = kf.extract_KEGG_data('get/gn:T00700')
        self.assertTrue(darray)
    
    def test_get_taxid_gb_info(self):
        """test function get_taxid_gb_info"""
 
        darray = kf.extract_KEGG_data('get/gn:T00007')
        gb = kf.get_taxid_gb_info(darray)
        self.assertIn('GCF_000005845.2', gb)
        self.assertEqual(len(gb), 1)

        darray = kf.extract_KEGG_data('get/gn:T00700')
        gb = kf.get_taxid_gb_info(darray)
        self.assertIn('GCA_000020185.1', gb)
        self.assertEqual(len(gb), 1)

if __name__ == '__main__':
    unittest.main(exit=False)
