from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Tests genbak_4_orgs module'

import os
import re
import pickle
import unittest
from retrievegeneseqs import genbank_4_orgs as go
PATH = os.path.dirname(os.path.abspath(__file__))

orgs = [['536056.3', 'Escherichia coli DH1'],
        ['953739.5', 'Streptomyces venezuelae ATCC 10712'],
        ['eco', 'eco_KEGG'], ['sbc', 'sbc_KEGG'], ['ecj', 'ecj_KEGG']]

class RunGenbank4orgs(unittest.TestCase):
    def setUp(self):
        """Initialize before every test."""
        print ("Initializing tests")

    def tearDown(self):
        """Clean up after each test."""

    def test_get_database_organism_genbank_ids(self):
        """test function get_database_organism_genbank_ids"""
        orgs_gb =  go.get_database_organism_genbank_ids(database=False, listoforganisms=orgs, outputfile=PATH+'/temp.list')

        self.assertEqual(orgs_gb['953739.5'], ['GCA_000253235.1'])
        self.assertEqual(orgs_gb['536056.3'], ['GCA_000023365.1'])
        self.assertEqual(orgs_gb['eco'], ['GCF_000005845.2'])
        self.assertEqual(orgs_gb['sbc'], ['GCA_000020185.1'])
        self.assertEqual(orgs_gb['ecj'], ['GCA_000010245.1'])

if __name__ == '__main__':
    unittest.main(exit=False)
    os.remove(PATH+'/temp.list')