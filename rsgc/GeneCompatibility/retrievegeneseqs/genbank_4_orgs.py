__author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'Gets genome for chassis organism from NCBI and generates BLAST database for genome'

import os
import re
import pickle
from tqdm import tqdm
from rsgc.Database import query as Q
from rsgc.Database import mackinac
from rsgc.GeneCompatibility.KEGG import kegg_functions as kf

PATH = os.path.dirname(os.path.abspath(__file__))

def get_database_organism_genbank_ids(database=False, listoforganisms=False, output_directory='', outputfile=''):
    '''get genbank accessions for organisms in database 
    so phylogenetic analysis can be done'''

    db_org_gbs={}

    if database:
        genbank_file = re.sub('db$', 'gb', database)
        DB = Q.Connector(database)
        db_orgs = DB.get_all_models()
    elif listoforganisms:
        db_orgs = listoforganisms
        genbank_file = outputfile

    if os.path.isfile(genbank_file) is False:
        print ('STATUS:\tGetting genbank accessions for organisms in RetSynth database ...')
        for org in tqdm(db_orgs):
            ###get genebank accession for genomes acquired from PATRIC/mackinac 
            if re.search('^\d+', org[0]):
            
                tempdict = mackinac.get_genome_summary(org[0])
            
                try:
                    db_org_gbs[org[0]] = [tempdict['assembly_accession']]
            
                except:
                    print ('WARNING:\tNo genbank accession for {}'.format(org[0]))        
            ###get genebank accession for genomes acquired KEGG
            elif org[1].endswith('_KEGG'):

                darray = kf.extract_KEGG_data('get/gn:'+org[0])
                gb = kf.get_taxid_gb_info(darray)
                db_org_gbs[org[0]] = gb
        
        with open(os.path.join(output_directory, genbank_file, 'wb')) as fout:
            pickle.dump(db_org_gbs, fout)
        
        return(db_org_gbs)

    else:

        with open(os.path.join(output_directory, genbank_file), 'rb') as fin:
            db_org_gbs = pickle.load(fin)

        return(db_org_gbs)