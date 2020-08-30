__author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'KEGG functions'

import urllib.request, urllib.error, urllib.parse, http.client, re
from tqdm import tqdm

KEGG = 'http://rest.kegg.jp/'

def extract_KEGG_data(url, trycount=0):
    '''Pull data from KEGG url'''
    url = KEGG+url
    try:
        data = urllib.request.urlopen(url).read().decode('utf-8')
        darray = data.split('\n')
        return darray

    except:
        if trycount == 3:
            return None
        else:
            return extract_KEGG_data(url, trycount+1)

def extract_KEGG_orgIDs(db_org_gbs, taxid_info=False, kegg_id=True, kegg_org_types='all'):
    '''Retrieve bacteria organism IDs from KEGG'''

    orgIDs_gb = db_org_gbs
    keggorganisms = {}

    if kegg_org_types == 'all':
        targetorganisms = ['bacteria','plants', 'fungi']
    else:
        targetorganisms = [kegg_org_types]
    # targetorganisms  = ['fungi']
    darray_org = extract_KEGG_data('list/organism')

    for line in tqdm(darray_org):

        array = line.split('\t')

        for to in targetorganisms:
            try:
                if to in array[3].lower():
                    
                    if taxid_info:
                        ###get taxonomic id and genebank accession number for kegg ID
                        darray_genome = extract_KEGG_data('get/gn:'+array[0])
                        gb = get_taxid_gb_info(darray_genome)


                        if array[1] not in orgIDs_gb:
                            orgIDs_gb.setdefault(array[1], gb)
                    
                    if kegg_id:
                        keggorganisms.setdefault(array[0], {})

                        if taxid_info:
                            keggorganisms[array[0]].setdefault(array[1], gb)

                        else:
                            keggorganisms[array[0]].setdefault(array[1], [])

                        keggorganisms[array[0]].setdefault('org_type', to)

                    break
 
            except IndexError:
                pass

    return orgIDs_gb, keggorganisms

def get_taxid_gb_info(darray):
    '''get taxonomy id and genbank accession for kegg genome'''

    gb = []

    if darray:

        for line in darray:

            # if line.startswith('TAXONOMY'):

            #     larray = line.split()
            #     taxid = re.sub('TAX:', '', larray[1])
            
            match = re.search('GenBank', line)
            if match:

                line = line.strip()
                match_search = re.search('Assembly:(\S+)\)', line)

                try:
                    gb.append(match_search.group(1))

                except AttributeError:
                    pass
        if not gb:
            for line in darray:

                # if line.startswith('TAXONOMY'):

                #     larray = line.split()
                #     taxid = re.sub('TAX:', '', larray[1])
                
                match = re.search('RefSeq', line)
                if match:

                    line = line.strip()
                    match_search = re.search('Assembly:(\S+)\)', line)

                    try:
                        gb.append(match_search.group(1))

                    except AttributeError:
                        pass


    return gb
