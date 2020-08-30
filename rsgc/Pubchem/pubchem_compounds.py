__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Gets compound ID from pubchem IDs or general names'

import re
import http.client
import urllib.request, urllib.error, urllib.parse
from copy import deepcopy
import pubchempy
class PubchemConnector(object):
    """Retrieves a compound's database ID from from pubchem ID or compound name"""
    def __init__(self, DB):
        '''Initalize class'''
        self.DB = DB

    def _check_synonyms(self, synonyms):
        '''
        Check all synonyms associated with a given a pubchem ID or general name
        '''
        count = 0
        for synonym in synonyms[0]['Synonym']:
            db_cpdID = self.DB.get_compound_ID(synonym)
            if db_cpdID == 'None':
                synonym = re.sub('-', '_', synonym)
                synonym = re.sub('\s+', '_', synonym)
                db_cpdID = self.DB.get_compound_ID(synonym)
                if db_cpdID != 'None':
                    break
                elif db_cpdID == 'None':
                    db_cpdID = self.DB.get_compound_ID(synonym+'_c0')
                    if db_cpdID != 'None':
                        break
                    elif db_cpdID == 'None':
                        db_cpdID = self.DB.get_compound_ID(synonym+'_c')
                        if db_cpdID != 'None':
                            break
                        else:
                            count += 1

            else:
                break
        if count == len(synonyms):
            return 'None'
        else:
            if isinstance(db_cpdID, list):
                db_cpdID.sort()
                if len(db_cpdID) == 1:
                    return db_cpdID[0][0]
                elif len(db_cpdID) > 1:
                    print ('WARNING:\tMultiple database IDs were retrieved')
                    IDs = []
                    for ID in db_cpdID:
                        IDs.append(ID[0])
                    return IDs
            else:
                return db_cpdID

    def get_ID_from_pubchemID(self, pubchemID, inchidb):
        '''
        Retrieve db_cpdID from pubchemID
        '''
        if not inchidb:
            try:
                synonyms = pubchempy.get_synonyms(pubchemID, 'cid')
                if synonyms:
                    return self._check_synonyms(synonyms)
                else:
                    return 'None'
            except (pubchempy.PubChemHTTPError, http.client.BadStatusLine, urllib.error.URLError):
                return 'None'
        else:
            try:
                c = pubchempy.get_compounds(pubchemID, 'cid')
                if c:
                    return c[0].inchi+'_c0'
                else:
                    return 'None'
            except (pubchempy.PubChemHTTPError, http.client.BadStatusLine, urllib.error.URLError):
                return 'None'

    def get_ID_from_name(self, name):
        '''
        Retrieve ID from general name
        '''
        db_cpdID = self.DB.get_compound_ID(name)
        orginal_name = deepcopy(name)
        if db_cpdID == 'None':
            name = re.sub('-', '_', name)
            name = re.sub('\s+', '_', name)
            name = re.sub(',', '_', name)
            db_cpdID = self.DB.get_compound_ID(name)
            if db_cpdID == 'None':
                db_cpdID = self.DB.get_compound_ID(name+'_c0')
                if db_cpdID == 'None':
                    db_cpdID = self.DB.get_compound_ID(name+'_c')

            if db_cpdID == 'None':
                synonyms = pubchempy.get_synonyms(orginal_name, 'name')
                if synonyms:
                    return self._check_synonyms(synonyms)
                else:
                    return 'None'
            elif len(db_cpdID) == 1:
                return db_cpdID[0][0]

            elif len(db_cpdID) > 1:
                print ('WARNING:\tMultiple database IDs were retrieved, automatically use first ID')
                return db_cpdID[0][0]
        else:
            return db_cpdID[0][0]
