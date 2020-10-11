from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'leanne382@gmail.com'
__description__ = 'Code contains database functions'

import pubchempy
import re
import os
from copy import deepcopy
try:
    import urllib.request as urllib2
    import http.client as httplib
except:
    import urllib2
    import httplib
PATH = os.path.dirname(os.path.abspath(__file__))
print (PATH)
KEGG = 'http://rest.kegg.jp/'
translate_dict = {}
if os.path.exists(PATH+"/KEGG2INCHI.txt") is True:
    
    with open(PATH+"/KEGG2INCHI.txt", "r") as fin:
        for line in fin:
            line = line.strip()
            larray = line.split('\t')
            translate_dict[larray[0]] = str(larray[1])+"\t"+str(larray[2])
translatefile = open(PATH+"/KEGG2INCHI.txt", "a")
    

def verbose_print(verbose, line):
    '''verbose print function'''
    if verbose:
        print(line)

def extract_KEGG_data(url, verbose):
    '''Extract Kegg db info'''
    verbose_print(verbose, "STATUS: Getting data for url "+url)
    try:
        verbose_print(verbose, url)
        data = urllib2.urlopen(url).read()
        darray = str(data).split('\\n')
        return darray
    except urllib2.HTTPError:
        return None

def kegg2pubcheminchi(cpd, verbose):
    '''Convvert kegg ID to InChI value'''
    if cpd in translate_dict:
        larray = translate_dict[cpd].split("\t")
        for i, x in enumerate(larray):
            larray[i] = None if larray[i] == 'None' else larray[i]
        return(larray[0], larray[1])
    
    else:
        darray = extract_KEGG_data(KEGG+'get/'+cpd, verbose)

        inchicpd = None
        cas = None
        if darray:
            for value in darray:
                array = value.split()
                if 'PubChem:' in array:
                    index = array.index('PubChem:')
                    sid = array[index+1]
                    try:
                        substance = pubchempy.Substance.from_sid(sid)
                        substance_cids = substance.cids
                        if substance_cids:
                            try:
                                compounds = pubchempy.get_compounds(substance_cids[0])
                                if compounds:
                                    inchicpd = compounds[0].inchi
                            except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                                pass
                    except (pubchempy.PubChemHTTPError, httplib.BadStatusLine, urllib2.URLError):
                        verbose_print(verbose, 'WARNING: Could not get substance for {} {}'.format(sid, cpd))
                        pass
                if 'CAS:' in array:
                    index = array.index('CAS:')
                    cas = array[index+1]
        translatefile.write(cpd+"\t"+str(inchicpd)+"\t"+str(cas)+"\n")
        return (inchicpd, cas)

def get_KEGG_IDs(ID, compartment, KEGGdict):
    '''Retrieve KEGG IDs'''
    new_ID = re.sub('_'+compartment+'$', '', ID)
    try:
        KEGG_ID = KEGGdict[new_ID]
        original_KEGG_ID = str(deepcopy(KEGG_ID))
        if not KEGG_ID:
            KEGG_ID = ID+'0'
        else:
            KEGG_ID = str(KEGG_ID)+'_'+str(compartment)+'0'
    except KeyError:
        KEGG_ID = ID+'0'
        original_KEGG_ID = None
    return(KEGG_ID, original_KEGG_ID)

def open_translation_file(file_name):
    '''opens and stores KEGG translation files '''
    dictionary = {}
    with open(file_name) as fin:
        for line in fin:
            line = line.strip()
            larray = line.split('\t')
            try:
                KEGGIDS = larray[1].split('|')
                dictionary[larray[0]] = KEGGIDS[0]
            except IndexError:
                dictionary[larray[0]] = None
    return dictionary