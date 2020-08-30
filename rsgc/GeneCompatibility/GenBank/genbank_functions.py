__author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'GenBank functions'

import urllib.request, urllib.error, urllib.parse, http.client, re

GenBank = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&id='

def extract_GenBank_data(GenBankID, trycount=0):
    '''Pull data from GenBank GenBankID'''
    GenBankID = GenBank+GenBankID
    try:
        data = urllib.request.urlopen(GenBankID).read().decode('utf-8')
        darray = data.split('\n')
        return darray

    except:
        if trycount == 3:
            return None
        else:
            return extract_GenBank_data(GenBankID, trycount+1)

def extract_GenBank_sequence(GenBankID):
    '''Extract gene sequence from GenBank ID'''

    GenBankData = extract_GenBank_data(GenBankID)

    if GenBankData is None:
        print('WARNING:\tNCBI has no information for GenBank ID %s' % GenBankID)
        return 'no_gene_info'
    else:
        idxstart = 0
        try:
            # Searching for
            # ORIGIN
            #         1 agcaatatta caactaacaa ...
            #        61 ...

            # while not re.search('translation', GenBankData[idxstart]):
            #     idxstart += 1
            # idxstart += 1
            # idxend = idxstart + 1
            # while GenBankData[idxend].startswith(' '):
            #     idxend += 1
            # indent = len(re.search('\s+\d+\s', GenBankData[idxstart]).group(0))

            # protseq = ''.join([''.join(d[indent:].split(' ')) for d in GenBankData[idxstart:idxend]])
            # print (protseq)
            protseq = ''
            START = False
            for line in GenBankData:
                if re.search('translation=', line):
                    line = line.strip()
                    protseq = re.sub('\/translation=','',line)
                    START=True
                elif line.startswith('ORIGIN'):
                    START=False
                elif START is True:
                    line = line.strip()
                    protseq=protseq+line 

            while not GenBankData[idxstart].startswith('ORIGIN'):
                idxstart += 1
            idxstart += 1
            idxend = idxstart + 1
            while GenBankData[idxend].startswith(' '):
                idxend += 1
            indent = len(re.search('\s+\d+\s', GenBankData[idxstart]).group(0))

            geneseq = ''.join([''.join(d[indent:].split(' ')) for d in GenBankData[idxstart:idxend]])
            return geneseq, protseq




        except IndexError:
            print('WARNING:\tNCBI has no sequence information for GenBank ID %s' % GenBankID)
            return 'no_gene_info'

