_author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'Gets 16/18S info from NCBI ftp site for organisms in KEGG and RS database'

import re
import os
import ftplib
import urllib
import gzip
import shutil
import zipfile
import glob
from tqdm import tqdm
from GeneCompatibility.KEGG import kegg_functions as kf
PATH = os.path.dirname(os.path.abspath(__file__))

def verbose_print(verbose, line):
    if verbose:
        print(line)

class NCBI_SSU(object):
    def __init__(self, orgs_gb, outputfile_name, wipe_folder=True, verbose=False, run=True):
        '''initialize'''

        self.verbose = verbose
        self.orgs_gb = orgs_gb

        ##CHECK TO SEE IF NECESSARY .fa FILE HAS BEEN GENERATED###
        if run:
            if os.path.isfile(PATH+'/data/'+outputfile_name) is False:
    
                print ('STATUS:\tretrieving 16S data for orgs in database and KEGG ...')
    
                self.outputfile = open(PATH+'/data/'+outputfile_name, 'w')

                self.ASSEM_SUM = {}
                self.ASSEM_SUM.setdefault('GCA', {})
                self.ASSEM_SUM.setdefault('GCF', {})

                if os.path.isdir(PATH+'/data/assemblysummaries/') is False:
                    self.unzip_files(PATH+'/data/assemblysummaries.zip')

                if os.path.isdir(PATH+'/data/defaultfiles/') is False:
                    self.unzip_files(PATH+'/data/defaultfiles.zip')

                self.default_SSU = {}
                self.load_SSU_from_default_file()

                self.process_assembly_summary(PATH+'/data/assemblysummaries/assembly_summary_genbank.txt')
                self.process_assembly_summary(PATH+'/data/assemblysummaries/assembly_summary_refseq.txt')

                self.get_genome_zip_files()

                self.outputfile.close()

                if wipe_folder:
                    self.wipe_ncbi_gn_data_folder()

            else: 
                print ('STATUS:\toutput file {} already exists therefore skipping loading 16S sequences from NCBI FTP '.format(outputfile_name))


    def fill_assembly_dictionary(self, larray, index_gca, index_gcf, index_ftp):
        '''load information from assembly file into dictionary'''

        self.ASSEM_SUM['GCA'].setdefault(larray[index_gca], {})
        self.ASSEM_SUM['GCA'][larray[index_gca]].setdefault('ftp_path', larray[index_ftp])
        
        self.ASSEM_SUM['GCF'].setdefault(larray[index_gcf], {})
        self.ASSEM_SUM['GCF'][larray[index_gcf]].setdefault('ftp_path', larray[index_ftp])

    def process_assembly_summary(self, infile):
        '''pull genome IDS from assembly summary files and load into dictionary'''

        print ('STATUS:\tprocess {}'.format(infile))

        with open(infile) as fin:
    
            header = fin.readline()

            for line in fin:

                if line.startswith('#'):
                    lheader = line.split('\t')
                    index_gca = lheader.index('# assembly_accession')
                    index_gcf = lheader.index('gbrs_paired_asm')
                    index_ftp = lheader.index('ftp_path')

                else:
                    larray = line.split('\t')
                    self.fill_assembly_dictionary(larray, index_gca, index_gcf, index_ftp)

    def get_genome_zip_files(self):
        '''implement getting genome 16S RNA for organisms in database and kegg'''

        print ('STATUS:\tgetting rna files from NCBI ftp site ...')
        for org in tqdm(self.orgs_gb):
            if self.orgs_gb[org]:
                
                accessionid = self.orgs_gb[org][0]

                if accessionid in self.default_SSU:
                    if self.default_SSU[accessionid] and self.default_SSU[accessionid] != '':
                        self.outputfile.write('>'+accessionid+'\n')
                        self.outputfile.write(self.default_SSU[accessionid]+'\n')
                
                else:
                    if accessionid.startswith('GCA'):
                        self.extract_rna_genomic_file_NCBI(accessionid, self.ASSEM_SUM['GCA'])

                    elif accessionid.startswith('GCF'):
                        self.extract_rna_genomic_file_NCBI(accessionid, self.ASSEM_SUM['GCF'])

                    else:
                        print (self.orgs_gb[org])
                        print ('WARNING:\t{} does not start with GCF or GCA so it is getting skipped ... '.format(accessionid))

    def extract_rna_genomic_file_NCBI(self, genomeID, ASSEMDICT):
        '''Pull data from NCBI ftp'''

        FTPERROR = False

        ###CHECKS TO SEE IF GENOME ID IS IN OUR ASSEMBLY DICTIONARY###
        try:

            ###GET FTP SITE###
            ftp_url = ASSEMDICT[genomeID]['ftp_path']
            lftp_url = ftp_url.split('/')
            filename = lftp_url[-1]+'_rna_from_genomic.fna.gz'
            url = ftp_url+'/'+filename

            ###PATHWAYS TO .gz and .fna FILES###
            source_filepath = PATH+'/ncbi_gn_data/'+filename
            dest_filepath = PATH+'/ncbi_gn_data/'+filename
            dest_filepath = re.sub('.gz', '', dest_filepath)

            ###CHECKS TO SEE IF FNA FILE FOR GENOME IS PRESENT IN FOLDER###
            if os.path.isfile(dest_filepath) is False:
 
                ###CHECKS TO SEE IF GZ FILE FOR GENOME IS PRESENT IN FOLDER###
                if os.path.isfile(source_filepath) is False:

                    ###ATTEMPTS TO PULL FILE FROM FTP SITE###
                    try:
                        file_retrieved = urllib.request.urlretrieve(url, source_filepath)

                    except urllib.error.URLError:
                        verbose_print(self.verbose, 'WARNING:\t{} could not be loaded cause file does not exist'.format(url))
                        self.fill_default_noinfo_file(genomeID)
                        FTPERROR = True

                ###IF PULLING gz FILE WAS SUCESSFUL PULL 16/18S INFO FROM FROM .fna FILE###
                if FTPERROR is False:
                    block_size = 65536

                    ###GENERATES .fna FILE FROM .gz FILE###
                    with gzip.open(source_filepath, 'rb') as s_file, \
                            open(dest_filepath, 'wb') as d_file:
                        shutil.copyfileobj(s_file, d_file, block_size)

                    ###PULL INFO FROM .fna FILE###
                    self.pull_SSU_from_fna_file(dest_filepath, genomeID)

                    ###REMOVE .gz FILE###
                    os.remove(source_filepath)
            else:
                self.pull_SSU_from_fna_file(dest_filepath, genomeID)

        except KeyError:
            print ('STATUS:\t{} not in assembly summary file'.format(genomeID))
            self.fill_default_noinfo_file(genomeID)


    def pull_SSU_from_fna_file(self, filename, genomeID):
        '''extract SSU information from loaded fasta files from the FTP load'''

        fasta = {}
        with open(filename) as file_one:

            for line in file_one:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    active_sequence_name = line[1:]

                    if active_sequence_name not in fasta:
                        fasta[active_sequence_name] = []
                    continue

                sequence = line
                fasta[active_sequence_name].append(sequence)

        for k in fasta:

            if (re.search('\[product=16S ribosomal RNA\]', k) or re.search('\[product=18S ribosomal RNA\]', k)):

                self.outputfile.write('>'+str(genomeID)+'\n')
                self.outputfile.write(''.join(fasta[k])+'\n')
                self.append_data_to_default_file('>'+str(genomeID)+'\n', ''.join(fasta[k])+'\n')
                break

    def wipe_ncbi_gn_data_folder(self):
        '''empty folder of genome rna files'''

        list_files = glob.glob(PATH+'/ncbi_gn_data/*')
        
        for file in list_files:

            if file.endswith('gz') or file.endswith('fna'):

                os.remove(file)
    
    def append_data_to_default_file(self, header, seq):
        '''append new SSU data to current default file (will speed up
        subsequent runs after initial)'''

        with open(PATH+'/data/defaultfiles/default_SSU.fa', 'a') as fout:
            fout.write(header)
            fout.write(seq)

    def load_SSU_from_default_file(self):
        '''load in SSUs stored in default file'''

        with open(PATH+'/data/defaultfiles/default_SSU.fa') as fin:
 
            for line in fin:
 
                line = line.strip()
 
                if line.startswith('>'):
                    key = re.sub('>', '', line)
                
                else:
                    self.default_SSU[key] = line
        
        with open(PATH+'/data/defaultfiles/default_noinfo.txt') as fin:
 
            for line in fin:
 
                line = line.strip()

                self.default_SSU[line] = None
 
    def fill_default_noinfo_file(self, genomeID):
        '''fill a defualt file stores all ids that have no rna_genome file
        (will save time in future runs)'''

        with open(PATH+'/data/default_noinfo.txt', 'a') as fout:
            fout.write(genomeID+'\n')

    def unzip_files(self, file_to_unzip):

        zipref = zipfile.ZipFile(file_to_unzip, 'r')
        zipref.extractall(PATH+'/data/.')