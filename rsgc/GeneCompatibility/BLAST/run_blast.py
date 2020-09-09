__author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'Code runs blast on gene sequences for specified enzyme against chassis organism genome'

import os
import re
import subprocess
from sys import platform
PATH = os.path.dirname(os.path.abspath(__file__))

if platform == 'darwin' or platform == "linux" or platform == "linux2":
    _BLAST_PATH = PATH+'/ncbi-blast-2.9.0+_mac/bin/'
elif platform == "win32" or platform == "win64" or platform == "cygwin":
    _BLAST_PATH = PATH+'/ncbi-blast-2.9.0+_win/bin/'

def make_directory(folder):
    try:
        os.mkdir(folder)
    except:
        pass

class ImplementBLAST(object):
    def __init__(self, blastdb, Seq16S_filename, output_directory):
        '''Initialize'''
        
        ##PATH TO STORED CHASSIS ORGANISM BLAST DATABASE
        self.blastdb = blastdb
        self.Seq16S_filename = Seq16S_filename

        make_directory(os.path.join(output_directory, "query"))
        make_directory(os.path.join(output_directory, "blastoutput"))

        if platform=="cygwin":
            self.Seq16S_filename = "C:/cygwin64"+self.Seq16S_filename
            self.blastdb = "C:/cygwin64"+self.blastdb
            self._BLAST_QUERY_PATH="C:/cygwin64"+os.path.abspath(os.path.join(output_directory,'query'))
            self._BLAST_OUTPUT_PATH="C:/cygwin64"+os.path.abspath(os.path.join(output_directory,'blastoutput'))
        else:
            self._BLAST_QUERY_PATH=os.path.abspath(os.path.join(output_directory,'query'))
            self._BLAST_OUTPUT_PATH=os.path.abspath(os.path.join(output_directory,'blastoutput'))            

    def read_16S_fasta_file(self):
        '''read in 16S sequences'''
        
        self.fasta16S = {}

        with open(self.Seq16S_filename) as file_one:

            for line in file_one:

                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    active_sequence_name = line[1:]

                    if active_sequence_name not in self.fasta16S:
                        self.fasta16S[active_sequence_name] = []
                    continue

                sequence = line
                self.fasta16S[active_sequence_name].append(sequence)

    def generate_temp_query_file(self, chassis_org):
        '''generate temporary query file'''

        with open(self._BLAST_QUERY_PATH+'temp_input_query.fa', 'w') as fout:

            fout.write('>'+chassis_org+'\n')
            fout.write(''.join(self.fasta16S[chassis_org])+'\n')

    def generate_blast_dbs(self):
        '''make blast database'''

        if os.path.isfile(self.blastdb+'.nin') is False:
            print("STATUS: Making blastdbs")

            args_makedbs = [_BLAST_PATH+'makeblastdb', '-in', self.Seq16S_filename, '-dbtype',
                            'nucl', '-out', self.blastdb]

            process = subprocess.run(args_makedbs)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # stdoutdata, stderrdata = process.communicate()
        
        else:
            print ('INFO:\t{} blast databases already exist'.format(self.blastdb))

    def blastn_sequences(self, chassis_org):
        '''perform blastn'''

        ##Run blast 
        ### -db = BLAST db (chassis organism)
        ### -query = sequence to search against database 
        ### -task = type of blast to do 
        ### -out = output blast file
        ### -strand = strand to search against (we do both plus and minus)
        args_blastn = [_BLAST_PATH+"blastn", "-db", self.blastdb, 
                       "-query", self._BLAST_QUERY_PATH+'temp_input_query.fa', 
                       "-out", self._BLAST_OUTPUT_PATH+"blast_out_"+chassis_org+'.txt']

        process = subprocess.run(args_blastn)
        # stdoutdata, stderrdata = process.dd()

    def process_blast_output(self,chassis_org):
        '''read in blast output file'''

        self.blast_results = {}
        self.blast_bitscore = {}
        tempset = set()

        with open(self._BLAST_OUTPUT_PATH+"blast_out_"+chassis_org+'.txt') as fin:
            count = 0

            for line in fin:

                line = line.strip()

                if line.startswith('GCA') or line.startswith('GCF'):
                    
                    larray = line.split()
                    if larray[0] not in tempset:
                        count+=1
                        tempset.add(larray[0])
                        self.blast_results[count]=larray[0]
                        self.blast_bitscore[larray[0]]=larray[1]