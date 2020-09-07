import urllib.request, urllib.error, urllib.parse, http.client
from tqdm import tqdm
import re
import os
import pickle
from rsgc.GeneCompatibility.KEGG import kegg_functions as kf
from rsgc.GeneCompatibility.GenBank import genbank_functions as gbf

codontable = {  
         'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',  
         'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  
         'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',  
         'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                   
         'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',  
         'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',  
         'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',  
         'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',  
         'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',  
         'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',  
         'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',  
         'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',  
         'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',  
         'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',  
         'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',  
         'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',  
    } 
stop_codons  = ["TAA", "TAG", "TGA"]

class GeneSeqKEGG(object):
	def __init__(self, orgs2pullgene=False, get_orgs=False):
		'''initialize'''

		self.orgs2pullgene = orgs2pullgene
		self.get_orgs = get_orgs
		self.genesequences = {}
		self.GenBankIDs = {}

	def get_geneseq_for(self, enzyme):
		'''Get genes for target enzyme for all bacterial organisms'''
		# getting_sequences = 'STATUS:\tgetting gene sequences for EC: %s' % enzyme
		# if self.orgs2pullgene:
		# 	getting_sequences += '\nSTATUS:\t\tfor KEGG organisms %s' % ','.join(self.orgs2pullgene)		
		# print (getting_sequences)
		
		if self.get_orgs:
			self.orgs_with_enzyme = list()

		if enzyme in self.genesequences:
			return self.genesequences[enzyme]
			
		geneseq = {}
		darray = kf.extract_KEGG_data('get/ec:'+enzyme)
		if darray:	
			idxstart = 0
			try:
				while not darray[idxstart].startswith('GENES'):
					idxstart += 1
				idxend = idxstart + 1
				while darray[idxend].startswith(' '):
					idxend += 1
				indent = len(re.search('GENES\s+', darray[idxstart]).group(0))
				darray = [d[indent:] for d in darray[idxstart:idxend]]

				for line in darray:
					# line looks like
					# ORGID: geneID(...) geneID geneID(...)
					data = line.lower().split(' ')
					# data looks like
					# [orgID:, geneID(...), geneID, geneID(...)]
					orgID = data[0][:-1]
					if self.get_orgs:
						self.orgs_with_enzyme.append(orgID)
				
					elif self.get_orgs is False and orgID in self.orgs2pullgene:
						genes = data[1:]
						geneseq[orgID] = self.extract_KEGG_geneseq(genes, orgID)
				self.genesequences[enzyme] = geneseq

			except IndexError:
				print ('WARNING:\t{} in KEGG has no gene information'.format(enzyme))				
				print('STATUS:\tgetting sequence information for EC:%s...' % enzyme)
				if self.get_orgs is False:
					geneseq = self.get_seq_for(enzyme, darray)
					self.genesequences[enzyme] = geneseq
		else:
			print(('STATUS:\tFailed to get KEGG data for EC:%s' % enzyme))
		return geneseq

	def extract_KEGG_geneseq(self, genes, orgID):
		'''Pull gene sequence for target enzyme from KEGG'''
		geneseq = {}
		for gene in genes:
			gene = gene.split('(')[0]
			darray = kf.extract_KEGG_data('get/%s:%s/ntseq' % (orgID, gene))
			if darray:
				geneseq[gene] = ''.join(darray[1:])
			else:
				print(('WARNING:\tFailed to get gene sequence %s for organism %s' % (gene, orgID)))
		return geneseq

	def get_seq_for(self, enzyme, darray=None):
		'''Pull gene sequence for target enzyme from KEGG from all bacterial organisms'''
		geneseq = {}
		if not darray:
			darray = kf.extract_KEGG_data('get/ec:'+enzyme)

		if darray:
			idxstart = 0
			try:
				while not darray[idxstart].startswith('  SEQUENCE'):
					idxstart += 1
				idxend = idxstart + 1
				while darray[idxend].startswith(' '):
					idxend += 1
				indent = len(re.search('\s+SEQUENCE\s+', darray[idxstart]).group(0))
				darray = [d[indent:] for d in darray[idxstart:idxend]]

				for line in darray:
					# line looks like
					# [orgID:geneID]
					data = line.lower()[1:-1]
					# data looks like
					# orgID:geneID
					if data:
						print('STATUS:\tfound sequence information for EC:%s' % enzyme)
						[orgID, geneID] = data.split(':')
						print('STATUS:\tgetting GenBank ID for %s:%s' % (orgID, geneID))
						GenBankID = self.get_GenBank_id(enzyme, data)
						genbank_nt_seq, genbank_pt_seq = gbf.extract_GenBank_sequence(GenBankID)

						if genbank_nt_seq != '' and genbank_pt_seq !='':
							print ('INFO:\tpulled nucleotide and protein sequence for genbank ID {} going to pull coding sequence for nucleotide sequence'.format(GenBankID))

							cds = self.extract_cds(genbank_nt_seq, genbank_pt_seq)
							geneseq[orgID] = {geneID:cds}

						else: 
							print ('WARNING:\tdo not have gene and protein sequence so unable to accurately pull out cds')

							if GenBankID and genbank_nt_seq != '':
								geneseq[orgID] = {geneID:genbank_nt_seq}

							else:
								geneseq[orgID] = {geneID:'no_gene_info'}

			except IndexError:
				print ('WARNING:\t{} in KEGG has no sequence information'.format(enzyme))
		else:
			print(('WARNING:\tFailed to get KEGG data for EC:%s' % enzyme))
		return geneseq

	def get_GenBank_id(self, enzyme, orgIDgeneID):
		'''Pull GenBank ID for sequence associated with target enzyme'''
		darray = kf.extract_KEGG_data('get/'+orgIDgeneID)

		gbid = None
		if darray:	
			idxstart = 0
			try:
				while not darray[idxstart].startswith('DBLINKS'):
					idxstart += 1
				idxend = idxstart + 1
				while darray[idxend].startswith(' '):
					idxend += 1
				indent = len(re.search('DBLINKS\s+', darray[idxstart]).group(0))
				darray = [d[indent:] for d in darray[idxstart:idxend]]

				for line in darray:
					# line looks like
					# DATABASE: ID
					data = line.split(': ')
					# data looks like
					# [DATABASE, ID]
					if data[0] == 'GB':
						gbid = data[1]
						self.GenBankIDs[enzyme] = gbid
						break
			except IndexError:
				print ('WARNING:\t{} in KEGG has no GenBank information'.format(enzyme))
		else:
			print(('WARNING:\tFailed to get KEGG data for %s' % orgIDgeneID))
		return gbid


	def extract_cds(self, genseq, protseq):
		'''extract coding seqeuence from mrna sequence'''

		cds = ''
		derived_protseq = ''
		
		start_codon_loc = genseq.find('atg')
		new_genseq = genseq[int(start_codon_loc):]
		new_genseq_array = [new_genseq[i:i+3] for i in range(0, len(new_genseq), 3)]

		for codon in new_genseq_array:
			codonUP = codon.upper()

			if codonUP not in stop_codons:
				derived_protseq=derived_protseq+codontable[codonUP]
				cds = cds+codon

			elif codonUP in stop_codons:
				cds = cds+codon
				break
		
		if derived_protseq == re.sub('"', '', protseq):
			return cds
		
		else:
			print ('WARNING:\tderived protein sequence does not match actual protein sequence so have not identified correct ORF attempting again ')
			new_genseq = re.sub('^atg', '', new_genseq)
			self.extract_cds(new_genseq, protseq)

