__author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'lwhitmo@sandia.gov, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'main run genecompatability'

import os
import re
import pickle
import glob
from copy import deepcopy
from rsgc.GeneCompatibility.retrievegeneseqs import genbank_4_orgs as g4o
from rsgc.GeneCompatibility.retrievegeneseqs import geneSeqs_KEGG as gse
from rsgc.GeneCompatibility.BLAST import run_blast as rb
from rsgc.GeneCompatibility.NCBI_SSU import pull_16S_seq_NCBI as psN
from rsgc.GeneCompatibility.KEGG import kegg_functions as kf
try:
    from bio.Data import IUPACData
except ImportError:
    from Bio.Data import IUPACData
from rsgc.GeneCompatibility.D_Tailor import SequenceAnalyzer as sa
from rsgc.GeneCompatibility.D_Tailor import Solution as so 
from rsgc.GeneCompatibility.D_Tailor.Features import CAI 
# from GeneCompatibility.D_Tailor.runDTailor import CAIEcoliDesigner
from rsgc.GeneCompatibility.D_Tailor.runDTailor import CAIDesigner
from rsgc.GeneCompatibility.D_Tailor.DesignOfExperiments.Design import RandomSampling,Optimization,FullFactorial
from rsgc.GeneCompatibility.D_Tailor.RunningExamples.Designer.TranslationFeaturesEcoliDesigner import TranslationFeaturesEcoliDesigner

###LOAD IUPAC DICTIONARY### 
iupacdict = IUPACData.ambiguous_dna_values
del iupacdict['A']
del iupacdict['C']
del iupacdict['T']
del iupacdict['G']

PATH = os.path.dirname(os.path.abspath(__file__))

def verbose_print(verbose, line):
    if verbose:
        print(line)

def make_tmp_folders():    
    def make_folder(fp):
        try:
            os.mkdir(fp)
        except FileNotFoundError:
            make_folder('/'.join(fp.split('/')[:-1]))
            make_folder(fp)
        except FileExistsError:
            pass

    folders = [PATH+'/D_Tailor/tmp/transterm_files',
                PATH+'/D_Tailor/tmp/structures',
                PATH+'/D_Tailor/tmp/unafold_files',
                PATH+'/.temp_fasta_files',
                PATH+'/BLAST/blastoutput',
                PATH+'/BLAST/query',
                PATH+'/BLAST/blastdbs',
                PATH+'/NCBI_SSU/ncbi_gn_data']
    for folder in folders:
        make_folder(folder)



def empty_tmp_folders():
    '''empty temporary folders'''

    def clear_folder(files):
        for file in files: 
            os.remove(file)

    folders =[PATH+'/D_Tailor/tmp/transterm_files/*.*', PATH+'/.temp_fasta_files/*.fa',
              PATH+'/D_Tailor/tmp/structures/*.*',PATH+'/D_Tailor/tmp/unafold_files/*.*',
              PATH+'/BLAST/blastoutput/*.*', PATH+'/BLAST/query/*.fa']

    for folder in folders:
        files = glob.glob(folder)
        clear_folder(files) 

def reverse_org_gbs_dict(orgs_gbs):
    '''generates new dictionary where the genebank ids 
    are the keys and values are database and kegg ids'''

    gbs_orgs = {}

    for key, gbs in orgs_gbs.items():
        for gb in gbs:
            gbs_orgs.setdefault(gb, []).append(key)

    return gbs_orgs


def gc_main(database,  output_directory='', default_db=''):
    ###GETS GENEBANK IDS FOR ORGANISMS IN OUR DATABASE COLLECTED FROM REPOSITORIES PATRIC AND KEGG###

    try:
        output_genecompdb = os.path.abspath(os.path.join(output_directory,'genecomp_databases'))
        os.mkdir(output_genecompdb)
    except:
        output_genecompdb = os.path.abspath(os.path.join(output_directory,'genecomp_databases'))
        pass
    db_org_gbs = g4o.get_database_organism_genbank_ids(database)

    DB_NAME = database.split('/')[-1]
    if default_db and DB_NAME == default_db:
        USING_DEFAULTS = True
        DB_NAME = 'defaultdb'
    else:
        USING_DEFAULTS = False
        

    ###GETS ADDITIONAL ORGANISMS (PLANT AND FUNGI) NOT IN OUR DATABASE###
    if not USING_DEFAULTS and os.path.isfile(os.path.join(output_genecompdb,'kegg_bac_plant_fungi_%s.list' % DB_NAME)) is False:
        print ('STATUS:\tDo not have pre-saved list of kegg bacteria, plant and fungi species therefore generating new file...')
        
        orgs_gbs, keggorganisms = kf.extract_KEGG_orgIDs(db_org_gbs, taxid_info=True, kegg_id=True)
        
        with open(os.path.join(output_genecompdb, 'kegg_bac_plant_fungi_%s.list' % DB_NAME), 'wb') as fout:
             pickle.dump(keggorganisms, fout)

        ##GET REVERSE DICTIONARY OF orgs_gbs###        
        gbs_orgs = reverse_org_gbs_dict(orgs_gbs)

    else:
        print ('STATUS:\tHave pre-saved list of kegg bacteria, plant and fungi species therefore opening file...')

        orgs_gbs = db_org_gbs
        print(output_genecompdb)
        with open(os.path.join(output_genecompdb, 'kegg_bac_plant_fungi_%s.list' % DB_NAME), 'rb') as fin:
             keggorganisms = pickle.load(fin)

    ###ADDING KEGG ORGANISMS TO LIST OF ORGANISMS IN OUR DATABASE###
    keggorganisms_ids = {}

    orgs_gbs_bac = deepcopy(orgs_gbs)
    for kegg, value_dict in keggorganisms.items():

        for key, value in value_dict.items():

            keggorganisms_ids[key] = value_dict['org_type']

            if key not in orgs_gbs and key != 'org_type':

                if value_dict['org_type'] == 'bacteria':
                    orgs_gbs_bac[key] = value

                orgs_gbs[key] = value

    ##GET REVERSE DICTIONARY OF orgs_gbs###
    gbs_orgs = reverse_org_gbs_dict(orgs_gbs)

    ###GET 16S SEQUENCE TO CALCULATE EVOLUTIONARY DISTANCES FROM CHASSIS ORGANISM###
    ###NOTE: GENOME RNA SEQUENCES ARE STORED IN THE FOLDER NCBI_SSU/ncbi_gn_data/ WHEN
    ###COMPLETE THIS FOLDER IS ABOUT 5GBS IN SIZE FOR DEFAULT DB, THE wipe_folder OPTION IF SET TO 
    ###True WILL REMOVE ALL SEQUENCES WHEN DONE.
    print('NCBI')
    psN.NCBI_SSU(orgs_gbs_bac, os.path.join(output_genecompdb, 'kegg_bac_16S_%s.fa' % DB_NAME), output_genecompdb, wipe_folder=False)
 
    ###CALCULATE EVOLUTIONARY DISTANCES AND PULL DISTANCES FOR CHASSIS ORGANISM, NOTE DISTANCES ARE SAVED IN A DISTMAT FILE IN phylo/data FOLDER 
    ###IF THIS FILE IS NOT ALREADY GENERATED THIS WILL TAKE SOME TIME TO BUILD###
    try: 
        blast_path = os.path.join(output_genecompdb, 'BLAST')
        os.mkdir(blast_path)
    except:
        pass
    try:
        blastdb_path = os.path.join(blast_path, 'blastdbs')
        os.mkdir(blastdb_path)
    except:
        pass

    R = rb.ImplementBLAST(os.path.abspath(os.path.join(blastdb_path,'kegg_bac_16S_%s' % DB_NAME)),
                          os.path.abspath(os.path.join(output_genecompdb, 'kegg_bac_16S_%s.fa' % DB_NAME)), os.path.abspath(blast_path))
    R.generate_blast_dbs()
    return (orgs_gbs, gbs_orgs, R, keggorganisms_ids, output_genecompdb)

def gc_enzyme(enzyme, orgs_gbs, gbs_orgs, target_org, R, keggorganisms_ids, output_genecompdb, output_directory, user_cai_table=None, cai_optimal_threshold=0.50):
    ###COLLECT TYPE OF ORGANISMS (BACTERIA, PLANT OR FUNGI) THAT HAVE THE GENE###
    R.read_16S_fasta_file()
    R.generate_temp_query_file(orgs_gbs[target_org][0])
    R.blastn_sequences(orgs_gbs[target_org][0])
    R.process_blast_output(orgs_gbs[target_org][0])

    GE = gse.GeneSeqKEGG(orgs2pullgene=False, get_orgs=True)
    GE.get_geneseq_for(enzyme)

    org_types = []

    for k in GE.orgs_with_enzyme:
        try:
            org_types.append(keggorganisms_ids[k])
        except KeyError:
            pass

    ###TYPE OF ORGANISMS FOUND TO HAVE ENZYME###
    org_types = set(org_types)

    ###NUMBER OF ORGANISMS TO PULL GENE OUT OF###
    top_num = 10

    ###IF BACTERIA HAVE THE ENZYME OF INTEREST PULL OUT SEQUENCES FROM CLOSELY RELATED ORGANISMS
    if 'bacteria' in org_types:
        # print ('STATUS:\tbacterial species have this enzyme so running evolutionary distances to get most evolutionary close organisms to chassis')

        ###CHECK TO MAKE SURE GENES WERE AQUIRED###
        num_orgs = 0

        ###RANGE TO GET NUMBER OF TOP ORGANISMS###
        start = 1
        end = start + top_num
        # print(gbs_orgs)
        all_blast_results = [gbs_orgs[R.blast_results[i]][0] for i in R.blast_results.keys()]
        intersect = set(all_blast_results).intersection(set(GE.orgs_with_enzyme))

        if intersect:
            while num_orgs == 0:

                # print('\nSTATUS:\tattempting to retrieve gene sequences for closely related organisms (ranked %d-%d)' % (start, end-1))

                top_distances = [i for i in range(start, end)]
                top_distance_orgs = [gbs_orgs[R.blast_results[i]][0] for i in top_distances]

                GE = gse.GeneSeqKEGG(orgs2pullgene=top_distance_orgs, get_orgs=False)
                GE.get_geneseq_for(enzyme)

                if len(GE.genesequences[enzyme]) > 0:
                    num_orgs = len(GE.genesequences[enzyme])
                    # print('STATUS:\tfound gene sequences for closely related organisms')

                else:
                    if end == len(R.blast_results):
                        # print('STATUS:\tno gene info found for any closely related organisms')
                        break

                    # print ('STATUS:\tno gene info found for orgs {}-{}'.format(start,end-1))                

                    start += top_num
                    end += top_num

                    if end > len(R.blast_results):
                        end = len(R.blast_results)
        else:
            print ('STATUS: No overlap between blast results and bacterial organisms with enzyme subsequently just pull 10 gene sequences')    
            top_organisms = GE.orgs_with_enzyme[:10]
            GE = gse.GeneSeqKEGG(orgs2pullgene=top_organisms, get_orgs=False)
            GE.get_geneseq_for(enzyme)


    else:
        print ('STATUS:\tno bacterial species have enzyme {} therefore pulling out {} gene sequences for enzyme from fungi or plant species that have it'.format(enzyme,top_num))
        orgs_with_gene = GE.orgs_with_enzyme[:top_num]
        GE = gse.GeneSeqKEGG(orgs2pullgene=orgs_with_gene, get_orgs=False)
        GE.get_geneseq_for(enzyme)

    if enzyme in GE.genesequences and len(GE.genesequences[enzyme]) > 0:
        try:
            os.mkdir(os.path.join(output_genecompdb,'.temp_fasta_files'))
        except:
            pass
        with open(os.path.join(output_genecompdb,'.temp_fasta_files/tempseqfile._{}.fa'.format(enzyme)), 'w') as fout:
            for key, value_dict in GE.genesequences.items():
                for key2, value_dict2 in value_dict.items():
                    for key3 in value_dict2:
                        fout.write('>'+key2+':'+key3+'\n')
                        tempdna = []
                        for s in value_dict2[key3]:
                            # print (s)
                            S = s.upper()
                            if S in iupacdict:
                                print ('WARNING:\tiupac base in dict')
                                s = iupacdict[S][0].lower()
                                tempdna.append(s)
                            else:
                                tempdna.append(s)           
                            tempdnastring = ''.join(tempdna)
                            
                        fout.write(tempdnastring+'\n')

        #D_Tailor
        SA = sa.SequenceAnalyzer(os.path.join(output_genecompdb,'.temp_fasta_files/tempseqfile._{}.fa'.format(enzyme)), 'FASTA')
        seqs = SA.list_of_input_sequences
        cai_scores = []
        design_param = {"cdsCAI" : { 'type' : 'REAL' , 
                                    'thresholds' : { '1': (0.13,0.29), 
                                                        '2': (0.29,0.33), 
                                                        '3': (0.33,0.37),
                                                        '4': (0.37,cai_optimal_threshold), 
                                                        '5': (cai_optimal_threshold,0.86)}}
                        }

        design = Optimization(["cdsCAI"],design_param, '5')  

        org_cai_table = get_cai_table(target_org, user_cai_table)
        max_cai_index, max_cai, cai_scores = get_max_cai(seqs, cai_scores, design, org_cai_table, initial=True)

        outputfile = open('%s/geneseqs_%s_%s.txt' % (output_directory,enzyme,target_org), 'w')


        if max_cai < cai_optimal_threshold:
            dtailorbool = run_DTailor(seqs, max_cai_index, max_cai, cai_scores, design, org_cai_table, outputfile)

            if dtailorbool is None:
                outputfile.write('>No Sequence found \n')

        else:
            outputfile.write('>original_sequence (best) CAI: {}\n'.format(max_cai))
            outputfile.write(seqs[max_cai_index]['sequence']+'\n')
        
        outputfile.close()
    else:
        print('STATUS:\tFailed to retrieve gene sequences for EC:%s' % enzyme)
    empty_tmp_folders()

def get_max_cai(seqs, cai_scores, design, cai_table, initial=True):
    '''get max cai from list of sequences'''

    if initial:
        for value in seqs:
            solution = so.Solution(sol_id=value['name'], sequence=value['sequence'],design=design, cai_table=cai_table)
            cai_obj = CAI.CAI(solution=solution, label="cds", cai_table=cai_table,
                                args= {'cai_range' : (0,len(solution.sequence)),
                                        'mutable_region' : list(range(3,len(solution.sequence)))})
            cai_scores.append(cai_obj.scores['cdsCAI'])

        cai_scores, max_cai, max_cai_index=new_max_cai(cai_scores, '')

    return max_cai_index, max_cai, cai_scores

def new_max_cai(cai_scores, max_cai_index, initial=True):
    '''get sequence with max cai score'''
    
    if initial is False:
        del cai_scores[max_cai_index]
        if len(cai_scores) == 0:
            return [],0,''

    max_cai = max(cai_scores)
    max_cai_index = cai_scores.index(max_cai)    

    return cai_scores, max_cai, max_cai_index

def run_DTailor(seqs, max_cai_index, max_cai, cai_scores, design, org_cai_table, outputfile):
    '''Run D-Tailor analysis'''
    print ('STATUS:\tCAI is {} so implementing D-Tailor Sequence Designer to identify sequence with better CAI'.format(max_cai))

    # tfec_designer = CAIEcoliDesigner(seqs[max_cai_index]['name'], seqs[max_cai_index]['sequence'],
    #                                     design, PATH+'/genecompatability', outputfile, createDB=True)
    tfec_designer = CAIDesigner(seqs[max_cai_index]['name'], seqs[max_cai_index]['sequence'],
                                        design, PATH+'/genecompatability', outputfile, org_cai_table, createDB=True)

    solution = so.Solution(seqs[max_cai_index]['name'], seqs[max_cai_index]['sequence'],design=design,cai_table=org_cai_table)
    tfec_designer.configureSolution(solution)
    valid = tfec_designer.validateSolution(solution)

    if valid:
        outputfile.write('>original_sequence CAI: {}\n'.format(max_cai))
        outputfile.write(seqs[max_cai_index]['sequence']+'\n')
        res = tfec_designer.run(selection="directional")
        return True
    
    else:
        cai_scores, max_cai, max_cai_index = new_max_cai(cai_scores, max_cai_index, initial=False)
        if len(cai_scores) > 0: 
           dtailorbool =  run_DTailor(seqs, max_cai_index, max_cai, cai_scores, design, org_cai_table, outputfile)
           return dtailorbool
        else:
            print ('WARNING:\tNone of the sequences initially pulled are valid therefore no better sequence found')
            return None

def read_cai_table(cai_table_path):
    codon_table=dict()
    with open(cai_table_path) as fin:
        
        lines = [line for line in fin.read().split('\n')]
        lines = [re.sub("\(\s*\d+\)", "", line) for line in lines]
        linesx = [line.split() for line in lines]
        for line in linesx:
            if len(line) > 0:
                it = iter(line)
                for x in it: 
                    codon_table.setdefault(x.lower(), round(float(next(it))*0.01, 2))
    return codon_table

def get_cai_table(organism, user_cai_table):

    try:
        cai_table = read_cai_table(PATH+'/D_Tailor/CAI_Tables/%s_cai.txt' % organism)
    except:
        try:
            cai_table = read_cai_table(user_cai_table)
        except:
            print('ERROR:\tNo CAI table found for organism %s.' % organism)
    return cai_table

if __name__ == '__main__':

    enzyme = '1.1.1.80'
    # enzyme = '4.2.3.110'
    database = '/home/leann/sandia_work/DevelopedDatabases/PATRICrhodoONLY_mc_inchi.db'
    # database = '/Users/bernguy/Documents/RetSynth_lt/rs/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC_cas_SPRESI_reduced1x.db'
    # default_db = 'DBINCHIECOLIDH1_CP_MC_cas_SPRESI_reduced1x.db'
    target_org = '269796.9'
    # target_org = '953739.5'
    # target_org = '536056.3'
    # target_org = '1348662.3'

    orgs_gbs, gbs_orgs, R, keggorganisms_ids, output_genecompdb = gc_main(database, output_directory="../../../testgc/")
    gc_enzyme(enzyme, orgs_gbs, gbs_orgs, target_org, R, keggorganisms_ids, output_genecompdb, output_directory="../../../testgc/", user_cai_table=None, cai_optimal_threshold=0.30)