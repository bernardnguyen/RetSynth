from __future__ import print_function
__author__ = 'Leanne Whitmore, Bernard Nguyen and Corey Hudson'
__email__ = 'leanne382@gmail.com, bernguy@sandia.gov and cmhudso@sandia.gov'
__description__ = 'Main code to RetSynth (RS)'

from multiprocessing import Process
import argparse
try:
    import pickle as cPickle
except:
    import cPickle
import os
import re
import glob
import time
import shutil
import zipfile
from sys import platform
PATH = os.path.dirname(os.path.abspath(__file__))

def unzip_necessary_files_and_libraries():
  '''Unzip Default Databases and indigo libraries'''
  def unzip_folder(foldername):
      if os.path.isdir(foldername) is False:

        zipref = zipfile.ZipFile(foldername+'.zip', 'r')
        zipref.extractall('.')
  try: 
      unzip_folder(PATH+"/ConstructedDatabases")
  except FileNotFoundError:
      print("WARNING: No ConstructedDatabase")

  if platform == 'darwin':
      unzip_folder(PATH+'/indigopython130_mac')

  elif platform == "linux" or platform == "linux2":
      unzip_folder(PATH+'/indigopython130_linux')

  elif platform == "win32" or platform == "win64" or platform == "cygwin":
      unzip_folder(PATH+'/indigopython130_win')

unzip_necessary_files_and_libraries()

from timeit import default_timer as timer
from rsgc.Parser import read_startcompounds as rtsc
from rsgc.Parser import read_targets as rt
from rsgc.Parser import generate_output as go
from rsgc.Parser import structure_similarity as ss
from rsgc.Parser import generate_html as gh
from rsgc.Visualization_chemdraw import reaction_files as rf
from rsgc.ShortestPath import extractinfo as ei
from rsgc.ShortestPath import constraints as co
from rsgc.ShortestPath import integerprogram_pulp as ip_pulp
from rsgc.Database import initialize_database as init_db
from rsgc.Database import build_kbase_db as bkdb
from rsgc.Database import build_modelseed as bms
from rsgc.Database import build_metacyc_db as bmcdb
from rsgc.Database import build_user_rxns_db as burdb
from rsgc.Database import build_ATLAS_db as batlasdb
from rsgc.Database import build_MINE_db as bminedb
from rsgc.Database import build_KEGG_db as bkeggdb
from rsgc.Database import build_SPRESI_db as bspresidb
from rsgc.Database import query as Q
from rsgc.Database import remove_duplicate_cpds as rdc
from rsgc.FBA import build_model as bm
from rsgc.FBA import optimize_target as ot
from rsgc.FBA import compare_results as cr
from rsgc.FBA import retrieve_producable_mets as rpm
from rsgc.FBA import compareKO_results as crko
from rsgc.Visualization_graphviz import SP_Graph_dot as spgd
from rsgc.GeneCompatibility.gc import gc_main as gc
from rsgc.GeneCompatibility.gc import gc_enzyme as gce



def verbose_print(verbose, line):
    if verbose:
        print(line)

def get_compartmentID_from_db(DB, compartment):
    '''Retrieves specified compartment ID'''
    compartment = compartment.lower()
    compartmentID_array = DB.get_compartment(compartment)
    if compartmentID_array is None or len(compartmentID_array) == 0 or compartmentID_array[0] == '':
        print ('WARNING: Could not retrieve a compartment ID from the database')
        if compartment == 'cytosol':
            compartmentID = 'c0'
        elif compartment == 'extracellular':
            compartmentID = 'e0'
        else:
            compartmentID = 'c0'
    else:
        compartmentID = compartmentID_array[0]
    return (compartmentID)

def get_new_temp_imgs_folder(PATH, count):
    '''Check if folder to store images is already present if so new temp folder is made'''
    count+=1
    try:
        os.mkdir(PATH+'/temp_imgs_'+str(count))
        return PATH+'/temp_imgs_'+str(count)
    except OSError:
        PATH_NEW = get_new_temp_imgs_folder(PATH, count)
        return PATH_NEW

def _specific_target(target_id):
    '''Determines if there was a specified organism'''

    if target_id in ['', 'NA', 'N/A']:
        return False

    else:
        return True


def run_flux_balance_analysis(target_info, ex_info, incpds_active,
                              inrxns, media, ko,
                              output, DB, RP, verbose):
    '''
    Run flux balance analysis on target organism with added reactions
    necessary to produce target compound
    '''
    fba = bm.BuildModel(target_info[2], incpds_active, inrxns, DB, verbose, media)
    opt_fba = ot.OptimizeTarget(target_info[0], target_info[2], fba.model, ex_info.temp_rxns,
                                ex_info.temp_exmets, fba.compounds_dict, incpds_active,
                                inrxns, DB, verbose, ko, RP)
    # print (opt_fba.fbasol.fluxes)
    comparisonresults = cr.Compare(target_info[0], fba.solution, opt_fba.fbasol,
                                   ex_info.temp_rxns, DB)
    output.output_FBA(target_info, fba.solution, opt_fba, comparisonresults, ex_info.temp_rxns)
    output.output_theoretical_yield(target_info[0], target_info[2], opt_fba.fbasol,
                                    opt_fba.compounds_dict)

    if ko:

        output.output_essential_reactions(target_info[0], target_info[2], opt_fba.essentialrxns)
        comparisonKOresults = crko.CompareKO(target_info[0], opt_fba.compounds_dict, opt_fba.fbasol,
                                             opt_fba.KOsolutions, ex_info.temp_rxns, DB)
        output.output_FBA_KOs(target_info, opt_fba.fbasol, opt_fba.compounds_dict, comparisonKOresults, ex_info.temp_rxns)

    return opt_fba
 
def retrieve_shortestpath(target_info, IP, LP, LPchem, database, output, temp_imgs_PATH, CRV, timer_output,
                    media_for_FBA, flux_balance_analysis, knockouts, images, figures_graphviz, figures_chemdraw,
                    evaluate_reactions, rankpathways_logP_solvent, rankpathways_logP, rankpathways_boilingpoint,
                    show_rxn_info, output_path, multiple_solutions, start_compounds, gene_compatibility,
                    cai_optimal_threshold, user_cai_table,orgs_gbs, gbs_orgs, RGC, keggorganisms_ids,
                    output_genecompdb, verbose):
    '''Retrieve the shortest path for target organism'''

    start = timer()
    DB = Q.Connector(database)
    SP = None
    ranktype = None

    verbose_print(verbose, "STATUS: getting path for {}".format(target_info))

    if images == False:
        _images = False
    else:
        _images = True

    if not _specific_target(target_info[2]) and not start_compounds:
        print ('WARNING: No organism given therefore target {} compound will be skipped ... '.format(target_info[0]))

    else:

        if start_compounds:

            incpds_active = rtsc.readfile_startcompounds(start_compounds)
            inrxns_active = []

        else:

            incpds_active = DB.get_compounds_in_model(target_info[2])
            inrxns_active = DB.get_reactions_in_model(target_info[2])

        if target_info[0] in incpds_active: #Check if compound exists in organism

            output.output_compound_natively_present_in_target_organism(target_info)

        else:

            optimal_pathways = IP.run_glpk(LP, incpds_active, inrxns_active, target_info[0],
                                        multiplesolutions=multiple_solutions)

            if optimal_pathways:                    
                ex_info = ei.Extract_Information(optimal_pathways, incpds_active, inrxns_active, DB)
                R = rf.ReactionFiles(output_path, DB, ex_info.temp_rxns,
                                    target_info[0], target_info[2], incpds_active,
                                    figures_chemdraw)

                output.output_shortest_paths(target_info, ex_info.temp_rxns)
                output.output_raw_solutions(target_info[0], target_info[2], R.ordered_paths,
                                            ex_info.temp_rxns, ex_info.temp_external, incpds_active, SP)


                if flux_balance_analysis:

                    opt_fba = run_flux_balance_analysis(target_info, ex_info,
                                                        incpds_active, inrxns_active,
                                                        media_for_FBA, knockouts,
                                                        output, DB, SP, verbose)

                    if figures_graphviz:

                        G = spgd.GraphDot(DB, output_path, incpds_active, inrxns_active,
                                        temp_imgs_PATH, opt_fba.fbasol.fluxes)
                        G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)

                    if figures_chemdraw:

                        R.generate_cdxml_files(RP=None, ranktype=None, fba_fluxes=opt_fba.fbasol.fluxes, show_rxn_info=show_rxn_info)

                elif (figures_graphviz or figures_chemdraw) and not flux_balance_analysis:

                    if figures_graphviz:

                        G = spgd.GraphDot(DB, output_path, incpds_active, inrxns_active, temp_imgs_PATH)
                        G.sc_graph(target_info[0], target_info[2], ex_info.temp_rxns, _images)

                    if figures_chemdraw:

                        R.generate_cdxml_files(RP=SP, ranktype=ranktype, show_rxn_info=show_rxn_info)

                if gene_compatibility:
                    verbose_print(args.verbose, 'STATUS:\tOptimizing gene sequences for optimal pathway reaction enzymes...')
                    enzymes = list(R.enzyme_set)
                    if len(enzymes) != 0:
                        output.generate_gc_directory()
                        for enzyme in enzymes:
                            verbose_print(args.verbose,'\nSTATUS:\tOptimizing gene compatibility for EC: %s' % enzyme)
                            if os.path.isfile(os.path.join(output.GC_output_path, "geneseqs_{}_{}.txt".format(enzyme, target_info[2]))):
                                print ("STATUS: already have sequence information for {} in {}".format(enzyme, target_info[2]))
                            else:
                                gce(enzyme, orgs_gbs, gbs_orgs, target_info[2], RGC, keggorganisms_ids, output_genecompdb,
                                    output_directory=output.GC_output_path, user_cai_table=user_cai_table,
                                    cai_optimal_threshold=cai_optimal_threshold)

            else:

                output.output_shortest_paths(target_info, [])


                if flux_balance_analysis:
                    verbose_print(verbose, 'WARNING: No optimal path for %s in species %s therefore no flux balance will be performed' % (target_info[0], target_info[2]))
    end = timer()

    if timer_output:

        output.output_timer('Time to find all paths for {}\t{}\t{}\n'.format(target_info[0], (end-start), (end-start)/60))

    verbose_print(verbose, "Time to find all paths for "+str(target_info[0])+' '+str(end - start))


class RetSynthGC(object):
    """
    Attributes:
        targets = file to list of target compounds
        output_path = directory to output paths
        output_xlsx_format = output results to xlsx file (otherwise output in text file)
        start_compounds = list of compounds in starting material if you do not want to produce compounds in a particular organism
        verbose = output all status updates
        generate_database = generate new database
        generate_database_constraint file = generate necessary processing file for constraints (filename.constraints) necessary when generating new database
        databse = use preconstructed database 
        database_constraint = use preconstructed database file 
        inchidb = use inchi strings for compound IDs in database, takes longer but results in more accurate pathways 
        patric_models = build database with reactions from patric database
        patric_username = to use patric database need to set up an account and put specify username when adding information to a new database 
        patric_password = to use patric database need to set up an account and put specify password when adding information to a new database 
        patric_media = when building models in patric specify media that networks are built off of 
        patricfile = file that tells patric which organisms to include in database (defualt file includes Escherichia coli DH1, Streptomyces venezuelae ATCC 10712 and Corynebacterium argentoratense DSM 44202)
        metacyc = build database with metacyc information 
        metacyc_addition = xml file with metacyc information (Currently we can not process RDF formats of this information)
        kegg = build database with kegg information 
        kegg_orgainsim_type = type of organisms to pull reaction information from 
        kegg_number_of_organisms = how much information to pull from kegg for database
        kegg_number_of_organism_pathways = number of pathways to from kegg for database
        atlas = build database with atlas repository information 
        atlas_dump_directory = directory with raw data atlas 
        mine = build database with MINE repository information 
        mine_dump_directory = directory with raw data from mine
        SPRESI = build database with SPRESI repository information
        spresi_dump_directory = directory wiht raw data from spresi (this has to be purchased not provided with RetSynth)
        user_rxns_2_database = database in a text file that user wants to add to the database
        gene_compatability = run gene compatabiity module (defualt is set to True)
        cai_threshold = amount of codon similarity to be met for gene sequence for enzyme (default = 0.5)
        user_cai_table = path to table of codon usage frequencies 
        flux_balance_analysis = specifies RetSynth to run flux balance analysis 
        media_for_FBA = specifies media to use for FBA, should be the same as patric_media for best results 
        knockouts = specifies RetSynth to perform knockout analysis (takes time, knockouts out each reaction in an organism to and examines if we get more target compound production)
        limit_reactions = limits the number of reactions that can be in novel production pathway (defualt is 10)
        k_number_of_paths = number of sub optimal pathways to identify (defualt is 0, therefore only identifying optimal pathways)
        multiple_solutions = identify multiple optimal suboptimal solutions (defualt is true)
        run_tanimoto_threshold = identify other target compounds similar to your specified target compound as determined by the tanimoto value (defualt is false but if the user wants to run this they must specify a user tanimoto threshold, database must have been generated with --inchidb option)
        figures_chemdraw = generate figures that can viewed in chemdraw only 
        figures_graphviz = generate pathway pngs with graphviz 
        show_rxn_info = shows all reaction information which can be alot if database includes information from SPRESI
    """
    def __init__(self, targets=None, output_path=".", output_html=True, output_xlsx_format=False, processors=4, start_compounds=None, verbose=False,
     generate_database=False, generate_database_constraints=False, database=False, database_constraints=False, inchidb=True, patric_models=False,
     patric_username=None, patric_password=None, patric_reaction_type="bio", patric_media="Complete", patric_sbml_output=False,
     patricfile=PATH+'/Database/data/PATRIC_genome_complete_07152018.csv', patric_models_already_built=False,
     metacyc=False, metacyc_addition=None, metacyc_reaction_type="bio", kegg=False, kegg_reaction_type="bio", kegg_organism_type="bacteria",
     kegg_number_of_organisms="all", kegg_number_of_organism_pathway="all", atlas=False, atlas_dump_directory=None,
     atlas_reaction_type="bio", mine=False, mine_dump_directory=None, mine_reaction_type=True, SPRESI=False,
     spresi_dump_directory=None, spresi_reaction_type="chem", user_rxns_2_database=False, gene_compatability=True, cai_threshold=0.5, user_cai_table=None,
     user_rxns_2_database_type="bio", flux_balance_analysis=False, media_for_FBA="Complete", knockouts=False, limit_reactions=10,
     limit_cycles='None', solver_time_limit=30, evaluate_reactions="all",
     k_number_of_paths=0, multiple_solutions=str(True), cycles=str(True), run_tanimoto_threshold=False, figures_chemdraw=False,
     show_rxn_info=False, figures_graphviz=False, images=True, timer_output=False):
        self.targets = targets
        self.output_path = output_path
        self.output_xlsx_format = output_xlsx_format
        self.processors = processors
        self.start_compounds = start_compounds
        self.verbose = verbose
        self.generate_database = generate_database
        self.generate_database_constraints = generate_database_constraints
        self.database = database
        self.database_constraints = database_constraints
        self.inchidb = inchidb
        self.patric_models = patric_models
        self.patric_username = patric_username
        self.patric_password = patric_password
        self.patric_reaction_type = patric_reaction_type
        self.patric_media = patric_media
        self.patricfile = patricfile
        self.patric_sbml_output = patric_sbml_output
        self.patric_models_already_built = patric_models_already_built
        self.metacyc = metacyc
        self.metacyc_addition = metacyc_addition
        self.metacyc_reaction_type = metacyc_reaction_type
        self.kegg = kegg
        self.kegg_reaction_type = kegg_reaction_type
        self.kegg_organism_type = kegg_organism_type
        self.kegg_number_of_organisms = kegg_number_of_organisms
        self.kegg_number_of_organism_pathway = kegg_number_of_organism_pathway
        self.atlas = atlas
        self.atlas_dump_directory = atlas_dump_directory
        self.atlas_reaction_type = atlas_dump_directory
        self.mine = mine
        self.mine_dump_directory = mine_dump_directory
        self.mine_reaction_type = mine_reaction_type
        self.SPRESI = SPRESI
        self.spresi_dump_directory = spresi_dump_directory
        self.spresi_reaction_type = spresi_reaction_type
        self.user_rxns_2_database = user_rxns_2_database
        self.user_rxns_2_database_type = user_rxns_2_database_type
        self.flux_balance_analysis = flux_balance_analysis
        self.media_for_FBA = media_for_FBA
        self.knockouts = knockouts
        self.limit_reactions = limit_reactions
        self.limit_cycles = limit_cycles
        self.solver_time_limit = solver_time_limit
        self.evaluate_reactions = evaluate_reactions
        self.k_number_of_paths = k_number_of_paths
        self.multiple_solutions = multiple_solutions
        self.cycles = cycles
        self.run_tanimoto_threshold = run_tanimoto_threshold
        self.figures_chemdraw = figures_chemdraw
        self.figures_graphviz = figures_graphviz
        self.images = images
        self.show_rxn_info = show_rxn_info
        self.timer_output = timer_output
        self.rankingpathways_seperation_file = rankingpathways_seperation_file
        self.rankpathways_boilingpoint = rankpathways_boilingpoint
        self.rankpathways_logP = rankpathways_logP
        self.rankpathways_logP_solvent = rankpathways_logP_solvent
        self.gene_compatability = gene_compatability
        self.cai_threshold = cai_threshold
        self.user_cai_table = user_cai_table
        #Generate new output path if need be 
        try:
            verbose_print(self.verbose, "STATUS: generating output folder "+self.output_path)
            os.mkdir(self.output_path)
        except:
            verbose_print(self.verbose, "STATUS: output folder already exists "+self.output_path)
            pass

        #Main program
        self.check_arguments()
        all_db_compounds, all_db_reactions, database = self.retrieve_database_info()
        targets, ignore_reactions, include_rxns, output, temp_imgs_PATH = self.read_in_and_generate_output_files(database)

        LP, LPchem = self.retrieve_constraints(all_db_reactions, all_db_compounds, ignore_reactions, include_rxns, database)
        IP, CRV = self.construct_and_run_integerprogram(targets, output, database, rankingfile=self.rankingpathways_seperation_file)
        if gene_compatibility:
            orgs_gbs, gbs_orgs, R, keggorganisms_ids, output_genecompdb = gc(database, 
                                                                            output_directory=output_path,
                                                                            default_db='%s.db' % DEFAULT_DB_NAME)
        else:
            orgs_gbs=False
            gbs_orgs=False
            R=False
            keggorganisms_ids=False
            output_genecompdb=False

        args_targets = [targets[i:i+self.processors]
                      for i in range(0, len(targets), self.processors)]
        for targets in args_targets:
            processes = []
            for target in targets:
                processes.append(Process(target=retrieve_shortestpath, args=(target, IP, LP, LPchem, database, output,
                                                            temp_imgs_PATH, CRV, self.timer_output, self.media_for_FBA,
                                                            self.flux_balance_analysis, self.knockouts, self.images,
                                                            self.figures_graphviz, self.figures_chemdraw, self.evaluate_reactions,
                                                            self.rankpathways_logP_solvent, self.rankpathways_logP,
                                                            self.rankpathways_boilingpoint, self.show_rxn_info,
                                                            self.output_path, self.multiple_solutions, self.start_compounds,
                                                            self.gene_compatability, self.cai_threshold, self.user_cai_table,
                                                            orgs_gbs, gbs_orgs, R, keggorganisms_ids, output_genecompdb,
                                                            self.verbose)))
                                        
            for p in processes:        
                p.start()
            for p in processes:
                p.join()

        if self.output_xlsx_format:
            output.convert_output_2_xlsx()
        
        if self.output_html:
            print("STATUS: writing html file")
            gh.HtmlOutput(len(targets), self.output_path, self.flux_balance_analysis, self.figures_graphviz, self.output_path+"/Results.html")
        
        '''Remove all temporary images'''
        shutil.rmtree(temp_imgs_PATH)

        '''Removes all dot files if they exist'''
        for filename in glob.glob(self.output_path+"/solution_figures/*.dot"):
            os.remove(filename)
    
    def check_arguments(self):
        '''Checks and makes sure all required arguments are provided'''
        if self.knockouts and not self.flux_balance_analysis:
            raise ValueError('--knockouts option requires that \
                        --flux_balance_analysis option be also specified')
        # if self.database and (not self.generate_database_constraints and not self.database_constraints):
        #     raise ValueErrorprint ('WARNING: User specified specific database but not constraint file therefore default constraint file will be used but may not match user specified database')
    
        if self.metacyc and not self.metacyc_addition:
            raise ValueError('--metacyc requires use of parameters metacyc_addition')        

        if self.atlas and not self.atlas_dump_directory:
            raise ValueError('--atlas requires use of --atlas_dump_directory')

        if self.SPRESI and not self.spresi_dump_directory:
            raise ValueError('--SPRESI requires use of --spresi_dump_directory')

        if self.mine and  not self.mine_dump_directory:
            raise ValueError('--mine requires use of --mine_dump_directory')

        if self.patric_models and not self.patric_password:
            raise ValueError('--patric_models requires options --patric_models, --patric_username, and --patric_password be specified')
        
        if self.patric_models and not self.patric_username:
            raise ValueError('--patric_models requires options --patric_models, --patric_username, and --patric_password be specified')
        
        if self.patric_username and not self.patric_password:
            raise ValueError('--patric_username requires options --patric_models, --patric_username, and --patric_password be specified')

        if not self.patric_models and self.patric_models_already_built:
            raise ValueError('--self.previously_built_patric_models requires options --patric_models, --patric_username, and --patric_password be specified')

        # if not self.patric_models and self.patricfile:
        #    raise ValueError('--self.patricfile requires options --patric_models, --patric_username, and --patric_password be specified')     

        if self.start_compounds and self.flux_balance_analysis:
            raise ValueError('Flux balance cannot be performed on a set of starting compounds\
                        would need to use an organisms metabolism to simulate flux')

        if not self.multiple_solutions and self.k_number_of_paths:
            raise ValueError('Cannot find k_number_of_paths correctly \
                        without finding all multiple_solutions')
        if not self.targets:
            raise ValueError('Requires an input file of target compounds')
        
        if not self.rankingpathways_seperation_file and (self.rankpathways_boilingpoint or self.rankpathways_logP):
            raise ValueError('Requires an boiling point input file for database')

    def retrieve_database_info(self):
        '''
        Generates database or uses previously generated database.
        Can also add metacyc database to a Kbase metabolic database
        '''
        if self.generate_database:
            '''
            Generate a database
            '''
            init_db.Createdb(self.generate_database)
            database = self.generate_database
            new_db = True

        elif self.database:
            new_db = False
            database = self.database

        if self.patric_models:
            #Add PATRIC repository to database
            bms.BuildModelSeed(username=self.patric_username, password=self.patric_password, rxntype=self.patric_reaction_type,
                                inchidb=self.inchidb, DBpath=self.generate_database, output_folder=self.output_path, media=self.patric_media, 
                                patricfile=self.patricfile, newdb=new_db, sbml_output=self.patric_sbml_output,
                                previously_built_patric_models=self.patric_models_already_built)
            
        if self.metacyc:
            #Add metacyc repository to database
            bmcdb.Translate(database, self.metacyc_addition,
                            self.inchidb, self.metacyc_reaction_type, self.verbose)

        if self.kegg and (self.patric_models or self.kbase or self.metacyc):
            #Add kegg repository database
            BKD = bkeggdb.CompileKEGGIntoDB(database, self.kegg_organism_type,
                                            self.inchidb, self.processors, self.kegg_number_of_organisms,
                                            self.kegg_number_of_organism_pathways,
                                            self.kegg_reaction_type, True)

        elif self.kegg and not self.kbase and not self.patric_models and not self.metacyc:
            #Add kegg repository database
            print ('STATUS: Add only KEGG to RetSynth database')
            BKD = bkeggdb.CompileKEGGIntoDB(database, self.kegg_organism_type,
                                            self.inchidb, self.processors,
                                            self.kegg_number_of_organisms, self.kegg_number_of_organism_pathways,
                                            self.kegg_reaction_type, False)

        if self.user_rxns_2_database:
            #Add user identified reactions
            burdb.AddUserRxns2DB(database, self.user_rxns_2_database,
                                model_id='UserAdded', rxntype=self.user_rxns_2_database_type)
        if self.SPRESI:
            #Translate synthetic rdf files from SPRESI into database
            DB = Q.Connector(database)
            cytosol_compartmentID = get_compartmentID_from_db(DB, 'cytosol')
            bspresidb.RDF_Reader(self.spresi_dump_directory,
                                database,
                                self.spresi_reaction_type,
                                cytosol_compartmentID, self.processors)

        if self.mine:
            #Add MINE repository to database
            bminedb.BuildMINEdb(self.mine_dump_directory, database,
                                self.inchidb, self.mine_reaction_type)

        if self.atlas:
            #Add ATLAS repository to database
            batlasdb.build_atlas(self.atlas_dump_directory, database, self.inchidb,
                                    self.processors, self.atlas_reaction_type)
            
        if self.inchidb and (self.patric_models or self.metacyc or self.kegg or self.SPRESI or self.mine or self.atlas):
            #Remove duplicate compounds from database
            rdc.OverlappingCpdIDs(database)

        ##IF DATABASE IS NOT SPECIFIED USE DEFUALT DATABASE IN ./ConstructedDatabases FOLDER##
        if not self.generate_database and not self.database:

            if self.media_for_FBA == 'Carbon-D-Glucose':

                print ('WARNING: No database specified using pre constructed database with media {}'.format(self.media_for_FBA))

                database = PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_GL_MC.db'
                DB = Q.Connector(database)

            elif self.media_for_FBA == 'Complete':

                print ('WARNING: No database specified using pre constructed database with media {}'.format(self.media_for_FBA))

                database = PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_CP_MC_cas_SPRESI_reduced1x.db'            

                DB = Q.Connector(database)
        else:
            DB = Q.Connector(database)

        ###GET TYPE OF REACTIONS TO BE USED IN THE ANALYSIS SO CORRECT CONTRAINT FILE IS GENERATED##
        allcpds = DB.get_all_compounds()
        if self.evaluate_reactions == 'all':
            allrxns = DB.get_all_reactions()
        elif self.evaluate_reactions == 'bio':
            allrxns = DB.get_reactions_based_on_type('bio')
        elif self.evaluate_reactions == 'chem':
            allrxns = DB.get_reactions_based_on_type('chem')
        return(allcpds, allrxns, database)

    def read_in_and_generate_output_files(self, database):
        '''Read in target input file and generate output files'''
        DB = Q.Connector(database)
        R = rt.Readfile(self.targets, DB, self.inchidb)
        if not R.targets:
            raise ValueError('ERROR: No targets, try different compounds')
        temp_imgs_PATH = get_new_temp_imgs_folder(self.output_path, 0)
        OUTPUT = go.Output(DB, self.output_path, self.verbose, self.rankpathways_boilingpoint, self.flux_balance_analysis, self.knockouts, self.timer_output)
        if self.run_tanimoto_threshold:
            verbose_print(self.verbose, 'STATUS: {} tanimoto threshold being used'.format(float(self.tanimoto_threshold)*100))
            cytosol_compartmentID = get_compartmentID_from_db(DB, 'cytosol')
            extracell_compartmentID = get_compartmentID_from_db(DB, 'extracellular')
            SIM = ss.TanimotoStructureSimilarity(R.targets, DB.get_all_compounds(),
                                                cytosol_compartmentID, extracell_compartmentID,
                                                self.verbose, self.tanimoto_threshold)
            OUTPUT.output_final_targets(SIM.finaltargets, self.tanimoto_threshold)
            return(SIM.finaltargets, R.ignorerxns, R.includerxns, OUTPUT, temp_imgs_PATH)

        else:
            return(R.targets, R.ignorerxns, R.includerxns, OUTPUT, temp_imgs_PATH)

    def retrieve_constraints(self, allrxns, allcpds, ignore_reactions, include_rxns, database):
        '''
        Generates database constraints or uses previously generated
        database constraints (.constraints) file
        '''
        DB = Q.Connector(database)

        def store_constraint_file(filename, LP):
            '''store generated constraints into .constraints'''
            
            print ('STATUS: Dumping newly generated variables...')

            with open(filename, 'wb') as fout1:
                
                print ('STATUS: Dumping LP structure...')
                cPickle.dump(LP.lp, fout1)
                
                print ('STATUS: Dumping all compounds...')
                cPickle.dump(LP.allcpds, fout1)

                print ('STATUS: Dumping all reaction variables...')
                cPickle.dump(LP.variables, fout1)

                print ('STATUS: Dumping all reactions 1...')
                cPickle.dump(LP.allrxnsrev_dict_rev, fout1)
    
                print ('STATUS: Dumping all reactions 2...')
                cPickle.dump(LP.allrxnsrev_dict, fout1)

                print ('STATUS: Dumping all reactions 3...')
                cPickle.dump(LP.allrxnsrev, fout1)

        def unload_constraint_file(filename):
            '''unload constraints from a .constraints'''

            print ('STATUS: Loading pre-stored variables...')

            with open(filename, 'rb') as fin1:
                
                print ('STATUS: Loading LP structure...')
                lp = cPickle.load(fin1)

                print ('STATUS: Loading all compounds...')
                allcompounds4matrix = cPickle.load(fin1)

                print ('STATUS: Loading all reaction variables...')
                variables = cPickle.load(fin1)

                print ('STATUS: Loading all reactions 1...')
                allrxnsrev_dict_rev = cPickle.load(fin1)

                print ('STATUS: Loading all reactions 2...')
                allrxnsrev_dict = cPickle.load(fin1)

                print ('STATUS: Loading all reactions 3...')
                allrxnsrev = cPickle.load(fin1)

            return (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev)

        def load_preconstructed_constraint_files(media, mediatype, args):
            '''load defualt .constraint files'''

            print ('WARNING: No database constraint file specified using pre constructed database constraint file for database with media {}'.format(self.media_for_FBA))
            if self.evaluate_reactions == 'all':

                (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = unload_constraint_file(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC_cas_SPRESI_reduced1x.constraints'.format(mediatype))
                
                return (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev)
            
            elif  self.evaluate_reactions == 'bio':
                
                (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = unload_constraint_file(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC_cas_SPRESI_reduced1x_bio.constraints'.format(mediatype))
                
                return (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev)            
            
            elif self.evaluate_reactions =='chem':

                (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = unload_constraint_file(PATH+'/ConstructedDatabases/DBINCHIECOLIDH1_{}_MC_cas_SPRESI_reduced1x_chem.constraints'.format(mediatype))

                
                return (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev)  



        ###RETRIEVE SPECIFIED BY USER CONSTRAINTS###
        if self.generate_database_constraints:
            LP = co.ConstructInitialLP(allrxns, allcpds, DB, ignore_reactions, include_rxns)
            store_constraint_file(self.generate_database_constraints, LP)

        elif self.database_constraints:
            (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = unload_constraint_file(self.database_constraints) 
            LP = co.ConstructInitialLP(allrxns, allcompounds4matrix, DB,
                                    ignore_reactions, include_rxns,
                                    lp, variables, allrxnsrev_dict_rev,
                                    allrxnsrev_dict, allrxnsrev)
        else:
            if self.media_for_FBA=='Carbon-D-Glucose':
                
                (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = load_preconstructed_constraint_files(self.media_for_FBA, 'GL', args)
            
            elif self.media_for_FBA=='Complete':
                
                (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = load_preconstructed_constraint_files(self.media_for_FBA, 'CP', args)
            
                LP = co.ConstructInitialLP(allrxns, allcompounds4matrix, DB,
                                        ignore_reactions, include_rxns,
                                        lp, variables, allrxnsrev_dict_rev,
                                        allrxnsrev_dict, allrxnsrev)
            else: 
                print ('ERROR: No identified pre constraint file...stopping run')

        LPchem=None

        if self.rankpathways_boilingpoint or self.rankpathways_logP:
            print ('STATUS: Load necessary chemical constraints for ranking pathways based on boiling point')

            constraint_file = re.sub('.db$', '_chem.constraints', database)
            allrxns_chem = DB.get_reactions_based_on_type('chem')

            if not  os.path.isfile(constraint_file):
                LPchem = co.ConstructInitialLP(allrxns_chem, allcpds, DB, [], [])
                store_constraint_file(constraint_file, LPchem)         

            else:
                (lp, allcompounds4matrix, variables, allrxnsrev_dict_rev, allrxnsrev_dict, allrxnsrev) = unload_constraint_file(constraint_file)  

                LPchem = co.ConstructInitialLP(allrxns, allcompounds4matrix, DB,
                                                [], [], lp, variables, allrxnsrev_dict_rev,
                                                allrxnsrev_dict, allrxnsrev)

        return LP, LPchem

    def construct_and_run_integerprogram(self, targets, output, database, rankingfile=None):
        '''
        Constructs ILP and solves it identifying shortest path to the target
        '''
        
        DB = Q.Connector(database)
        
        if self.timer_output:

            IP = ip_pulp.IntergerProgram(DB, self.limit_reactions,
                                        self.limit_cycles, self.k_number_of_paths,
                                        self.cycles, self.verbose, self.solver_time_limit, output)
        else:

            IP = ip_pulp.IntergerProgram(DB, self.limit_reactions,
                                        self.limit_cycles, self.k_number_of_paths,
                                        self.cycles, self.verbose, self.solver_time_limit, self.timer_output)
        if self.rankpathways_boilingpoint or self.rankpathways_logP:

            CRV = rp.ConstructRankingVariables(database, rankingfile, self.verbose)

            return (IP, CRV)
        
        else: 
        
            return (IP, None)