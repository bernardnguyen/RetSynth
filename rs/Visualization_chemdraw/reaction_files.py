__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Generate Reaction smiles file (can be used in ChemDraw)'
import re
import os
from sys import platform
import pubchempy as pcp
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_renderer
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_renderer
elif platform == "win32" or platform == 'win64' or platform=="cygwin":
    from indigopython130_win import indigo
    from indigopython130_win import indigo_renderer
from Visualization_chemdraw.exclude_cpds import retrieve_promiscuous_mets
from Visualization_chemdraw.cdxml_editor import CDXML_Editor


''' Tree Structure for Ordered Path '''
class Node(object):
    def __init__(self, root=None):
        self.root = root
        self.children = []

    def addChild(self, child):
        self.children += [child]

class Tree(object):
    def __init__(self, root=None):
        self.root = root

    def set_root(self, new_root):
        self.root = new_root

    def list_nodes(self):
        ''' Depth first search, outputs list of all nodes'''
        def traverse(root):
            output = []
            for child in root.children:
                output += traverse(child)
            output += [root.root]
            return output

        return traverse(self.root)


class ReactionFiles(object):
    def __init__(self, output_path, DB, reactions, target,
                 target_organism, incpds, cdxmlfiles=True):
        self.output_path = output_path
        self.DB = DB
        self.incpds = incpds
        self.incpds_set = set(incpds)
        self.reactions = reactions
        self.target = target
        self.target_organism = target_organism
        self.target_organism_name = self.DB.get_organism_name(target_organism)
        self.cdxmlfiles = cdxmlfiles
        self.ordered_paths = self.order_of_paths()
        self.promiscuous = retrieve_promiscuous_mets(DB)
        self.enzyme_set = set()

    def get_cdxml(self, compounds, cdxml_cpds, promiscuous_cpds):
        '''Get cdxml files for relevant compounds'''
        IN = indigo.Indigo()
        IR = indigo_renderer.IndigoRenderer(IN)

        for cpd in compounds:
            cpd_name = compounds[cpd]           
            if cpd.startswith('InChI'):                
                cpd = re.sub('_\w+\d+$', '', cpd)
                if cpd_name == 'None':
                    cpd_name = '_'.join(cpd.split('/'))

                mol = IN.loadMolecule(cpd)                              
                if self.promiscuous.get(cpd_name) is None:
                    cdxml_cpds.append(cpd_name)
                    IN.setOption("render-comment", cpd_name)
                    IN.setOption("render-output-format", "cdxml")
                    IR.renderToFile(mol, self.output_path+'/compounds/'+cpd_name+'.cdxml')
                else:
                    promiscuous_cpds.append(cpd_name)
            else:
                if cpd_name == 'None':
                    cpd_name = cpd
                cdxml_cpds.append(cpd_name)
        return cdxml_cpds, promiscuous_cpds

    def generate_output_folders(self, folder):
        '''Check to see if output folders are present
            if not set up new output folder'''
        try:
            os.mkdir(folder)
        except OSError:
            pass

    def order_of_paths(self):
        ''' Ordered pathway with target as tree root '''
        ordered_paths = {}

        for count_pathway, os_dict in list(self.reactions.items()):
            ordered_paths[count_pathway] = Tree()
            rxn_set = set(os_dict)
            node_dict = {}

            for rxn in rxn_set:
                if os_dict[rxn]['direction'] == 'forward' and self.target in os_dict[rxn]['products']:
                    node_dict[rxn] = Node(rxn)
                    ordered_paths[count_pathway].set_root(node_dict[rxn])
                    rxn_set.remove(rxn)
                    break
                elif os_dict[rxn]['direction'] == 'reverse' and self.target in os_dict[rxn]['reactants']:
                    node_dict[rxn] = Node(rxn)
                    ordered_paths[count_pathway].set_root(node_dict[rxn])
                    rxn_set.remove(rxn)
                    break
            

            def add_all_children(node, rxn_set):
                ''' Adds all children of root and removes rxn from rxn_set '''
                if len(rxn_set) == 0:
                    return

                store_reactants = []
                if os_dict[node.root]['direction'] == 'forward':
                    store_reactants = list(os_dict[node.root]['reactants'].keys())
                else:
                    store_reactants = list(os_dict[node.root]['products'].keys())

                for sr in store_reactants:
                    for rxn in rxn_set:
                        if os_dict[rxn]['direction'] == 'forward' and sr in os_dict[rxn]['products']:
                            node_dict[rxn] = Node(rxn)
                            node.addChild(node_dict[rxn])
                            rxn_set.remove(rxn)
                            add_all_children(node_dict[rxn], rxn_set)
                            break
                        elif os_dict[rxn]['direction'] == 'reverse' and sr in os_dict[rxn]['reactants']:
                            node_dict[rxn] = Node(rxn)
                            node.addChild(node_dict[rxn])
                            rxn_set.remove(rxn)
                            add_all_children(node_dict[rxn], rxn_set)
                            break
                    if len(rxn_set) == 0:
                        break

            add_all_children(ordered_paths[count_pathway].root, rxn_set)

        return ordered_paths


    def alter_name_length(self, path_to_figure, cpdname):
        '''Shorten compound name if it is too long'''
        if len(path_to_figure) > 250:
            remove_variable = len(path_to_figure) - 250
            cpdname = cpdname[:-remove_variable]
        return cpdname

    def get_target_name(self):
        name = self.DB.get_compound_name(self.target)
        if name != 'None':
            return re.sub(':|/', '_', name)
        else:
            return re.sub(':|/', '_', self.target)

    def generate_cdxml_files(self, RP=None, ranktype=None, fba_fluxes=None, show_rxn_info=True):
        '''Produce cdxml files for each pathway'''   
        _FBA = fba_fluxes is not None
        if RP:
            correct_type_paths = RP.correctly_typed_paths
        else:
            correct_type_paths = [x for x in self.reactions]
        target_reformat_orig = self.get_target_name()
        self.generate_output_folders(self.output_path+'/solution_figures')
        self.generate_output_folders(self.output_path+'/solution_figures/'+self.target_organism_name)                        
        self.generate_output_folders(self.output_path+'/solution_figures/'+ self.target_organism_name + '/' +
            target_reformat_orig)
        self.generate_output_folders(self.output_path+'/compounds')
        if RP:
            best_paths_folder = self.output_path+'/solution_figures/'+ self.target_organism_name + '/' + target_reformat_orig+'/best_paths_for_sep_bp'
            self.generate_output_folders(best_paths_folder)
            other_paths_folder = self.output_path + '/solution_figures/' + self.target_organism_name + '/' + target_reformat_orig + '/other_nonideal_solutions'
            self.generate_output_folders(other_paths_folder)


        # For each pathway
        for count_pathway, os_dict in self.reactions.items():
            if count_pathway in correct_type_paths:
                if RP:
                    if count_pathway in RP.best_sp_paths:
                        MAIN_CE = CDXML_Editor(output_path=best_paths_folder + '/solution_' + str(count_pathway) + '.cdxml')
                    else:
                        MAIN_CE = CDXML_Editor(output_path=other_paths_folder + '/solution_' + str(count_pathway) + '.cdxml')
                else:
                    MAIN_CE = CDXML_Editor(output_path=self.output_path + '/solution_figures/' + self.target_organism_name + '/' + target_reformat_orig +
                                          '/solution_' + str(count_pathway) + '.cdxml')                
                path = self.ordered_paths[count_pathway]

                cdxml_reactants = {}
                cdxml_products = {}
                promiscuous_reactants = {}
                promiscuous_products = {}
                misc_products = {}

                reaction_proteins = {}
                reaction_solvents = {}
                reaction_catalysts = {}
                reaction_SPRESI_info = {}

                fba_values = {}

                BP = {}
                BP_predicted = {}

                logP = {}
                logP_predicted = {}

                def get_boiling_points(compounds):
                    for cpd in compounds:
                        cpd_name = compounds[cpd]
                        if cpd_name == 'None':
                            cpd_name = re.sub('_\w+\d+$', '', cpd)
                            cpd_name = '_'.join(cpd_name.split('/'))
                        try:
                            BP[cpd_name] = float(RP.CRV.sp_processed[cpd]['bp'])
                        except:
                            try:
                                BP_predicted[cpd_name] = float(RP.predicted_bps[cpd])
                            except:
                                pass

                def get_logP(compounds):
                    for cpd in compounds:
                        cpd_name = compounds[cpd]
                        if cpd_name == 'None':
                            cpd_name = re.sub('_\w+\d+$', '', cpd)
                            cpd_name = '_'.join(cpd_name.split('/'))
                        try:
                            logP[cpd_name] = float(RP.CRV.sp_processed[cpd]['logP'])
                        except:
                            try:
                              logP_predicted[cpd_name] = float(RP.predicted_logps[cpd])  
                            except:
                                pass

                def get_info(root, parent=None):
                    for child in root.children:
                        get_info(child, parent=root)

                    rxn = root.root
                    cdxml_reactants[rxn] = []
                    cdxml_products[rxn] = []
                    promiscuous_reactants[rxn] = []
                    promiscuous_products[rxn] = []
                    misc_products[rxn] = []

                    reaction_proteins[rxn] = None
                    reaction_solvents[rxn] = None
                    reaction_catalysts[rxn] = None
                    reaction_SPRESI_info[rxn] = None

                    if os_dict[rxn]['direction'] == 'forward':
                        cdxml_reactants[rxn], promiscuous_reactants[rxn] = self.get_cdxml(os_dict[rxn]['reactants'],
                                cdxml_reactants[rxn], promiscuous_reactants[rxn])
                        cdxml_products[rxn], promiscuous_products[rxn] = self.get_cdxml(os_dict[rxn]['products'],
                                cdxml_products[rxn], promiscuous_products[rxn])
                    else:
                        cdxml_reactants[rxn], promiscuous_reactants[rxn] = self.get_cdxml(os_dict[rxn]['products'],
                                cdxml_reactants[rxn], promiscuous_reactants[rxn])
                        cdxml_products[rxn], promiscuous_products[rxn] = self.get_cdxml(os_dict[rxn]['reactants'],
                                cdxml_products[rxn], promiscuous_products[rxn])

                    if _FBA:
                        fba_values[rxn] = float(max(0,fba_fluxes[rxn]))

                    if RP is not None and ranktype == 'bp':
                        get_boiling_points(os_dict[rxn]['reactants'])
                        get_boiling_points(os_dict[rxn]['products'])

                    if RP is not None and ranktype == 'logP':
                        get_logP(os_dict[rxn]['reactants'])
                        get_logP(os_dict[rxn]['products'])

                    org = self.reactions[count_pathway][rxn]['organisms'][0]
                    protein = self.DB.get_proteins(rxn, org)
                    if protein != "None":
                        try:
                            ec_number = re.search('\d+\.\d+\.\d+\.\d+',protein)[0]
                            self.enzyme_set.add(ec_number)
                        except:
                            pass
                        reaction_proteins[rxn] = protein                    
                    elif show_rxn_info:
                        # RXN SOLVENT
                        # solvents = self.DB.get_solvents(rxn)                        
                        # if solvents != "None" and len(solvents) > 0:
                        #     solvents = [re.sub('_\w+\d+$', '', s) for s in solvents]
                        #     if self.use_iupac_names:
                        #         for idx,s in enumerate(solvents):
                        #             try:
                        #                 solvents[idx] = pcp.get_properties('IUPACName',identifier=s,namespace='inchi')[0]['IUPACName']
                        #             except:
                        #                 pass
                        #     reaction_solvents[rxn] = "Solvent: %s" % ','.join(solvents)

                        # # RXN CATALYST
                        # catalysts = self.DB.get_catalysts(rxn)
                        # if catalysts != "None" and len(catalysts) > 0:
                        #     catalysts = [re.sub('_\w+\d+$', '', c) for c in catalysts]
                        #     if self.use_iupac_names:
                        #         for idx,c in enumerate(catalysts):
                        #             try:
                        #                 catalysts[idx] = pcp.get_properties('IUPACName',identifier=c,namespace='inchi')[0]['IUPACName']
                        #             except:
                        #                 pass
                        #     reaction_catalysts[rxn] = "Catalyst: %s" % ','.join(catalysts)
                        
                        # RXN SPRESI
                        spresi_info = []
                        # Reaxys RXN
                        if rxn.endswith('_RY'):
                            spresi_info += ['Reaxys Reaction ID: %s' % rxn.split('_')[1]]
                        else:
                            rxn_temperatures = self.DB.get_temperature(rxn)
                            rxn_pressures = self.DB.get_pressure(rxn)
                            rxn_times = self.DB.get_time(rxn)
                            rxn_yields = self.DB.get_yield(rxn)
                            if rxn_temperatures != "None" and len(rxn_temperatures) > 0:
                                spresi_info += ["Temperature: %s" % ','.join(rxn_temperatures)]
                            if rxn_pressures != "None" and len(rxn_pressures) > 0:
                                spresi_info += ["Pressure: %s" % ','.join(rxn_pressures)]
                            if rxn_times != "None" and len(rxn_times) > 0:
                                spresi_info += ["Time: %s" % ','.join(rxn_times)]
                            if rxn_yields != "None" and len(rxn_yields) > 0:
                                spresi_info += ["Yield: %s" % ','.join(rxn_yields)]
                        if len(spresi_info) > 0:
                            reaction_SPRESI_info[rxn] = '\n'.join(spresi_info)
                    

                    if os_dict[rxn]['direction'] == 'forward':
                        cur_prods = os_dict[rxn]['products']
                    else:
                        cur_prods = os_dict[rxn]['reactants']
                    
                    difference = []
                    if parent:
                        next_rxn = parent.root 
                        if os_dict[next_rxn]['direction'] == 'forward':
                            next_reacts = set(os_dict[next_rxn]['reactants'].keys())
                        else:
                            next_reacts = set(os_dict[next_rxn]['products'].keys())
                        difference = set(cur_prods.keys()).intersection(set(cur_prods.keys()).symmetric_difference(next_reacts))
                    else:
                        difference = set(cur_prods.keys()).intersection(set(cur_prods.keys()).symmetric_difference(set([self.target])))

                    for cpd_id in difference:
                        cpd = cur_prods[cpd_id]
                        if cpd == 'None':
                            cpd_name = re.sub('_\w+\d+$', '', cpd)
                            cpd_name = '_'.join(cpd_name.split('/'))
                            if cpd_name not in promiscuous_products[rxn]:
                                misc_products[rxn] += [cpd_id]
                        elif cpd not in promiscuous_products[rxn]:
                            misc_products[rxn] += [cpd]

                get_info(path.root)

                if _FBA:
                    max_flux = max(fba_values.values())
                    if max_flux != 0:
                        fba_values = {fv_key: float(fv_val/max_flux) for fv_key, fv_val in fba_values.items()}

                
                
                last_id=0
                color_index = 0
                cdxml_compounds_path = self.output_path + '/compounds/'

                def generate_cdxml(root, color_index=color_index,last_id=last_id):
                    previous = []
                    for child in root.children:
                        child_box, last_id = generate_cdxml(child, last_id=last_id)
                        previous += [child_box]

                    rxn = root.root
                    CE = CDXML_Editor(cdxml_files_path=cdxml_compounds_path, BP=BP, BP_predicted=BP_predicted, logP=logP, logp_predicted=logP_predicted)
                    
                    # If first reaction and only reactants are promiscuous compounds
                    if len(cdxml_reactants[rxn]) == 0 and len(previous) == 0:
                        cdxml_reactants[rxn] += [promiscuous_reactants[rxn][0]]
                        promiscuous_reactants[rxn] = promiscuous_reactants[rxn][1:]
                    

                    last_id = CE.add_reactants(cdxml_reactants[rxn], previous, last_id)

                    if reaction_proteins[rxn]:
                        CE.add_transition(promiscuous_reactants[rxn],promiscuous_products[rxn],misc_products[rxn], 
                                        reaction_proteins=reaction_proteins[rxn])
                    else:
                        CE.add_transition(promiscuous_reactants[rxn],promiscuous_products[rxn],misc_products[rxn], 
                                        reaction_solvents=reaction_solvents[rxn], reaction_catalysts=reaction_catalysts[rxn], reaction_SPRESI_info=reaction_SPRESI_info[rxn])
                    
                    if _FBA:
                        MAIN_CE.add_color(fba_values[rxn])
                        CE.set_FBA(color_index)
                        color_index += 1
                    
                    CE.set_products(cdxml_products[rxn])                
                    return CE, last_id

                pathway_cdxml, last_id = generate_cdxml(path.root)

                # FINALLY ADD TARGET
                target_CE = CDXML_Editor(cdxml_files_path=cdxml_compounds_path, BP=BP, BP_predicted=BP_predicted, logP=logP, logp_predicted=logP_predicted)
                target_name = self.DB.get_compound_name(self.target)
                if target_name == 'None':
                    target_name = re.sub('_\w+\d+$', '', self.target)
                    target_name = '_'.join(target_name.split('/'))       
                target_CE.add_product(target_name, last_id)
                pathway_cdxml.append(target_CE.container, arrange="right")

                MAIN_CE.append(pathway_cdxml.container)
                MAIN_CE.generate_file()
