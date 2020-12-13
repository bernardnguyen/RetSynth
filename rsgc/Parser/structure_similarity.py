__author__ = 'Leanne Whitmore'
__email__ = 'lwhitmo@sandia.gov'
__description__ = 'Check based on compound structure if it is in the database'

import re
from tqdm import tqdm
from copy import deepcopy
from sys import platform
from rsgc.Database import query as Q
if platform == 'darwin':
    from rsgc.indigopython130_mac import indigo
    from rsgc.indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from rsgc.indigopython130_linux import indigo
    from rsgc.indigopython130_linux import indigo_inchi
elif platform == "win32" or platform == "win64" or platform == "cygwin":
    from rsgc.indigopython130_win import indigo
    from rsgc.indigopython130_win import indigo_inchi

def verbose_print(verbose, line):
    if verbose:
        print(line)

class TanimotoStructureSimilarity(object):
    """Identifies if compounds are in the database but under a different ID
    based on structure similarity"""
    def __init__(self, targets, database, all_compounds, cytosol, extracellular, verbose, threshold_score = 1):
        """Initialize"""
        self.verbose = verbose
        self.targets = []
        self.track_final_cpd = []
        self.DB = Q.Connector(database)
        '''Remove duplicated targets'''
        for target in targets:
            if target not in self.targets:
                self.targets.append(target)
        self.individualtargets = []
        self.organisms = []
        for target in self.targets:
            self.individualtargets.append(target[0])
            organism = ','.join([target[1], target[2], target[3]])
            if organism not in self.organisms:
                self.organisms.append(organism)
        self.organisms = list(set(self.organisms))
        self.all_compounds = all_compounds
        self.cytosol = cytosol
        self.extracellular = extracellular
        self.threshold_score = float(threshold_score)
        self.calculate_tanimoto_score()

    def process_targets(self, targets):
        """Reformat targets to remove compartment information"""
        reformat_target = []
        for target in targets:
            targetinchi = self.DB.get_inchi_from_compoundID(target)
            if targetinchi is not None:
                reformat_target.append(targetinchi)
        return(reformat_target)

    # def get_original_target(self, tmol):
    #     """Obtain original target"""
    #     index = None
    #     for count, target in enumerate(self.targets):
    #         if tmol in target:
    #             index =  count
         
    #     if not index: 
    #         verbose_print(self.verbose, 'WARNING:\tCould not get target '+tmol)
    #         return (None)
    #     else:
    #         return index

    def get_tanimoto_score(self, tmol, threshold, IN):
        temp = {}
        # print (threshold)
        for db_cpd in self.db_cpds_fp:
            score = IN.similarity(self.individualtargets_p_fp[tmol], self.db_cpds_fp[db_cpd], 'tanimoto')
            temp[db_cpd] = score
        # print (temp)
        max_score_cpds = [i for i, v in list(temp.items()) if float(v) >= float(threshold)]
        return (max_score_cpds)

    def calculate_tanimoto_score(self):
        """Calculate similarity for compounds that are not in the database"""
           
        self.finaltargets = self.targets
        self.individualtargets_p = self.process_targets(set(self.individualtargets))
        IN = self.retrieve_fingerprints()
        for count ,tmol in enumerate(self.individualtargets_p_fp.keys()):
            # print (tmol)
            if tmol in self.all_compounds:
                print (tmol)
                max_score_cpds = self.get_tanimoto_score(tmol, self.threshold_score, IN)
                if max_score_cpds:
                    verbose_print(self.verbose, 'STATUS:\t{} compounds have {} or greater similarity to target compound {}'.format(len(set(max_score_cpds)), float(self.threshold_score)*100, tmol))
                    for max_score_cpd in set(max_score_cpds):
                        if max_score_cpd != tmol:
                            new_target = self.DB.get_compound_ID_from_inchi(max_score_cpd)
                            if new_target not in self.individualtargets and new_target not in self.track_final_cpd:
                                self.fill_final_targets(new_target)

            else:
                max_score_cpds = self.get_tanimoto_score(tmol, self.threshold_score, IN)
                if max_score_cpds:
                    # index = self.get_original_target(tmol)
                    # if index:
                    #     del self.finaltargets[index]
                    verbose_print(self.verbose, 'STATUS:\t{} compounds have {} or greater similarity to target compound {}'.format(len(set(max_score_cpds)), float(self.threshold_score)*100, tmol))
                    for max_score_cpd in set(max_score_cpds):
                        new_target = self.DB.get_compound_ID_from_inchi(max_score_cpd)
                        if new_target not in self.individualtargets and new_target not in self.track_final_cpd:
                            verbose_print(self.verbose, 'STATUS:\tAdding compound {} to target list'.format(max_score_cpd))
                            self.fill_final_targets(new_target)
                else:
                    verbose_print(self.verbose,'STATUS:\tNo compounds in the database that are {} percent similar to target {}'.format(float(self.threshold_score)*100, tmol+'_'+self.cytosol))

    def fill_final_targets(self, new_target):
        for organism in self.organisms:
            organisms = organism.split(',')
            self.finaltargets.append([new_target]+organisms)
            self.track_final_cpd.append(new_target)

    def retrieve_fingerprints(self):
        """Retrieve fingerprints for database compounds and targets"""
        IN = indigo.Indigo()
        INCHI = indigo_inchi.IndigoInchi(IN)

        print(self.individualtargets_p)
        self.individualtargets_p_fp = {}
        verbose_print(self.verbose, "STATUS:\tgetting fingerprints for target compounds")
        for count, tmol in enumerate(tqdm(self.individualtargets_p)):
            try:
                self.individualtargets_p_fp[tmol] = INCHI.loadMolecule(tmol).fingerprint('full')
            except indigo.IndigoException:
                verbose_print(self.verbose, 'WARNING:\tCould not get fingerprint for {}'.format(tmol))

        self.db_cpds_fp = {}
        db_cpds_set = set(self.all_compounds)
        verbose_print(self.verbose, "STATUS:\tgetting fingerprints for database compounds")        
        for db_cpd in tqdm(db_cpds_set):
            try:
                mol = INCHI.loadMolecule(db_cpd)
                self.db_cpds_fp[db_cpd] = mol.fingerprint('full')
            except indigo.IndigoException:
                pass
        return IN
