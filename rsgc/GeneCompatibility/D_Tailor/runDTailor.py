from GeneCompatibility.D_Tailor import Functions
import sys
from GeneCompatibility.D_Tailor.SequenceDesigner import SequenceDesigner
# from GeneCompatibility.D_Tailor.Features.Structure import Structure,StructureMFE
from GeneCompatibility.D_Tailor.Features import CAI #,RNADuplex
from GeneCompatibility.D_Tailor.DesignOfExperiments.Design import Optimization#,RandomSampling,FullFactorial

class CAIDesigner(SequenceDesigner):
    
    def __init__(self, name, seed, design, dbfile, outputfile, cai_table, createDB=True):
        SequenceDesigner.__init__(self, name, seed, design, dbfile, outputfile, createDB)
        self.cai_table = cai_table
        
    def configureSolution(self, solution):
        '''
        Solution configuration
        '''

        if solution.sequence == None:
            return 0
        
        #Populate solution with desired features
        solution.mutable_region=list(range(0,len(solution.sequence))) # whole region
        solution.cds_region = (49,len(solution.sequence))
        solution.keep_aa = True
        
        cai_obj = CAI.CAI(solution=solution,label="cds",cai_table=self.cai_table,
                        args= {'cai_range':(0,len(solution.sequence)),'mutable_region':list(range(3,len(solution.sequence)))})
        solution.add_feature(cai_obj)

    def validateSolution(self, solution):
        '''
        Solution validation tests
        '''
        if solution.sequence == None or ('?' in list(solution.levels.values())):
            sys.stderr.write("SolutionValidator: Level unknown - "+str(solution.levels)+"\n")                        
            solution.valid = False
            return 0
        
        #check if solution is valid
        valid = True
      
        designed_region = solution.sequence
                
        #No internal Promoters
        (score, _, _) = Functions.look_for_promoters(designed_region)
        if score >= 15.3990166: #0.95 percentile for Promoter PWM scores
            valid = False
            sys.stderr.write("SolutionValidator: High Promoter score: "+str(score)+"\n")                    
        
        #No internal Terminator
        score = Functions.look_for_terminators(designed_region)
        if score >= 90: #90% confidence from transtermHP
            valid = False
            sys.stderr.write("SolutionValidator: High Terminator score\n")
            
        #No restriction enzymes
        if 'ggtctc' in designed_region or 'gagacc' in designed_region:
           sys.stderr.write("SolutionValidator: Restriction enzyme found\n")
           valid = False        
        
        solution.valid = valid
        
        return valid


# class CAIEcoliDesigner(SequenceDesigner):
    
#     def __init__(self, name, seed, design, dbfile, outputfile, createDB=True):
#         SequenceDesigner.__init__(self, name, seed, design, dbfile, outputfile, createDB)
        
#     def configureSolution(self, solution):
#         '''
#         Solution configuration
#         '''

#         if solution.sequence == None:
#             return 0
        
#         #Populate solution with desired features
#         solution.mutable_region=list(range(0,len(solution.sequence))) # whole region
#         solution.cds_region = (49,len(solution.sequence))
#         solution.keep_aa = True
        
#         cai_obj = CAI.CAI(solution=solution,label="cds",args= {'cai_range' : (0,len(solution.sequence)),
#                                                                'mutable_region' : list(range(3,len(solution.sequence))) } )
#         solution.add_feature(cai_obj)

#     def validateSolution(self, solution):
#         '''
#         Solution validation tests
#         '''
#         if solution.sequence == None or ('?' in list(solution.levels.values())):
#             sys.stderr.write("SolutionValidator: Level unknown - "+str(solution.levels)+"\n")                        
#             solution.valid = False
#             return 0
        
#         #check if solution is valid
#         valid = True
      
#         designed_region = solution.sequence
                
#         #No internal Promoters
#         (score, _, _) = Functions.look_for_promoters(designed_region)
#         if score >= 15.3990166: #0.95 percentile for Promoter PWM scores
#             valid = False
#             sys.stderr.write("SolutionValidator: High Promoter score: "+str(score)+"\n")                    
        
#         #No internal Terminator
#         score = Functions.look_for_terminators(designed_region)
#         if score >= 90: #90% confidence from transtermHP
#             valid = False
#             sys.stderr.write("SolutionValidator: High Terminator score\n")
            
#         #No restriction enzymes
#         if 'ggtctc' in designed_region or 'gagacc' in designed_region:
#            sys.stderr.write("SolutionValidator: Restriction enzyme found\n")
#            valid = False        
        
#         solution.valid = valid
        
#         return valid

