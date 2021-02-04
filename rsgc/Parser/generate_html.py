from __future__ import print_function
__author__ = 'Leanne Whitmore'
__email__ = 'leanne382@gmail.com'
__description__ = 'Generate html output'

import os
import re
import glob
import base64
from rsgc.Database import query as Q
from sys import platform
if platform == "cygwin":
    header_path = 'C:\cygwin64'
else:
    header_path = ''

PATH = os.path.dirname(os.path.abspath(__file__))
legend=[]
with open(os.path.join(PATH, "_static", "legend.txt")) as fin:
    for line in fin:
        line = line.strip()
        legend.append(line)
legend = "".join(legend)


class HtmlOutput(object):
    def __init__(self, eval_targets, output_path, fba, figures, database, outputfile):
        self.eval_targets = eval_targets
        self.output_path = output_path
        self.fba = fba
        self.figures = figures
        self.DB = Q.Connector(database)
        self.compoundpaths = glob.glob(os.path.join(self.output_path, "raw_compound_solutions", "*.txt"))
        self.output = open(outputfile, 'w')
        self.output.write("<!DOCTYPE html>\n\n")
        self.output.write("<html>\n")
        self.load_head_info()
        self.read_css_info()
        self.read_js_info()
        self.load_body_info()
        self.load_general_stats()
        self.load_pathways()
        self.load_individual_path_info()
        self.write_collapsible()
        self.output.write("</head>\n")
        self.output.write("</html>\n")
        self.output.close()

    def load_fba(self):
        temp_dict=dict()
        with open(os.path.join(self.output_path, "theoretical_yield.txt")) as fin:
            for line in fin:
                line = line.strip()
                larray = line.split("\t")
                cpd = larray[0].split("---")
                org = larray[1].split("---")
                temp_dict.setdefault(cpd[0], {})
                temp_dict[cpd[0]].setdefault(org[1], {})
                temp_dict[cpd[0]][org[1]]["bio"]=larray[-1]
                temp_dict[cpd[0]][org[1]]["tar"]=larray[-3]
        return(temp_dict)

    def write_collapsible(self):
        self.output.write("<script>")
        self.output.write("var coll = document.getElementsByClassName(\"collapsible\");\n")
        self.output.write("var i;\n")

        self.output.write("for (i = 0; i < coll.length; i++) {\n")
        self.output.write("coll[i].addEventListener(\"click\", function() {\n")
        self.output.write("this.classList.toggle(\"active\");\n")
        self.output.write("var content = this.nextElementSibling;\n")
        self.output.write("if (content.style.display === \"block\") {\n")
        self.output.write("content.style.display = \"none\";\n")
        self.output.write("} else {\n")
        self.output.write("content.style.display = \"block\";\n")
        self.output.write("}\n")
        self.output.write("});\n")
        self.output.write("}\n")
        self.output.write("</script>")

    def retrieve_cpd_name(self, cpdID):
        if self.DB.get_compound_name(cpdID) is None:
            return i
        else:
            return self.DB.get_compound_name(cpdID)        

    def load_individual_path_info(self):
        if self.fba:
            self.fba_dict = self.load_fba()
        self.output.write("<h3>Novel pathway information for targets:</h3>\n")
        for i in self.opt_path:
            cpdname = self.retrieve_cpd_name(i)
            for o in self.opt_path[i]:
                self.output.write("<button type=\"button\" class=\"collapsible\"><a id=\"{}\">{}_{}</a></button>\n".format(cpdname, cpdname, o))
                self.output.write("<div class=\"content\">\n")
                self.output.write("<p>{} pathways in {}</p>\n".format(cpdname, o))
                if self.figures:
                    itemp = re.sub("-|/|,|\)|\(|<|>|:|\]|\[|\/", "_", i)
                    itemp = re.sub(' ', '_', itemp)
                    encoded = base64.b64encode(open(os.path.join(self.output_path, "solution_figures", "SC_graph_{}_{}.png".format(itemp , o)), "rb").read()).decode('utf-8')
                    self.output.write("<img src=\"data:image/png;base64, {}\" style=\"width:100%\"></img>\n".format(encoded))
                if self.figures and self.fba:
                    self.output.write("<p>Flux through pathways legend</p>\n")
                    self.output.write("<img src=\"data:image/png;base64, {}\" style=\"width:50%\"></img>\n".format(legend))

                self.output.write("<h4>Pathway information:</h4>\n")
                self.output.write("<table style=\"width:100%\">\n")
                self.output.write("<tr>\n")
                self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Pathway</th>\n")
                self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Reaction</th>\n")
                self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Reaction Name</th>\n")
                self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Enzymes</th>\n")
                self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Reactants</th>\n")
                self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Products</th>\n")
                self.output.write("</tr>\n")
                for s in self.opt_path[i][o]['sol']:
                    for r in self.opt_path[i][o]['sol'][s]['rxn']:
                        self.output.write("<tr>\n")
                        self.output.write("<td style=\"text-align: center; vertical-align: middle; font-size: 13px;\">%s</td>\n" % s)
                        self.output.write("<td style=\"text-align: center; vertical-align: middle;font-size: 13px\">%s</td>\n" % r)
                        self.output.write("<td style=\"text-align: center; vertical-align: middle; font-size: 13px\">%s</td>\n" % self.opt_path[i][o]['sol'][s]['rxn'][r]["name"])
                        if os.path.exists(self.output_path+'/gene_compatibility/geneseqs_{}_{}.txt'.format(self.opt_path[i][o]['sol'][s]['rxn'][r]["enz"],  self.DB.get_model_ID(o))):
                            file_path=header_path+os.path.abspath(self.output_path+"gene_compatibility/geneseqs_{}_{}.txt".format(self.opt_path[i][o]['sol'][s]['rxn'][r]["enz"], self.DB.get_model_ID(o)))
                            self.output.write("<td style=\"text-align: center; vertical-align: middle;font-size: 13px\"><a href=file:///{} target=\"popup\"\">{}</a></td>\n".format(file_path, self.opt_path[i][o]['sol'][s]['rxn'][r]["enz"]))
                        else:
                            self.output.write("<td style=\"text-align: center; vertical-align: middle;font-size: 13px\">%s</td>\n" % self.opt_path[i][o]['sol'][s]['rxn'][r]["enz"])
                        self.output.write("<td style=\"text-align: center; vertical-align: middle;font-size: 13px\">%s</td>\n" % ",".join(self.opt_path[i][o]['sol'][s]['rxn'][r]["react"]))
                        self.output.write("<td style=\"text-align: center; vertical-align: middle;font-size: 13px\">%s</td>\n" % ",".join(self.opt_path[i][o]['sol'][s]['rxn'][r]["prod"]))
                        self.output.write("</tr>\n")
                self.output.write("</table>\n")
                if self.fba:
                    self.output.write("<h3>Theoretical yield information (FBA) results</h3>\n")
                    self.output.write("<p>%s</p>\n" % self.fba_dict[i][o]["tar"])
                    self.output.write("<p>%s</p>\n" % self.fba_dict[i][o]["bio"])
                self.output.write("</div>\n")

    def open_optimal_pathways(self):
        self.opt_path = dict()
        rxncounter=0
        rxnid = ""
        with open(os.path.join(self.output_path, "optimal_pathways.txt"), 'r') as fin:
            for line in fin:
                line = line.strip("\n")
                if line.startswith("SHORTEST PATH FOR "):
                    line = re.sub("SHORTEST PATH FOR ", "", line)
                    larray = line.split(" ")
                    if larray[1] != "in":
                        self.opt_path.setdefault(larray[0], {})
                        self.opt_path[larray[0]].setdefault(larray[5], {})
                        self.opt_path[larray[0]][larray[5]].setdefault("sol", {})
                        orgid = larray[5]
                    else:
                        self.opt_path.setdefault(larray[0], {})
                        self.opt_path[larray[0]].setdefault(larray[4], {})
                        self.opt_path[larray[0]][larray[4]].setdefault("sol", {}) 
                        orgid = larray[4]
                elif line.startswith("Solution"):
                    rxncounter=0
                    temparray = line.split(" ")
                    self.opt_path[larray[0]][orgid]["sol"].setdefault(temparray[1], {})
                    self.opt_path[larray[0]][orgid]["sol"][temparray[1]].setdefault("rxn", {})
                elif line.endswith("reactants") or line.endswith("products"):
                    line = line.strip()
                    if line.endswith("reactants"):
                        line = re.sub(" reactants", "", line)
                        react = line.split("\t")
                        if len(react) >2:
                            self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid]["react"].append(react[2])
                        else:
                            self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid]["react"].append(react[1])
                    if line.endswith("products"):
                        line = re.sub(" products", "", line)
                        prod = line.split("\t")
                        if len(react) >2:
                            self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid]["prod"].append(prod[2])
                        else:
                            self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid]["prod"].append(prod[1])
                elif line.startswith("No paths"):
                    pass
                elif line != "":
                    rxncounter+=1
                    rxnarray = line.split("\t")
                    rxnid=str(rxncounter)+"_"+rxnarray[0]
                    self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"].setdefault(rxnid, {})
                    self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid].setdefault("react", [])
                    self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid].setdefault("prod", [])
                    self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid].setdefault("enz", str(rxnarray[3]))
                    self.opt_path[larray[0]][orgid]["sol"][temparray[1]]["rxn"][rxnid].setdefault("name", str(rxnarray[1]))

    def load_pathways(self):
        self.output.write("<div class=\"section\" id=\"retSynth-results\">\n")
        self.output.write("<h1>Pathways</h1>\n")
        self.output.write("<table style=\"width:90%\">\n")
        self.output.write("<tr>\n")
        self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Compounds</th>\n")
        self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\">Organism</th>\n")
        self.output.write("<th style=\"text-align: center; vertical-align: middle; background-color: #75c47c; color: black\"># of Pathways</th>\n")
        self.output.write("</tr>\n")
        self.open_optimal_pathways()
        for i in self.opt_path:
            cpdname = self.retrieve_cpd_name(i)
            for o in self.opt_path[i]:
                self.output.write("<tr>\n")
                self.output.write("<td style=\"text-align: center; vertical-align: middle;\"><a href=\"#{}\">{}</a></td>\n".format(cpdname, cpdname))
                self.output.write("<td style=\"text-align: center; vertical-align: middle;\">%s</td>\n" % o) 
                self.output.write("<td style=\"text-align: center; vertical-align: middle;\">%s</td>\n" % len(self.opt_path[i][o]['sol']))
                self.output.write("</tr>\n")
        self.output.write("</table>\n")
    def load_general_stats(self):
        self.output.write("<p>Number of targets evaluated = %s </p>\n" % self.eval_targets)
        self.output.write("<p>Number of targets with solutions = %s </p>\n" % len(self.compoundpaths))
        # self.output.write("<\div>\n")
    
    def load_body_info(self):
        self.output.write("<body>\n")
        self.output.write("<div class=\"relbar-top\">\n")
        self.output.write("<div class=\"related\" role=\"navigation\" aria-label=\"related navigation\">\n")
        self.output.write("<h3>Navigation</h3>\n")
        self.output.write("<ul>\n")
        self.output.write("<li class=\"nav-item nav-item-this\"><a href=\"\">Results for RetSynth</a></li>\n") 
        self.output.write("</ul>\n")
        self.output.write("</div>\n")
        self.output.write("</div>\n")
        self.output.write("<div class=\"document\">\n")
        self.output.write("<div class=\"documentwrapper\">\n")
        self.output.write("<div class=\"bodywrapper\">\n")
        self.output.write("<div class=\"body\" role=\"main\">\n")
        self.output.write("<div class=\"section\" id=\"retSynth-results\">\n")
        self.output.write("<h1>RetSynth General Stats</h1>\n")

    def read_css_info(self):
        cssfiles = glob.glob(os.path.join(PATH, "_static", "*.css"))
        for cssfile in cssfiles:
            self.output.write("<style>")
            with open(cssfile, 'r') as fin:
                for line in fin:
                    line = line.strip()
                    if re.search("\@import", line):
                        pass 
                    else:
                        self.output.write(line)
            self.output.write("</style>\n")
    
    def read_js_info(self):
        jsfiles = glob.glob(os.path.join(PATH, "_static", "*.js"))
        for jsfile in jsfiles:
            self.output.write("<script>\n")
            with open(jsfile, 'r') as fin:
                for line in fin:
                    if re.search("\@import", line):
                        pass 
                    else:
                        self.output.write(line)
            self.output.write("</script>\n")
    
    def load_head_info(self):
        self.output.write("<head>\n")
        self.output.write("<meta charset=\"utf-8\" />\n")
        self.output.write("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n")
        self.output.write("<title>RetSynth Results!</title>\n")