import os, sys, multiprocessing
from threading import Thread
from rs_gui.stdoutqueue import StdoutQueue
from rs_gui.tkintertools import *

import tkinter as tk
import tkinter.filedialog
import tkinter.ttk
from tkinter.scrolledtext import ScrolledText as tkst
from tkinter.messagebox import askokcancel

from sys import platform
if platform == 'darwin':
    from indigopython130_mac import indigo
    from indigopython130_mac import indigo_inchi
elif platform == "linux" or platform == "linux2":
    from indigopython130_linux import indigo
    from indigopython130_linux import indigo_inchi
elif platform == "win32" or platform == 'win64':
    from indigopython130_win import indigo
    from indigopython130_win import indigo_inchi
IN = indigo.Indigo()
II = indigo_inchi.IndigoInchi(IN)

from RetSynthIO import *
from Database import mackinac
HOMEDIR = os.path.expanduser('~')
DEFAULT_PATRICFILE = sys.path[0] + '/Database/data/PATRIC_genome_complete_07152018.csv'


# UI COLOR DEFINITIONS
green = 'chartreuse3'
blue = 'cyan2'
black = 'gray10'
white = 'floral white'

class Window(object):
	def __init__(self, title='', width=100, height=100, background=white):
		self.window = tk.Tk()
		self.window.title(title)
		self.window.geometry('%dx%d' % (width,height))

		self.style = tkinter.ttk.Style(self.window)
		self.style.theme_use('clam')
		self.window['background'] = background

		self.initialize()
		self.window.mainloop()

	def initialize(self):
		print("Initialize has not been implemented")		

	def init_colors(self):
		tkinter.ttk.Style().configure('.', font='Verdana 14')
		tkinter.ttk.Style().configure('green.TFrame', background=green)
		tkinter.ttk.Style().configure('green.TLabel', foreground=white, background=green)

		tkinter.ttk.Style().configure('white.TFrame', background=white)
		tkinter.ttk.Style().configure('white.TButton', background=white, foreground=black, border=0)
		tkinter.ttk.Style().map('white.TButton', background=[('active',white)])
		tkinter.ttk.Style().map('white.TEntry', foreground=[('disabled',black)], background=[('disabled',white)])

		tkinter.ttk.Style().configure('black.TFrame', background=black, foreground=green)
		tkinter.ttk.Style().configure('black.TLabel', foreground=black, background=green)
		tkinter.ttk.Style().configure('errorBlack.TLabel', foreground=green, background=black)
		tkinter.ttk.Style().configure('black.TButton', foreground=green, background=black, border=0, font='Verdana 14')
		tkinter.ttk.Style().map('black.TButton', background=[('active',black)])

		tkinter.ttk.Style().configure('white.TLabel', background=white, foreground=black)
		tkinter.ttk.Style().configure('white.TCheckbutton', background=white, foreground=black)
		tkinter.ttk.Style().map('white.TCombobox',fieldbackground=[('readonly',white)],selectforeground=[('readonly',black)],
			selectbackground=[('readonly',white)])

class TopLevel(object):
	def __init__(self, parent_window, background=white):
		self.parent = parent_window
		self.window = tk.Toplevel(self.parent)
		self.window.title('')		

		self.style = tkinter.ttk.Style(self.window)
		self.style.theme_use('clam')
		self.window['background'] = background		

		self.mainframe = tkinter.ttk.Frame(self.window, style='white.TFrame')		
		
		self.initialize()
		tkinter.ttk.Style().configure('.', font='Verdana 12')
		self.window.mainloop()

	def initialize(self):
		print('Initialize has not been implemented')


class ErrorWindow(TopLevel):
	def __init__(self,parent_window,message,background=black):
		self.message = message
		super(ErrorWindow,self).__init__(parent_window,background=background)		
		
	def initialize(self):
		self.window.resizable(0,0)
		self.titleframe = tkinter.ttk.Frame(self.window,style='green.TFrame')
		self.titleframe.pack(fill='x')

		self.mainframe.style = 'black.TFrame'
		self.mainframe.pack(fill='both',expand=True,pady=(15,20),padx=20)		

		tkinter.ttk.Label(self.titleframe,text='Error - Whoops!',font='Impact 40',style='green.TLabel').pack()
		tkinter.ttk.Label(self.mainframe,text=self.message,font='Verdana 14',style='errorBlack.TLabel',justify=tk.CENTER).pack()


class MainWindow(Window):
	def __init__(self, DB_connector, default_database, models, width=800, height=600):
		self.DB = DB_connector
		self.default_database = default_database
		self.default_models = models
		super(MainWindow,self).__init__(width=width, height=height)
		
	def initialize(self):
		self.init_colors()	

		### DEFAULT RETSYNTH PARAMETERS ###			
		self.window.target = tk.StringVar()
		self.window.target_options = ['Name','InChI','KEGG','CAS','Chemical formula']
		self.window.default_models = self.default_models

		self.window.output_path = tk.StringVar()
		self.window.output_path.set('')
		self.window.output_check = 1
		self.window.kparam = 0
		self.window.evalrxns = 'all'
		self.window.oxf_check = 0
		self.window.tanimoto_check = 0
		self.window.tanimoto_thresh = 0.90
		self.window.fba_check = 0
		self.window.ko_check = 0
		self.window.gene_compatibility_check = 0
		self.window.cai_optimal_threshold = 0.5
		self.window.iupac_name_check = 0
		self.window.chemical_rxn_info_check = 1		
		
		self.window.default_database = self.default_database
		self.window.database_file = self.default_database
		self.window.default_db_check = 1
		self.window.existing_db_check = 0
		self.window.generate_db_check = 0

		self.window.gdbname = ''
		self.window.inchi_check = 0
		self.window.patric_check = 0
		self.window.patric_un = ''
		self.window.patric_pw = ''
		self.window.patric_media = 'Complete'
		self.window.patric_file = ''
		self.window.metacyc_check = 0
		self.window.metacyc_addition = ''
		self.window.kegg_check = 0
		self.window.kegg_org_type = 'bacteria'
		self.window.kegg_nunorgs = 0
		self.window.kegg_nunorg_paths = 0
		self.window.proxy = ''


		# TITLE FRAME
		titleframe = tkinter.ttk.Frame(self.window,style='green.TFrame')
		titleframe.pack(fill='x',ipady=10)
		tkinter.ttk.Label(titleframe,text='RetSynth',font='Impact 60',style='green.TLabel').pack(expand=True)

		# MODEL/TARGET FRAME
		mainframe = tkinter.ttk.Frame(self.window,style='white.TFrame')
		mainframe.pack(fill='x',pady=(60,0))
		modelframe = tkinter.ttk.Frame(mainframe, style='white.TFrame')
		modelframe.pack(padx=50, fill='x')
		targetframe = tkinter.ttk.Frame(mainframe, style='white.TFrame')
		targetframe.pack(pady=20,padx=50, fill='x')
		maintargetframe = tkinter.ttk.Frame(targetframe, style='white.TFrame')
		maintargetframe.pack(fill='x')

		# Model Selection
		tkinter.ttk.Label(modelframe,text='Starting Model:',font='Verdana 14',style='white.TLabel').pack(anchor='center')
		self.window.start = tkinter.ttk.Combobox(modelframe,values=self.window.default_models,state='readonly',style='white.TCombobox')
		self.window.start.current(0)
		self.window.start.pack(fill='x', expand=True, ipady=5)

		# Target
		self.window.target_label = tk.StringVar()
		self.window.target_label.set('Target (%s):' % '/'.join(self.window.target_options))
		tkinter.ttk.Label(maintargetframe,textvariable=self.window.target_label,font='Verdana 14',style='white.TLabel').pack(anchor='center')		
		target_entry = tkinter.ttk.Entry(maintargetframe, textvariable=self.window.target, justify='left', style='white.TEntry')
		target_entry.pack(side='left',fill='x', expand=True, ipady=5)
		target_entry.bind('<KeyRelease>', self.on_keyrelease)
		tkinter.ttk.Button(maintargetframe,text='Browse',width=7,command=self.get_target_from_file,
			style='black.TButton').pack(side='right')
		self.autocomplete = tk.Listbox(targetframe)

		# BOTTOM FRAME
		bottomframe = tkinter.ttk.Frame(self.window,style='black.TFrame')
		bottomframe.pack(side='bottom',fill='x')
		tkinter.ttk.Button(bottomframe, text='Settings',width=10,command=self.settingspanel,
			style='black.TButton').pack(side='left',padx=20,pady=20)
		tkinter.ttk.Button(bottomframe, text='Run', width=10, command=self.run_rs,
			style='black.TButton').pack(side='right',padx=20,pady=20)

	def get_target_from_file(self):
		file = tkinter.filedialog.askopenfilename(title="Select Target Compound File (mol/sdf)",initialdir=HOMEDIR,
			filetypes=[("Mol files","*.mol"),("SDF files","*.sdf")])
		if file:
			try:
				m = IN.loadMoleculeFromFile(file)
				inchi = II.getInchi(m)		
			except:
				inchi = False

			if inchi:
				self.window.target_label.set('Target:')
				self.window.target.set('InChI from file -- %s' % inchi)
			else:
				self.window.target.set('ERROR: Could not retrieve InChI from file')

	### AUTOCOMPLETE FUNCTIONS ###
	def on_keyrelease(self,event):
		if self.window.gdbname:
			self.autocomplete.delete(0, 'end')
			self.autocomplete.pack_forget()
			return
		# get text from entry
		value = event.widget.get()
		if len(value) >= 4:
			data = find_matches(value, database=self.window.database_file)               
			# update data in listbox
			self.listbox_update(data)
		else:
			# delete previous data
			self.autocomplete.delete(0, 'end')
			self.autocomplete.pack_forget()

	def on_select(self,event):
		selection = event.widget.curselection()
		if selection:
			self.window.target.set(str(event.widget.get(selection)))
			value = self.window.target.get()
			data = find_matches(value, database=self.window.database_file)
			self.autocomplete.delete(0,'end')
			self.autocomplete.pack_forget()

	def listbox_update(self,data):
		if len(data) > 0:
			self.autocomplete.pack(fill='both', expand=True)
			self.autocomplete.bind('<<ListboxSelect>>', self.on_select)	
		# delete previous data
		self.autocomplete.delete(0, 'end')
		# put new data
		for cpd_name, cpd_id in data:
			if cpd_name == 'None':
				cpd_name = 'Name N/A'
			self.autocomplete.insert('end', '%s -- %s' % (cpd_name,cpd_id))

	def settingspanel(self):
		SettingsWindow(parent_window=self.window)

	def run_rs(self):
		request_form = {}

		if self.window.start.get() and ' -- ' in self.window.target.get() and self.window.target.get().split(' -- ')[1] != '':
			request_form["output_path"] = self.window.output_path.get()
			request_form["start"] = self.window.start.get()
			request_form["target"] = self.window.target.get()
			
			request_form["kparam"] = self.window.kparam
			request_form["evalrxns"] = self.window.evalrxns
			request_form["oxf_check"] = self.window.oxf_check
			request_form["tanimoto_check"] = self.window.tanimoto_check
			request_form["tanimoto_thresh"] = self.window.tanimoto_thresh
			request_form["fba_check"] = self.window.fba_check
			request_form["ko_check"] = self.window.ko_check
			request_form["gene_compatibility_check"] = self.window.gene_compatibility_check
			request_form["cai_optimal_threshold"] = self.window.cai_optimal_threshold
			request_form["iupac_name_check"] = self.window.iupac_name_check
			request_form["chemical_rxn_info_check"] = self.window.chemical_rxn_info_check			

			request_form["database_file"] = self.window.database_file
			request_form["gdbname"] = self.window.gdbname
			request_form["inchi_check"] = self.window.inchi_check
			request_form["patric_check"] = self.window.patric_check
			request_form["patric_un"] = self.window.patric_un
			request_form["patric_pw"] = self.window.patric_pw
			request_form["patric_media"] = self.window.patric_media
			request_form["patric_file"] = self.window.patric_file
			request_form["metacyc_check"] = self.window.metacyc_check
			request_form["metacyc_addition"] = self.window.metacyc_addition
			request_form["kegg_check"] = self.window.kegg_check
			request_form["kegg_org_type"] = self.window.kegg_org_type
			request_form["kegg_nunorgs"] = self.window.kegg_nunorgs
			request_form["kegg_nunorg_paths"] = self.window.kegg_nunorg_paths			

			RunWindow(request_form)
		else: 
			ErrorWindow(parent_window=self.window,message='Did you select a starting model and a target compound?')

class SettingsWindow(TopLevel):
	def __init__(self, parent_window):
		super(SettingsWindow,self).__init__(parent_window)
	
	def initialize(self):
		self.window.geometry('650x900')

		self.mainframe.pack(side='right', fill='both', expand=True)
		self.sideframe = tkinter.ttk.Frame(self.window,style='green.TFrame')
		self.sideframe.pack(fill='y',side='left')
		
		### MAIN WINDOW
		self.pathwayframe = tkinter.ttk.Frame(self.mainframe,style='white.TFrame')
		self.outputframe = tkinter.ttk.Frame(self.mainframe,style='white.TFrame')
		self.databaseframe = tkinter.ttk.Frame(self.mainframe,style='white.TFrame')
		self.networkframe = tkinter.ttk.Frame(self.mainframe,style='white.TFrame')

		self.pathwayframe.pack(padx=20,pady=20,fill='both',expand=True)

		### SIDE MENU
		self.pathway_label = tkinter.ttk.Label(self.sideframe,text='Pathway',font='Verdana 16',
			style='black.TLabel')
		self.pathway_label.pack(padx=10,pady=(20,10),anchor='e')
		self.pathway_label.bind('<Button-1>',lambda evt, whichframe=self.pathwayframe: 
                            self.frame_change(self.pathwayframe))
		
		self.output_label = tkinter.ttk.Label(self.sideframe,text='Output',font='Verdana 16',
			style='green.TLabel')
		self.output_label.pack(padx=10,pady=10,anchor='e')
		self.output_label.bind('<Button-1>',lambda evt, whichframe=self.outputframe: 
                            self.frame_change(self.outputframe))

		self.database_label = tkinter.ttk.Label(self.sideframe,text='Database',font='Verdana 16',
			style='green.TLabel')
		self.database_label.pack(padx=10,pady=10,anchor='e')
		self.database_label.bind('<Button-1>',lambda evt, whichframe=self.databaseframe: 
                            self.frame_change(self.databaseframe))

		self.network_label = tkinter.ttk.Label(self.sideframe,text='Network',font='Verdana 16',
			style='green.TLabel')
		self.network_label.pack(padx=10,pady=10,anchor='e')
		self.network_label.bind('<Button-1>',lambda evt, whichframe=self.networkframe: 
                            self.frame_change(self.networkframe))

		updatebutton = tkinter.ttk.Button(self.sideframe, text='Update',width=10,command=self.updatesettings,
			style='black.TButton')
		updatebutton.pack(side='bottom',anchor='s',padx=10,pady=10)



		### PATHWAY MENU
		## SELECTION
		tkinter.ttk.Label(self.pathwayframe,text='Selection',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10)

		# K parameter
		kframe = tkinter.ttk.Frame(self.pathwayframe, style='white.TFrame')
		kframe.pack(fill='x',padx=(10,0))
		kparam_val = tk.IntVar()
		kparam_val.set(self.parent.kparam)
		tkinter.ttk.Label(kframe,text='Select levels of shortest pathways',style='white.TLabel').pack(anchor='w')		
		self.kparam = tk.Scale(kframe,from_=0,to=10,orient=tk.HORIZONTAL,background=white,foreground=black, showvalue=0, variable=kparam_val)
		self.kparam.pack(side='left',fill='x',expand=True)
		CreateToolTip(self.kparam, text='0 chooses the paths with the least number of reaction steps')
		tkinter.ttk.Label(kframe, textvariable=kparam_val, style='white.TLabel',width=2).pack(side='right')

		# Evaluate Reactions
		eframe = tkinter.ttk.Frame(self.pathwayframe,style='white.TFrame')
		eframe.pack(fill='x',padx=(10,0))
		tkinter.ttk.Label(eframe,text='Reaction types',style='white.TLabel').pack(anchor='w')
		self.evalrxns = tkinter.ttk.Combobox(eframe,values=['all','bio','chem'],state='readonly',style='white.TCombobox')
		self.evalrxns.pack(fill='x',expand=True)
		CreateToolTip(self.evalrxns, text='Types of reactions to include in the final pathways')
		self.evalrxns.set(self.parent.evalrxns)


		## ANALYSIS
		tkinter.ttk.Label(self.pathwayframe,text='Analysis',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=(30,10))

		# Tanimoto 
		self.tanimoto_check = tk.IntVar()
		tanimotobutton = tkinter.ttk.Checkbutton(self.pathwayframe,text='Tanimoto structural similarity analysis',
			variable=self.tanimoto_check,style='white.TCheckbutton')
		tanimotobutton.pack(fill='x',padx=(10,0))
		CreateToolTip(tanimotobutton, text='Find pathways for compounds structurally similar to the target (similarity assessed using the tanimoto index)')
		self.tanimoto_check.set(self.parent.tanimoto_check)

		tanimoto_frame = tkinter.ttk.Frame(self.pathwayframe, style='white.TFrame')
		tanimoto_frame.pack(fill='x',padx=(30,0))		
		tanthresh = tk.DoubleVar()	
		tanthresh.set('%.2f' % self.parent.tanimoto_thresh)
		tkinter.ttk.Label(tanimoto_frame, text='Tanimoto index threshold', style='white.TLabel').pack(anchor='w')		
		self.tanimoto_thresh = tk.Scale(tanimoto_frame,resolution=0.01,from_=0,to=1,orient=tk.HORIZONTAL,background=white,foreground=black, showvalue=0, variable=tanthresh)
		self.tanimoto_thresh.pack(side='left',fill='x',expand=True)
		tkinter.ttk.Label(tanimoto_frame, textvariable=tanthresh, style='white.TLabel', width=4).pack(side='right')

		# FBA
		self.fba_check = tk.IntVar()
		self.fbabutton = tkinter.ttk.Checkbutton(self.pathwayframe,text='Flux balance analysis',
			variable=self.fba_check,style='white.TCheckbutton')
		self.fbabutton.pack(fill='x',padx=(10,0))
		CreateToolTip(self.fbabutton, text='Identify maximum theoretical yield of the target using flux balance analysis')
		self.fba_check.set(self.parent.fba_check)

		# Knockout	
		self.ko_check = tk.IntVar()
		kobutton = tkinter.ttk.Checkbutton(self.pathwayframe,text='With knockouts',style='white.TCheckbutton',
			variable=self.ko_check,command=self.ko_change)
		kobutton.pack(fill='x',padx=(30,0))
		CreateToolTip(kobutton, text='Identify reaction knockouts which increase maximum theoretical yield of the target')
		self.ko_check.set(self.parent.ko_check)
		if self.ko_check.get():
			self.fbabutton.state(['disabled'])

		## GENE COMPATIBILITY
		tkinter.ttk.Label(self.pathwayframe,text='Gene Compatibility',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=(30,10))
		self.gene_compatibility_check = tk.IntVar()
		gcbutton = tkinter.ttk.Checkbutton(self.pathwayframe,text='Optimize gene compatibility',
			variable=self.gene_compatibility_check,style='white.TCheckbutton')
		gcbutton.pack(fill='x',padx=(10,0))
		CreateToolTip(gcbutton, text='Retrieves and optimizes gene sequences based on codon usage for enzymes necessary for target compound  production (requires a valid network connection / proxy)')
		self.gene_compatibility_check.set(self.parent.gene_compatibility_check)

		cai_frame = tkinter.ttk.Frame(self.pathwayframe, style='white.TFrame')
		cai_frame.pack(fill='x',padx=(30,0))		
		caithresh = tk.DoubleVar()	
		caithresh.set('%.2f' % self.parent.cai_optimal_threshold)
		tkinter.ttk.Label(cai_frame, text='Codon Adaptation Index - Optimal Threshold', style='white.TLabel').pack(anchor='w')		
		self.cai_optimal_threshold = tk.Scale(cai_frame,resolution=0.01,from_=0.38,to=0.85,orient=tk.HORIZONTAL,background=white,foreground=black, showvalue=0, variable=caithresh)
		self.cai_optimal_threshold.pack(side='left',fill='x',expand=True)
		tkinter.ttk.Label(cai_frame, textvariable=caithresh, style='white.TLabel', width=4).pack(side='right')



		### OUTPUT MENU
		## FILE
		tkinter.ttk.Label(self.outputframe,text='File',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10)

		# Output directory
		self.output_path = tk.StringVar()
		self.output_path.set(self.parent.output_path.get())
		self.output_check = tk.IntVar()
		outputbutton = tkinter.ttk.Checkbutton(self.outputframe, text='Output to desktop (default)',
			variable=self.output_check,style='white.TCheckbutton',command=self.output_change)
		outputbutton.pack(fill='x',padx=(10,0))		
		self.output_check.set(self.parent.output_check)

		opframe = tkinter.ttk.Frame(self.outputframe, style='white.TFrame')
		opframe.pack(fill='x',padx=(30,0))		
		output = tkinter.ttk.Entry(opframe, textvariable=self.output_path, justify='left', style='white.TEntry', state='disabled')
		output.pack(side='left',fill='both', expand=True)
		self.opbutton = tkinter.ttk.Button(opframe,text='Browse',width=7,command=self.get_directory,
			style='black.TButton')		
		self.opbutton.pack(side='right')
		if self.output_check.get():
			self.opbutton.state(['disabled'])

		# XLSX Format
		self.oxf_check = tk.IntVar()
		tkinter.ttk.Checkbutton(self.outputframe,text='Output as .xlsx (instead of .txt)', variable=self.oxf_check,
			style='white.TCheckbutton').pack(fill='x',padx=(10,0))
		self.oxf_check.set(self.parent.oxf_check)	



		## OUTPUT FIGURES
		tkinter.ttk.Label(self.outputframe,text='Figures',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=(30,10))

		# IUPAC Name
		self.iupac_name_check = tk.IntVar()
		iupac_name_button = tkinter.ttk.Checkbutton(self.outputframe, text='Use IUPAC names in figures (instead of InChI)', 
			style='white.TCheckbutton',	variable=self.iupac_name_check)
		iupac_name_button.pack(fill='x',padx=(10,0))
		CreateToolTip(iupac_name_button, text='IUPAC name retrieval from PubChem requires a valid network connection / proxy address')
		self.iupac_name_check.set(self.parent.iupac_name_check)

		# Chemical RXN Info
		self.chemical_rxn_info_check = tk.IntVar()
		rxn_info_button = tkinter.ttk.Checkbutton(self.outputframe, text='Display chemical reaction information in figures', 
			style='white.TCheckbutton',	variable=self.chemical_rxn_info_check)
		rxn_info_button.pack(fill='x',padx=(10,0))
		self.chemical_rxn_info_check.set(self.parent.chemical_rxn_info_check)







		### DATABASE MENU
		## EXISTING
		tkinter.ttk.Label(self.databaseframe,text='Existing database',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10)

		# Default DB
		self.default_db_check = tk.IntVar()
		self.default_db_button = tkinter.ttk.Checkbutton(self.databaseframe, text='Use the default database',
			variable=self.default_db_check,style='white.TCheckbutton')
		self.default_db_button.pack(fill='x',padx=(10,0))		
		self.default_db_button.bind('<Button-1>',lambda evt, whichdb=self.default_db_button: 
                            self.db_change(self.default_db_button))
		self.default_db_check.set(self.parent.default_db_check)
		

		# Other DB
		self.existing_db_check = tk.IntVar()
		self.existing_db_button = tkinter.ttk.Checkbutton(self.databaseframe, text='Use another existing database',
			variable=self.existing_db_check,style='white.TCheckbutton')
		self.existing_db_button.pack(fill='x',padx=(10,0))
		self.existing_db_button.bind('<Button-1>',lambda evt, whichdb=self.existing_db_button: 
                            self.db_change(self.existing_db_button))
		self.existing_db_check.set(self.parent.existing_db_check)

		self.database_file = tk.StringVar()
		dbfileframe = tkinter.ttk.Frame(self.databaseframe, style='white.TFrame')
		dbfileframe.pack(fill='x',padx=(30,0))		
		dbfile = tkinter.ttk.Entry(dbfileframe, textvariable=self.database_file, justify='left', style='white.TEntry', state='disabled')
		dbfile.pack(side='left',fill='both', expand=True)
		self.database_file.set(self.parent.database_file)
		self.dbfilebutton = tkinter.ttk.Button(dbfileframe,text='Browse',width=7,command=self.get_database_file,
			style='black.TButton')		
		self.dbfilebutton.pack(side='right')
		if not self.existing_db_check.get():
			self.dbfilebutton.state(['disabled'])

		## GENERATE DB
		tkinter.ttk.Label(self.databaseframe,text='New database',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=(30,10))

		self.generate_db_check = tk.IntVar()
		self.generate_db_button = tkinter.ttk.Checkbutton(self.databaseframe, text='Generate a new database',
			variable=self.generate_db_check,style='white.TCheckbutton')
		self.generate_db_button.pack(fill='x',padx=(10,0))
		self.generate_db_button.bind('<Button-1>',lambda evt, whichdb=self.generate_db_button: 
                            self.db_change(self.generate_db_button))
		self.generate_db_check.set(self.parent.generate_db_check)

		# DB Name for -gdb and -gdbc
		tkinter.ttk.Label(self.databaseframe,text='Database name',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10,padx=(10,0))
		gdbframe = tkinter.ttk.Frame(self.databaseframe,style='white.TFrame')
		gdbframe.pack(fill='x',padx=(30,0))
		self.gdbname = tkinter.ttk.Entry(gdbframe, justify='left',style='white.TEntry')
		self.gdbname.pack(fill='both',expand=True)
		self.gdbname.delete(0,tk.END)
		self.gdbname.insert(0,self.parent.gdbname)

		# InChI IDs
		tkinter.ttk.Label(self.databaseframe,text='InChI',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10,padx=(10,0))
		self.inchi_check = tk.IntVar()
		inchibutton = tkinter.ttk.Checkbutton(self.databaseframe,text='Use InChI identifiers',variable=self.inchi_check,
			style='white.TCheckbutton')
		inchibutton.pack(fill='x',padx=(30,0))
		self.inchi_check.set(self.parent.inchi_check)
		if self.tanimoto_check.get():
			inchibutton.state(['disabled'])
		CreateToolTip(inchibutton, text='InChI IDs are required for Tanimoto similarity analysis')

		# PATRIC DB
		tkinter.ttk.Label(self.databaseframe,text='PATRIC',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10,padx=(10,0))
		self.patric_check = tk.IntVar()
		tkinter.ttk.Checkbutton(self.databaseframe,text='Build database using PATRIC data',variable=self.patric_check,
			style='white.TCheckbutton').pack(fill='x',padx=(30,0))
		self.patric_check.set(self.parent.patric_check)

		patricsubframe = tkinter.ttk.Frame(self.databaseframe,style='white.TFrame')
		patricsubframe.pack(fill='x',padx=(50,0))

		# Username and password
		unpwframe = tkinter.ttk.Frame(patricsubframe,style='white.TFrame')
		unpwframe.pack(fill='x')		
		self.patric_un = tkinter.ttk.Entry(unpwframe,style='white.TEntry')
		self.patric_un.pack(side='left',fill='x',expand=True,ipady=5)
		self.patric_pw = tkinter.ttk.Entry(unpwframe,style='white.TEntry',show='*')
		self.patric_pw.pack(side='right',fill='x',expand=True,ipady=5)
		labelframe = tkinter.ttk.Frame(patricsubframe,style='white.TFrame')
		labelframe.pack(fill='x')
		tkinter.ttk.Label(labelframe,text='Username',font='Verdana 10',style='white.TLabel').pack(side='left',fill='x')
		tkinter.ttk.Label(labelframe,text='Password',font='Verdana 10',style='white.TLabel').pack(side='right',fill='x')
		self.patric_un.delete(0,tk.END)
		self.patric_pw.delete(0,tk.END)
		self.patric_un.insert(0,self.parent.patric_un)
		self.patric_pw.insert(0,self.parent.patric_pw)

		# Media Type
		pmediaframe = tkinter.ttk.Frame(patricsubframe,style='white.TFrame')
		pmediaframe.pack(fill='x')
		tkinter.ttk.Label(pmediaframe,text='PATRIC media type',style='white.TLabel').pack(anchor='w')
		self.patric_media = tkinter.ttk.Combobox(pmediaframe,values=['Carbon-D-Glucose','LB','Complete'],state='readonly',style='white.TCombobox')
		self.patric_media.pack(fill='x',expand=True,ipady=5)
		self.patric_media.set(self.parent.patric_media)


		patricfileframe = tkinter.ttk.Frame(patricsubframe,style='white.TFrame')
		patricfileframe.pack(fill='x')
		tkinter.ttk.Label(patricfileframe,text='Select PATRIC file',style='white.TLabel').pack(anchor='w')
		self.patric_file = tk.StringVar()
		patricfile = tkinter.ttk.Entry(patricfileframe, textvariable=self.patric_file, justify='left',style='white.TEntry', state="disabled")
		patricfile.pack(side='left',fill='both',expand=True, ipady=5)
		tkinter.ttk.Button(patricfileframe, text='Browse',width=7,command=self.get_patric_file,
			style='black.TButton').pack(side='right')
		self.patric_file.set(self.parent.patric_file)

		# METACYC DB
		tkinter.ttk.Label(self.databaseframe,text='MetaCyc',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10,padx=(10,0))
		self.metacyc_check = tk.IntVar()
		tkinter.ttk.Checkbutton(self.databaseframe,text='Build database using MetaCyc data',variable=self.metacyc_check,
			style='white.TCheckbutton').pack(fill='x',padx=(30,0))
		self.metacyc_check.set(self.parent.metacyc_check)

		metacycsubframe = tkinter.ttk.Frame(self.databaseframe,style='white.TFrame')
		metacycsubframe.pack(fill='x',padx=(50,0))

		metacycfileframe = tkinter.ttk.Frame(metacycsubframe,style='white.TFrame')
		metacycfileframe.pack(fill='x')
		metacyc_label = tkinter.ttk.Label(metacycfileframe,text='Select MetaCyc File',style='white.TLabel')
		metacyc_label.pack(anchor='w')
		CreateToolTip(metacyc_label,text='File can be downloaded from the MetaCyc website')
		self.metacyc_addition = tk.StringVar()
		metacycaddition = tkinter.ttk.Entry(metacycfileframe, textvariable=self.metacyc_addition, justify='left',style='white.TEntry', state="disabled")
		metacycaddition.pack(side='left',fill='x',expand=True, ipady=5)
		tkinter.ttk.Button(metacycfileframe, text='Browse',width=7,command=self.get_metacyc_addition,
			style='black.TButton').pack(side='right')
		self.metacyc_addition.set(self.parent.metacyc_addition)

		# KEGG DB
		tkinter.ttk.Label(self.databaseframe,text='KEGG',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10,padx=(10,0))
		self.kegg_check = tk.IntVar()
		tkinter.ttk.Checkbutton(self.databaseframe,text='Build database using KEGG data',variable=self.kegg_check,
			style='white.TCheckbutton').pack(fill='x',padx=(30,0))
		self.kegg_check.set(self.parent.kegg_check)

		keggsubframe = tkinter.ttk.Frame(self.databaseframe,style='white.TFrame')
		keggsubframe.pack(fill='x',padx=(50,0))

		# Organism Type
		korgframe = tkinter.ttk.Frame(keggsubframe,style='white.TFrame')
		korgframe.pack(fill='x')
		tkinter.ttk.Label(korgframe,text='KEGG organism type:',style='white.TLabel').pack(anchor='w')
		self.kegg_org_type = tkinter.ttk.Combobox(korgframe,values=['bacteria','algae','plants','fungi','Eukaryotes','prokaryotes','all'],
			state='readonly',style='white.TCombobox')
		self.kegg_org_type.pack(fill='x',expand=True,ipady=5)
		self.kegg_org_type.set(self.parent.kegg_org_type)

		# Number of organisms
		knunframe = tkinter.ttk.Frame(keggsubframe, style='white.TFrame')
		knunframe.pack(fill='x')		
		knun = tk.IntVar()	
		knun.set(self.parent.kegg_nunorgs)
		tkinter.ttk.Label(knunframe, text='Maximum # of organisms', style='white.TLabel').pack(anchor='w')		
		self.kegg_nunorgs = tk.Scale(knunframe,from_=0,to=100,orient=tk.HORIZONTAL,background=white,foreground=black, showvalue=0, variable=knun)
		self.kegg_nunorgs.pack(side='left',fill='x',expand=True)
		CreateToolTip(self.kegg_nunorgs,text='0 organisms will be considered as "all"')
		tkinter.ttk.Label(knunframe, textvariable=knun, style='white.TLabel', width=3).pack(side='right')

		# Number of organism paths
		knunpframe = tkinter.ttk.Frame(keggsubframe, style='white.TFrame')
		knunpframe.pack(fill='x')		
		knunp = tk.IntVar()	
		knunp.set(self.parent.kegg_nunorg_paths)
		tkinter.ttk.Label(knunpframe, text='Maximum # of organism paths', style='white.TLabel').pack(anchor='w')		
		self.kegg_nunorg_paths = tk.Scale(knunpframe,from_=0,to=100,orient=tk.HORIZONTAL,background=white,foreground=black, showvalue=0, variable=knunp)
		self.kegg_nunorg_paths.pack(side='left',fill='x',expand=True)
		CreateToolTip(self.kegg_nunorg_paths,text='0 organism paths will be considered as "all"')
		tkinter.ttk.Label(knunpframe, textvariable=knunp, style='white.TLabel', width=3).pack(side='right')
		
		



		### NETWORK MENU
		## PROXY
		tkinter.ttk.Label(self.networkframe, text='Proxy',font='Verdana 14 bold',style='white.TLabel').pack(anchor='w',pady=10)

		# PROXY SETTINGS
		proxyframe = tkinter.ttk.Frame(self.networkframe,style='white.TFrame')
		proxyframe.pack(fill='x')
		
		tmplabelframe = tkinter.ttk.Frame(proxyframe,style='white.TFrame')
		tmplabelframe.pack(side='left')
		tkinter.ttk.Label(tmplabelframe,text='Network Proxy',style='white.TLabel').pack()
		tkinter.ttk.Label(tmplabelframe,text='<ADDRESS>:<PORT>',style='white.TLabel',font='Verdana 10').pack()
		self.proxyaddress = tk.StringVar()
		self.proxy = tkinter.ttk.Entry(proxyframe,style='white.TEntry',justify='left',textvariable=self.proxyaddress)
		self.proxyaddress.set(self.parent.proxy)
		self.proxy.pack(side='right',fill='x',expand=True,ipady=5)

	# Helper functions		
	def frame_change(self,whichframe):
		frame_label_dict = {self.pathwayframe:self.pathway_label,
							self.outputframe:self.output_label,
							self.databaseframe:self.database_label,
							self.networkframe:self.network_label}
		for frame,label in frame_label_dict.items():
			if frame != whichframe:
				frame.pack_forget()
				label.config(style = 'green.TLabel')
		
		whichframe.pack(padx=20,pady=20,fill='both',expand=True)
		frame_label_dict[whichframe].config(style = 'black.TLabel')

	def ko_change(self):
		if self.ko_check.get():
			self.fba_check.set(1)
			self.fbabutton.state(['disabled'])
		else:
			self.fbabutton.state(['!disabled'])

	def output_change(self):
		if self.output_check.get():
			self.opbutton.state(['disabled'])
		else:
			self.opbutton.state(['!disabled'])

	def get_directory(self):
		if self.parent.output_path.get():
			initialdir = self.parent.output_path.get()
		else:
			initialdir = HOMEDIR
		directory = tkinter.filedialog.askdirectory(title="Select Output Directory",initialdir=initialdir)
		self.parent.output_path.set(directory)

	def db_change(self,whichdb):
		btn_check_dict = {self.default_db_button:self.default_db_check,
						self.existing_db_button:self.existing_db_check,
						self.generate_db_button:self.generate_db_check}		
		for btn,check in btn_check_dict.items():
			if btn != whichdb:
				check.set(0)
				if btn == self.existing_db_button:
					self.dbfilebutton.state(['disabled'])
		
		if whichdb == self.existing_db_button:
			self.dbfilebutton.state(['!disabled'])

	def get_database_file(self):
		if self.database_file.get():
			initialdir = os.path.dirname(self.database_file.get())
		else:
			initialdir = HOMEDIR		
		file = tkinter.filedialog.askopenfilename(title="Select Database File",initialdir=HOMEDIR,filetypes=[("Database files","*.db")])
		self.database_file.set(file)

	def get_patric_file(self):
		if self.patric_file.get():
			initialdir = os.path.dirname(self.patric_file.get())
		else:
			initialdir = HOMEDIR
		file = tkinter.filedialog.askopenfilename(title="Select PATRIC File",initialdir=initialdir,filetypes=[('CSV Files','*.csv')])
		self.patric_file.set(file)

	def get_metacyc_addition(self):
		if self.metacyc_addition.get():
			initialdir = os.path.dirname(self.metacyc_addition.get())
		else:
			initialdir = HOMEDIR
		file = tkinter.filedialog.askopenfilename(title="Select MetaCyc File",initialdir=initialdir,filetypes=[("XML Files","*.xml")])
		self.metacyc_addition.set(file)

	# Update button
	def updatesettings(self):
		self.parent.kparam = self.kparam.get()
		self.parent.evalrxns = self.evalrxns.get()
		self.parent.oxf_check = self.oxf_check.get()
		self.parent.tanimoto_check = self.tanimoto_check.get()
		if self.tanimoto_check.get():
			self.parent.inchi_check = 1
		self.parent.tanimoto_thresh = self.tanimoto_thresh.get()			
		self.parent.fba_check = self.fba_check.get()
		self.parent.ko_check = self.ko_check.get()
		self.parent.iupac_name_check = self.iupac_name_check.get()
		self.parent.chemical_rxn_info_check = self.chemical_rxn_info_check.get()
		self.parent.gene_compatibility_check = self.gene_compatibility_check.get()
		self.parent.cai_optimal_threshold = self.cai_optimal_threshold.get()

		# Reset parent database files
		self.parent.database_file = ''
		self.parent.gdbname = ''

		# Generate a new database
		if self.generate_db_check.get():
			if not self.patric_check.get() and not self.metacyc_check.get() and not self.kegg_check.get():
				ErrorWindow(parent_window=self.window, message='Database generation requires at least one of PATRIC, MetaCyc, or KEGG.')
				return			
			
			# Double check database name
			if self.gdbname.get() == '':
				ErrorWindow(parent_window=self.window, message='Please name your custom database to be generated.')
				return
			else:
				self.parent.gdbname = self.gdbname.get()

			# InChI
			self.parent.inchi_check = self.inchi_check.get()
			
			# PATRIC
			self.parent.patric_check = self.patric_check.get()
			self.parent.patric_media = self.patric_media.get()
			self.parent.patric_file = self.patric_file.get()
			if self.patric_check.get():
				# Authentication
				try:
					token = mackinac.get_token(username=self.patric_un.get(), password=self.patric_pw.get())
					self.parent.patric_un = self.patric_un.get()
					self.parent.patric_pw = self.patric_pw.get()
				except:
					ErrorWindow(parent_window=self.window, message='PATRIC Authentication Error\nPlease double check your username, password, and proxy settings.')
					return			

			# MetaCyc
			if self.metacyc_check.get() and not self.metacyc_addition.get():
				ErrorWindow(parent_window=self.window, message='MetaCyc requires a corresponding MetaCyc XML file.\nThis can be downloaded from the MetaCyc website.')
				return
			else:
				self.parent.metacyc_check = self.metacyc_check.get()
				self.parent.metacyc_addition = self.metacyc_addition.get()

			# KEGG
			self.parent.kegg_check = self.kegg_check.get()
			self.parent.kegg_org_type = self.kegg_org_type.get()
			self.parent.kegg_nunorgs = self.kegg_nunorgs.get()
			self.parent.kegg_nunorg_paths = self.kegg_nunorg_paths.get()
			
		# Select an existing database
		elif self.existing_db_check.get():
			self.parent.database_file = self.database_file.get()			

		# Use default database
		else:			
			self.parent.default_db_check = 1
			self.parent.database_file = self.parent.default_database

		self.parent.default_db_check = self.default_db_check.get()
		self.parent.existing_db_check = self.existing_db_check.get()
		self.parent.generate_db_check = self.generate_db_check.get()

		self.update_parent_database()
		
		# Update proxy
		if self.proxyaddress.get():
			self.parent.proxy = self.proxyaddress.get()
			os.environ["PROXY"] = self.proxyaddress.get()
			os.environ["HTTP_PROXY"] = self.proxyaddress.get()
			os.environ["HTTPS_PROXY"] = self.proxyaddress.get()

		self.window.destroy()

	def update_parent_database(self):
		# Generate database
		if self.generate_db_check.get():
			self.parent.target.set('Custom target -- ')			
			if self.inchi_check.get():
				self.parent.target_options = ['InChI']
			else:
				self.parent.target_options = ['KEGG']
			self.parent.target_label.set('Target ID (%s):' % '/'.join(self.parent.target_options))

			if self.patric_file.get():
				patricfile = self.patric_file.get()
			else:
				patricfile = DEFAULT_PATRICFILE

			f = open(patricfile,'r')
			models = []
			f.readline() # Sacrificial readline - header row
			line = f.readline()
			while line:
				line = line.split(',')
				genome_id = line[0]
				genome_name = line[1]
				models += ['%s -- %s' % (str(genome_id),'_'.join(genome_name.split(' ')) + '_' + self.patric_media.get())]
				line = f.readline()
			f.close()
			self.parent.start['values'] = models

		# Select database
		elif self.existing_db_check.get():
			from Database import query as Q
			DB = Q.Connector(self.database_file.get())

			self.parent.target_options = ['Name','KEGG']
			if DB.custom_query("SELECT * FROM compound WHERE ID LIKE 'InChI%'") != 'None':
				self.parent.target_options += ['InChI','CAS','Chemical formula']
			self.parent.target_label.set('Target (%s):' % '/'.join(self.parent.target_options))
			
			models = ['%s -- %s' % (m[0],m[1]) for m in DB.get_all_fba_models()]
			self.parent.start['values'] = models

		# Default database
		else:
			self.parent.target_options = ['Name','Chemical formula','InChI','KEGG','CAS']
			self.parent.target_label.set('Target (%s):' % '/'.join(self.parent.target_options))
			self.parent.start['values'] = self.parent.default_models
			

class RunWindow(Window):
	def __init__(self,request_form,width=900,height=400,background=black):
		self.rqform = request_form
		super(RunWindow,self).__init__(width=width,height=height,background=background)
		
	def initialize(self):		
		self.autoscroll = True
		self.window.title('Retrieving pathways for %s' % self.rqform["target"])
		stframe = tkinter.ttk.Frame(self.window)
		stframe.pack(fill='both',expand=True)
		self.progresstext = tkst(stframe, bg=black, fg=green)
		self.progresstext.pack(fill='both',expand=True)
		self.progresstext.bind('<ButtonPress>',self.scroll_off)
		self.progresstext.bind('<ButtonRelease>',self.scroll_on)

		stdout_q = StdoutQueue()
		self.p = multiprocessing.Process(target=run_retsynth, args=(self.rqform,stdout_q,))
		monitor = Thread(target=self.text_catcher,args=(self.progresstext,stdout_q))
		monitor.daemon = True
		monitor.start()
		self.window.protocol("WM_DELETE_WINDOW", self.on_closing)		
		self.p.start()

	def scroll_on(self,event):
		self.autoscroll = True

	def scroll_off(self,event):
		self.autoscroll = False

	def text_catcher(self,text_widget,queue):
		while True:
			text_widget.insert(tk.END, queue.get())
			if self.autoscroll:
				text_widget.see(tk.END)
			if self.p.exitcode is not None:
				if self.p.exitcode != 0:
					exitmsg = '\nERROR:\tRetSynth failed with exit code: %s' % self.p.exitcode
					print(exitmsg,file=sys.stderr)
					text_widget.insert(tk.END, exitmsg)
					text_widget.see(tk.END)
				break
						

	def on_closing(self):
		if self.p.is_alive():
			if askokcancel("RetSynth is still running...","Are you sure you want to quit running RetSynth?"):
				self.p.terminate()
				self.window.destroy()
		else:
			self.window.destroy()