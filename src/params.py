import pandas as pd
import sys
import multiprocessing
import os

# Create param file
def create_default_param_file():
	# The user should have the option to output a reference param file and change it accordingly
	pass
	def _defaultparams():
		system = sys.platform
		input_dir = os.getcwd() + "/../input/" # The user should provide the input dir cmd? just like the orthologues pipeline
		output_dir = os.getcwd() + "/../Script_out/"
		threads_to_use =  multiprocessing.cpu_count() - 2
		# query_vectors = "test3vectors.fa" # TODO CHANGE THIS
		# blastn_db = "16ref_CG.fa" # TODO CHANGE THIS
		params = {
			"system": system,
			"in": input_dir,
			"out": output_dir,
			"evalue": 10, 
			"max_target_seqs": 10,
			"word_size": 4, # Minimum for blastn search
			"outfmt" : 5,
			"num_threads": threads_to_use,
			"query":query_vectors,
			"blastdb":blastn_db,
			"annot":None, # TODO Implement the annotation file
		}
	
		return params
	params = _defaultparams()
	param_df = pd.DataFrame.from_dict(params,orient='index')
	param_df.rename(columns={0:"Value"},inplace=True)
	param_df.to_excel(params['in'] + "params.xlsx")
	# This file is generated in the src directory. Should it be outputed in a different directory?
	# The user should choose the excel or csv output for this file
	""" To write in .md
	systems supported at the moment: Linux, WSL?
	The script can be installed in whatever folder that has user privileges
	input directory for files is script/../input by default
	Evalue notation accepts both scientific and floating point (due to pandas)
	Word_size should be int >= 4
	num_threads is system threads - 2 by default
	Currently using blast+ 2.11.0
	All the files should be put in the desired input directory, except for the params excel. which can be wherever
	The database and query should be only the file names. Not the paths
	TODO: the annotation should be given only as the file name, not the whole directory
	"""

def load_params(param_file):
	def _updateparams(param_file, params):
	
		param_csv = csv.reader(open(param_file, "r"), delimiter="\t")
		for lines in param_csv:
			key, value = lines
			params[key] = value
		
		if params["out"] == "": params["out"] = params["in"]
		
		# Gather the fasta files for the analysis
		fasta_files = os.listdir(params["in"] + "/Fasta_files/")
		params["fasta_files"] = fasta_files
		
		return params

	# These are returned to the main.py module
	# Implement error if no blast query and database names are not provided

create_default_param_file()