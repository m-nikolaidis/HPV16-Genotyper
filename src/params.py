import os
import sys
import pathlib
import multiprocessing
import pandas as pd

def default_param_file(outdir):
	"""Intialize the default parameters and write them to the param.xlsx file
	The user should have the option to output a reference param file and change it accordingly
	"""
	
	def _defaultparams():
		system = sys.platform
		input_dir = pathlib.Path(os.getcwd()) / pathlib.Path("../input") # The user should provide the input dir cmd? just like the orthologues pipeline
		query_f = pathlib.Path(os.getcwd()) / pathlib.Path("../input/Input.fasta") # The user should provide the input dir cmd? just like the orthologues pipeline
		output_dir = pathlib.Path(os.getcwd()) / pathlib.Path("../Script_out")
		threads_to_use =  multiprocessing.cpu_count() - 2
		params = {
			"system": system,
			"in": input_dir,
			"out": output_dir,
			"num_threads": threads_to_use,
			"query":query_f,
			"SNP_identification_Evalue":0.005,
			"SNP_identification_Word_Size":4,
			"SNP_annotation_file":None, # TODO Implement the annotation file
			"SNP_db_path":os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/NC_001526_probes.fa'))).replace(os.path.basename(__file__) + "/", ""),
			"Gene_identification_Evalue": 1e-5, 
			"Gene_identification_Word_size": 7,
			"GenesProfile_directory":os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/profiles'))).replace(os.path.basename(__file__) + "/", ""),
			"GenesRef_database":os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/16refs_gene_db.fa'))).replace(os.path.basename(__file__) + "/", "")		

		}
		return params

	params_f = pathlib.Path(outdir) / pathlib.Path("params.xlsx")
	params = _defaultparams()
	param_df = pd.DataFrame.from_dict(params,orient='index')
	param_df.rename(columns={0:"Value"},inplace=True)
	param_df.to_excel(params_f)