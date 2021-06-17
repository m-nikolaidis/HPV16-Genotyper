from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import logging
import pathlib
import pandas as pd
import os
import sys

"""Write in md:
oufmt columns can be altered based on user needs directly from the excel params file.
The addition of the text should be in the same cell and separated by spaces
The raw blast_results are automatically turned into an excel file for ease of interpretation (The txt file is also kept)
The script will read either the csv or the excel file, whichever is available. If none is, it will throw an error
TODO: Add list of all possibilities for non-cli users

Add something about the workflow on the MD file?
"""
def blastn_search(exe_params, workflow):
	"""
	Execute blast search
	Input query is ref vector fasta file, input database is genomes fasta for genotyping
	The params dict will be used for the parameters of blast
	"""
	threads = exe_params['num_threads']
	input_dir = exe_params['in']
	outdir = exe_params['out']
	query_f = exe_params['query']
	evalue = exe_params['evalue']
	# db_file = exe_params['blastdb'] # TODO: Remove this and from the params, or add two different fields, so the user can specify where they are located?
	# Better idea. Automatically generate the script_dir/../resources as input for the two databases, if the user wants can change it
	word_size = exe_params['word_size']
	blastres_f = outdir + query_f + "_Blastn_res.xml"
	system = exe_params['system']
	# system = sys.platform #TODO: remove after testing
	if  system == 'linux':
		makeblastdb_bin = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/Linux_bin/makeblastdb'))).replace(os.path.basename(__file__) + "/", "")
		blastn = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/Linux_bin/blastn'))).replace(os.path.basename(__file__) + "/", "")
	if system == 'win32' or system == 'win64':
		raise Exception("Windows BLAST not implemented yet")

	if workflow == "snp":
		db_path = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/NC_001526_probes.fa'))).replace(os.path.basename(__file__) + "/", "")
		db_file = os.path.basename(db_path)
		evalue = 0.005 # Fixed for SNP detection, based on previous analysis
		word_size = 4
	if workflow == "gene identification":
		db_path = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/16refs_gene_db.fa'))).replace(os.path.basename(__file__) + "/", "")
		db_file = os.path.basename(db_path)
	makedb = NcbimakeblastdbCommandline(cmd = makeblastdb_bin,
										dbtype ="nucl",
										input_file = input_dir + db_path,
										out = outdir + "/tmpBlastn_dir/" + db_file,
										)
	blastn_search = NcbiblastnCommandline(cmd = blastn,
											query = input_dir + query_f,
											db = outdir + "/tmpBlastn_dir/" + db_file,
											outfmt = 5,
											max_target_seqs = 100,
											out = blastres_f,
											evalue = evalue,
											num_threads = threads,
											word_size = word_size,
											dust="no"
											)
	makedb()
	# logging.debug(" Starting BLASTn search ")
	blastn_search()
	# logging.debug(" Finished BLASTn search ")

params_f = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../input/params.xlsx'))).replace(os.path.basename(__file__) + "/", "")
params_f_windows = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../input/params.xlsx'))).replace(os.path.basename(__file__) + "\\", "")
# If i have windows the path_cut should be "\\"
# In linux it is "/"
# Check if there is built in functionality to add them without specifying
paramsdf = pd.read_excel(params_f,index_col=0, engine="openpyxl") # TODO, add this from main script

params = paramsdf.to_dict()['Value']
workflow = "gene identification"
blastn_search(params, workflow)