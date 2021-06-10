from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import logging
import pathlib
import pandas as pd
import numpy as np
import os
"""Write in md:
oufmt columns can be altered based on user needs directly from the excel params file.
The addition of the text should be in the same cell and separated by spaces
The raw blast_results are automatically turned into an excel file for ease of interpretation (The txt file is also kept)
The script will read either the csv or the excel file, whichever is available. If none is, it will throw an error
TODO: Add list of all possibilities for non-cli users
"""

def blastn_search(exe_params):
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
	db_file = exe_params['blastdb']
	word_size = exe_params['word_size']
	blast_outfmt = exe_params['outfmt']
	blastres_f = outdir + query_f + "_Blastn_res.xml"
	if exe_params['system'] == 'linux':
		makeblastdb_bin = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/Linux_bin/makeblastdb'))).replace(os.path.basename(__file__) + "/", "")
		blastn = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/Linux_bin/blastn'))).replace(os.path.basename(__file__) + "/", "")
	if exe_params['system'] == 'win32' or exe_params['system'] == 'win64':
		print("Windows BLAST not implemented yet")

	makedb = NcbimakeblastdbCommandline(cmd = makeblastdb_bin,
										dbtype ="nucl",
										input_file = input_dir + db_file,
										out = outdir + "/tmpBlastn_dir/" + db_file,
										)
	blastn_search = NcbiblastnCommandline(cmd = blastn,
											query = input_dir + query_f,
											db = outdir + "/tmpBlastn_dir/" + db_file,
											outfmt = blast_outfmt,
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


paramsdf = pd.read_excel("../input/params.xlsx",index_col=0)
params = paramsdf.to_dict()['Value']
blastn_search(params)