import pandas as pd
from Bio import SeqIO
import re
from Bio.Align.Applications import MuscleCommandline
import numpy as np
import os
import pathlib
import sys

def prepare_gene_files_for_aln(params):
	"""
	Write a better docstring
	"""
	def _cut_genes(genomes_f, coords_df):
		org_regex = re.compile(r'^(\S+)')
		final_seqs = {}
		seqs = SeqIO.to_dict(SeqIO.parse(genomes_f, "fasta"))
		orgs = list(seqs.keys())
		for org in orgs:
			tmpdf = coords_df[coords_df["qaccver"] == org]
			indeces = tmpdf.index
			for idx in indeces:
				gene = tmpdf.loc[idx,"saccver"][-2:]
				if gene not in final_seqs:
					final_seqs[gene] = {}
				if org not in final_seqs[gene]:
					final_seqs[gene][org] = ""
				start = tmpdf.loc[idx,"qstart"]
				end = tmpdf.loc[idx,"qend"]
				seq = seqs[org].seq[int(start)-1:int(end)]
				final_seqs[gene][org] = str(seq)
		return final_seqs
		

	output_dir = params['out']
	# os.mkdir(output_dir + "tmp_dir/Batch_genes")
	db_path = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/16refs_gene_db.fa'))).replace(os.path.basename(__file__) + "/", "")
	dbseqs = SeqIO.to_dict(SeqIO.parse(db_path, "fasta"))
	
	query_f = params['query']
	coord_excel = output_dir + query_f + "Blastnres.xlsx"
	coords_df = pd.read_excel(coord_excel,engine="openpyxl")

	# https://pandas.pydata.org/pandas-docs/stable/user_guide/scale.html
	# Check with large dataset, to be more memory efficient

	batch_fa_dir = output_dir + "tmp_dir/" + "Batch_fasta_files/"
	batch_fa_files = os.listdir(batch_fa_dir)
	db_rec = SeqIO.to_dict(SeqIO.parse(db_path, "fasta"))
	db_orgs = [g+str(i) for g in ["A","B","C","D"] for i in range(1,5)]
	for batch_f in batch_fa_files:
		batchname = batch_f.split(".")[0]
		gene_seqs = _cut_genes(batch_fa_dir + batch_f, coords_df)
		for gene in gene_seqs:
			fout = output_dir + "tmp_dir/Batch_genes/" + batchname + "_" + gene + ".fa"   # TODO: use pathlib
			fhandle = open(fout,"w")
			for k in gene_seqs[gene]:
				str_to_write = ">" + k + "\n" + str(gene_seqs[gene][k]) + "\n"
				fhandle.write(str_to_write)
			for db_org in db_orgs:
				str_to_write = ">" + db_org + "\n" + str(db_rec[db_org + "_" + gene].seq) + "\n"
				fhandle.write(str_to_write)
			fhandle.close()
	#seqpath = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path(output_dir + 'tmp_dir/.fa'))).replace(os.path.basename(__file__) + "/", "")


def align_sequences(params,system, method="muscle"):
	output_dir = pathlib.Path(params['out'])
	# output_dir = "C:\\Users\\mario\\Desktop\\Lab\\HPV16\\HPV16_genotyping_tool\\Script_out"
	batch_gene_aln_dir = os.path.join(str(output_dir),str(pathlib.Path("tmp_dir/Alignments/")))
	# os.mkdir(batch_gene_aln_dir)
	muscle_bin_unix = os.path.join(str(pathlib.Path(__file__)), 
		str(pathlib.Path('../resources/Linux_bin/muscle3.8.31_i86linux64'))).replace(os.path.basename(__file__) + "/", "")
	muscle_bin_win = os.path.join(str(pathlib.Path(__file__)), 
		str(pathlib.Path('../resources/Win_bin/muscle3.8.31.exe'))).replace(os.path.basename(__file__) + "\\", "")
	batch_gene_dir = os.path.join(str(output_dir),str(pathlib.Path("tmp_dir/Batch_genes/")))
	files = os.listdir(batch_gene_dir)
	aln_files = []
	for f in files:
		fin = str(os.path.join(batch_gene_dir,pathlib.Path(f)))
		fout =str(os.path.join(batch_gene_aln_dir,pathlib.Path(f[:-3]))) + ".phy"
		aln_files.append(fout)
		if method == "muscle":
			if "win" in system:
				cline = MuscleCommandline(cmd = muscle_bin_win, input=fin, out=fout, phyi=True )
			if "linux" in system:
				cline = MuscleCommandline(cmd = muscle_bin_unix, input=fin, out=fout, phyi=True )
			# The phyi parameter will truncate the Fasta headers to 10 character long strings
			# Rename them afterwards. Does not interfere with alignment if the names are in variable lengths
		# cline()
	return aln_files

def _rename_phylip():
	pass


def build_tree(params,aln_files, method="PhyML", dist="Kimura", bootstrap_repl=1000):
	output_dir = pathlib.Path(params['out'])
	batch_gene_trees_dir = os.path.join(str(output_dir),str(pathlib.Path("tmp_dir/Trees/")))
	# os.mkdir(batch_gene_trees_dir)
	def _clean_phyml_output(phylip_aln_f, new_dir):
		# Move _phyml_tree.txt files to the appropriate directory and remove _phyml_stats.txt
		basename = os.path.basename(phylip_aln_f)
		phyml_tree_f = phylip_aln_f + "_phyml_tree.txt"
		phyml_tree_f_renamed = os.path.join(str(new_dir),str(pathlib.Path(basename + "_phyml_tree.txt")))
		phyml_stats_f = phylip_aln_f + "_phyml_stats.txt"
		try:
			pathlib.Path(phyml_tree_f).rename(phyml_tree_f_renamed)
		except:
			pass
		try:
			pathlib.Path(phyml_stats_f).unlink()
		except:
			pass

	seaview_lin_bin = os.path.join(str(pathlib.Path(__file__)), 
		str(pathlib.Path('../resources/Linux_bin/seaview5/seaview'))).replace(os.path.basename(__file__) + "/", "")
	phyml_lin_bin = os.path.join(str(pathlib.Path(__file__)), 
		str(pathlib.Path('../resources/Linux_bin/phyml'))).replace(os.path.basename(__file__) + "/", "")
	for f_phylip in aln_files:
		basename = os.path.basename(f_phylip)
		tree_file_out = os.path.join(str(batch_gene_trees_dir),
			str(pathlib.Path("tmp_dir/Trees/" + basename + "_NJ_tree.nwk")))
		if method == "BioNJ":
			cmd = seaview_lin_bin + " -build_tree -NJ -distance " \
			+ dist + " -replicates " + str(bootstrap_repl) + " -o " + tree_file_out + " " + f_phylip
			os.system(cmd)
			# Requires ubuntu 20 to run with latest seaview. 
			# (/lib/x86_64-linux-gnu/libm.so.6: version `GLIBC_2.29' not found)
			# Write in MD
			# Need to check if this works in Lubuntu vm, in WSL i get segmentation fault. I think is memory problem
			# It tries to acess forbidden memory adresses
		if method == "PhyML":
				cmd = phyml_lin_bin + " --quiet -o tl -s SPR -v estimated -m GTR -d nt -b -4 -f m -i " + f_phylip
				os.system(cmd)
				_clean_phyml_output(f_phylip,batch_gene_trees_dir)
				# TODO: should add check if files exist? If not skip





params_f_unix = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../input/params.xlsx'))).replace(os.path.basename(__file__) + "/", "")
params_f_win = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../input/params.xlsx'))).replace(os.path.basename(__file__) + "\\", "")
# If i have windows the path_cut should be "\\"
# In linux it is "/"
# Check if there is built in functionality to add them without specifying
# system = "Win32"
system = sys.platform
if "win" in system:
	paramsdf = pd.read_excel(params_f_win,index_col=0, engine="openpyxl") # TODO, add this from main script, check if correct
else:
	paramsdf = pd.read_excel(params_f_unix,index_col=0, engine="openpyxl") # TODO, add this from main script

params = paramsdf.to_dict()['Value']
# prepare_gene_files_for_aln(params)
aln_files = align_sequences(params,system)
build_tree(params, aln_files, method="BioNJ")