import os
import sys
import pathlib
import logging
import datetime
import argparse
import numpy as np
import pandas as pd
import multiprocessing
from Bio import SeqIO
# Tool modules
import annot
import blast
import env_setup
import phylogeny
import visualizations

def default_param_file(outdir):
	"""Intialize the default parameters and write them to the param.xlsx file
	The user should have the option to output a reference param file and change it accordingly
	"""
	
	def _defaultparams():
		system = sys.platform
		input_dir = pathlib.Path(os.getcwd()) / pathlib.Path("../input") # The user should provide the input dir cmd? just like the orthologues pipeline
		query_f = pathlib.Path(os.getcwd()) / pathlib.Path("../input/Input.fasta") # The user should provide the input dir cmd? just like the orthologues pipeline
		outdir = pathlib.Path(os.getcwd()) / pathlib.Path("../Script_out")
		threads_to_use =  multiprocessing.cpu_count() - 2
		params = {
			"system": system,
			"in": input_dir,
			"out": outdir,
			"num_threads": threads_to_use,
			"query": query_f,
			"SNP_identification_Evalue": 0.005,
			"SNP_identification_Word_Size": 4,
			"SNP_annotation_file":None, # TODO Implement the annotation file
			"SNP_db_path": pathlib.Path(os.getcwd()) / pathlib.Path('../resources/sequences/NC_001526_probes.fa'),
			"Gene_identification_Evalue": 1e-5, 
			"Gene_identification_Word_size": 7,
			"GenesProfile_directory": pathlib.Path(os.getcwd()) / pathlib.Path('../resources/sequences/profiles'),
			"GenesRef_database":pathlib.Path(os.getcwd()) / pathlib.Path('../resources/sequences/16refs_gene_db.fa')
		}
		return params
	if outdir.is_dir():
		params_f = outdir / pathlib.Path("params.xlsx")
	else:
		params_f = outdir
	params = _defaultparams()
	param_df = pd.DataFrame.from_dict(params,orient='index')
	param_df.rename(columns={0:"Value"},inplace=True)
	param_df.index.name = "Parameter"
	param_df.to_excel(params_f)

def main(paramsdf, exe=True):
	def _init_binaries(system):
		"""Initialize the paths for the needed binaries depending on the system
		input: System info
		return: List of pathlib paths
		"""
		if "win" in system:
			makeblastdb_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Win_bin/makeblastdb.exe'))).replace(os.path.basename(__file__) + "\\", "")
			blastn_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Win_bin/blastn.exe'))).replace(os.path.basename(__file__) + "\\", "")
			muscle_bin= os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Win_bin/muscle3.8.31.exe'))).replace(os.path.basename(__file__) + "\\", "")
			seaview_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Win_bin/seaview.exe'))).replace(os.path.basename(__file__) + "\\", "")
			phyml_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Win_bin/PhyML-3.1_win32.exe'))).replace(os.path.basename(__file__) + "\\", "")
			raise Warning("Windows BLAST fails for some reason when using biopython, but runs correctly from PowerShell")
			# TODO: Investigate the warning
		if "linux" in system:
			makeblastdb_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Linux_bin/makeblastdb'))).replace(os.path.basename(__file__) + "/", "")
			blastn_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Linux_bin/blastn'))).replace(os.path.basename(__file__) + "/", "")
			muscle_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Linux_bin/muscle3.8.31_i86linux'))).replace(os.path.basename(__file__) + "/", "")
			seaview_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Linux_bin/seaview4'))).replace(os.path.basename(__file__) + "/", "")
			phyml_bin = os.path.join(str(pathlib.Path(__file__)), 
				str(pathlib.Path('../resources/Linux_bin/PhyML-3.1_linux64'))).replace(os.path.basename(__file__) + "/", "")
		return [makeblastdb_bin,blastn_bin,muscle_bin,seaview_bin,phyml_bin]

	# Basic variables
	outdir = pathlib.Path(paramsdf.loc["out","Value"])
	env_setup.create_dirs(outdir, exist_ok=True)

	logfile = outdir / pathlib.Path("logfile.log")
	logger = logging.getLogger()
	fhandler = logging.FileHandler(filename=logfile)
	formatter = logging.Formatter("%(asctime)s\t%(message)s")
	fhandler.setFormatter(formatter)
	logger.addHandler(fhandler)
	logger.setLevel(logging.DEBUG)
	logging.info("Environment set successfully in dir %s" %(outdir))
	logging.info("Parameters used in this run: ")
	indeces = paramsdf.index
	for index in indeces:
		val = paramsdf.loc[index,"Value"]
		logging.info(F"param: {index},{val}")

	system = paramsdf.loc["system","Value"]
	indir = pathlib.Path(paramsdf.loc["in","Value"])
	query_f = pathlib.Path(paramsdf.loc["query","Value"])
	query_f_path = indir / query_f
	makeblastdb_bin,blastn_bin,muscle_bin,seaview_bin,phyml_bin = _init_binaries(system)
	annot_f = pathlib.Path(paramsdf.loc["SNP_annotation_file","Value"])
	
	# Environment setup
	env_setup.file_exists(query_f_path)
	env_setup.input_file_format(query_f_path)
	query_f_filt = env_setup.filter_input_file(query_f_path)
	query_f = query_f_filt # Posix | Windows Path
	query_f_path = indir / query_f
	env_setup.dir_priviledges(outdir)
	
	# BLASTn search
	workflow = "gene identification"
	geneid_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, exe=exe)
	geneid_blastn_res_path_xlsx = blast.parse_GeneID_results(geneid_blastn_res_path,paramsdf, exe=exe)
	workflow = "snp"
	snp_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, exe=exe)
	snp_blastn_res_path_xlsx = blast.parse_SNP_results(snp_blastn_res_path, exe=exe)
	

	# Gene alignments and trees
	aln_files = phylogeny.prepare_alns(outdir, query_f_path,geneid_blastn_res_path_xlsx, muscle_bin, exe=exe)
	profiledb_dir = paramsdf.loc["GenesProfile_directory","Value"]
	aln_files = phylogeny.profile_aln(outdir,aln_files,profiledb_dir, muscle_bin, exe=exe)
	trees_dir = phylogeny.build_trees(outdir, aln_files, phyml_bin, seaview_bin, exe=exe)
	
	# Clean tmp directory if user wants to save space before graphics
	# TODO: Implement
	
	geneid_res_df = pd.read_excel(geneid_blastn_res_path_xlsx, engine="openpyxl")
	color_dict = {
		"A":"#00ff00",
		"B":"#0080ff",
		"C":"#ff8000",
		"D":"#ff0000"
	}
	recombinants = blast.find_recombinants(geneid_res_df)
	visualizations.GeneID_fig(geneid_res_df, outdir, profiledb_dir, color_dict, recombinants)
	
	num_snps = annot.annotate_results(annot_f,snp_blastn_res_path_xlsx,exe=exe)
	
	snp_resultsdf = pd.read_excel(snp_blastn_res_path_xlsx,index_col=0, engine="openpyxl")
	snp_resultsdf = snp_resultsdf.drop("Nucleotides")
	snp_resultsdf = snp_resultsdf.T
	snp_resultsdf = snp_resultsdf.head(num_snps)
	snp_results_annot = outdir / "SNPs_results.xlsx"
	snp_resultsdf_annot = pd.read_excel(snp_results_annot,index_col=0, engine="openpyxl")
	snp_resultsdf_annot = snp_resultsdf_annot.T
	snp_resultsdf_annot = snp_resultsdf_annot.head(num_snps)
	color_dict = {
		"Lin_A":"#00ff00",
		"Lin_B":"#0080ff",
		"Lin_C":"#ff8000",
		"Lin_D":"#ff0000",
		"Lin_BCD":"#ebeb34",
		"Other":"#ffffff"
	}
	visualizations.SNP_fig(snp_resultsdf_annot, outdir, color_dict, snp_resultsdf, recombinants)

	visualizations.render_trees(trees_dir, outdir)

	return trees_dir, outdir

if __name__ == "__main__":
	# Arguments
	parser = argparse.ArgumentParser(description="HPV16 genotyping tool CLI")
	parser.add_argument('--p_in', metavar="F",type=str,help='Absolute or relative path to the existing parameters file')
	parser.add_argument('--p_out', metavar="F",type=str,help='Directory to write the default parameters file')
	parser.add_argument('--clean', metavar="bool",type=bool,help='Remove the tmp directories created during execution (True / False)')
	parser.add_argument('--vis', metavar="bool",type=bool,help='Visualize phylogenetic trees on ete3 library. If this is false the trees will be rendered as images')
	parser.add_argument('--tree_render', metavar="Format",type=str,help='Type of file to render the phylogenetic trees available types: jpg (Default), pdf, png')
	
	# TODO: Add a dual switch p_in and p_out. Either read or write parameters file, and make them positional
	args = parser.parse_args()
	
	params_o = args.p_out
	params_i = args.p_in
	params_clean = args.clean
	params_viz = args.vis
	params_tree_render = args.tree_render
	
	now = datetime.datetime.now()
	start_dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	
	if params_o:
		default_param_file(params_o)
		print("Created parameters file successfully, please rerun the script using the parameters file")
		sys.exit(0) # Kill the CLI when outputing parameters file

	create_annot = False
	if create_annot:
		fin = "../input/signature_nucl_data_20210602_for_params.xlsx" # TODO: Fix to Take it from params
		annot_f = pathlib.Path("../input/annot.xlsx") # TODO: Fix to Take it from params
		annot.create_annot_file(fin,annot_f)
		# Output this in the input directory?

	paramsdf = pd.read_excel(params_i,index_col=0, engine='openpyxl')
	trees_dir, outdir = main(paramsdf)
	
	if params_clean:
		env_setup.clean_tmp()
		# TODO: Implement this

	# Visualizations
	params_viz = False
	if params_viz:
		phylogeny.visualize_trees() # TODO: Needs implementation
	if params_tree_render:
		visualizations.render_trees(trees_dir, outdir, mode=params_tree_render)
	else:
		visualizations.render_trees(trees_dir, outdir)
	end_dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
	print("Application started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))
	# logging.info("The scripts started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))
