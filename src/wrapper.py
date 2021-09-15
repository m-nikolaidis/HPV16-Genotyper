import os
import sys
import pathlib
import logging
import pandas as pd
import multiprocessing
from Bio import SeqIO
# Tool modules
import annot
import blast
import env_setup
import phylogeny

def _init_binaries(system: sys.platform) -> list:
	"""
	Initialize the paths for the needed binaries depending on the system
	input: System info
	return: List of string
	"""
	if "win" in system:
		makeblastdb_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "makeblastdb.exe"
		blastn_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "blastn.exe"
		muscle_bin= pathlib.Path(__file__).parent / "resources" / "Win_bin" / "muscle3.8.31.exe"
		seaview_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "seaview.exe"
		phyml_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "PhyML-3.1_win32.exe"
	if "linux" in system:
		makeblastdb_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "makeblastdb"
		blastn_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "blastn"
		muscle_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "muscle3.8.31_i86linux"
		seaview_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "seaview4"
		phyml_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "PhyML-3.1_linux64"
	return [str(makeblastdb_bin), str(blastn_bin), str(muscle_bin), 
			str(seaview_bin), str(phyml_bin)
			]

def _defaultparams() -> dict:
	system = sys.platform
	input_dir = ""
	query_f = ""
	outdir = pathlib.Path(__file__).absolute().parent / pathlib.Path("Script_out")
	threads_to_use =  multiprocessing.cpu_count() - 2
	params = {
		"system": system,
		"in": input_dir,
		"out": outdir,
		"num_threads": threads_to_use,
		"query": query_f,
		"SNP_identification_Evalue": 0.005,
		"SNP_identification_Word_Size": 4,
		"SNP_annotation_file": pathlib.Path(__file__).parent / 'resources'/ 'SNPs_annotation.xlsx',
		"SNP_db_path": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'NC_001526_probes.fa',
		"cSNP_db_path": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'NC_001526_cancer_probes.fa',
		"Gene_identification_Evalue": 1e-5, 
		"Gene_identification_Word_size": 7,
		"GenesProfile_directory": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'profiles',
		"GenesRef_database":pathlib.Path(__file__).parent / 'resources' / 'sequences' / '16refs_gene_db.fa',
		"SimplotRef_database":pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'Reference_genomes_profile_mafft_GINSI.fa',
		"HPV_filter_db": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'alphapapillomavirus9.fa'
	}
	return params

def main(paramsdf: pd.DataFrame, exe: bool = True) -> pathlib.Path:
	# Basic variables
	hpv16error = False # To stop execution if no hpv16 sequences are found
	outdir = pathlib.Path(paramsdf.loc["out","Value"])
	env_setup.create_dirs(outdir, exist_ok=True)
	system = paramsdf.loc["system","Value"]
	indir = pathlib.Path(paramsdf.loc["in","Value"])
	query_f = pathlib.Path(paramsdf.loc["query","Value"])
	query_f_path = indir / query_f
	makeblastdb_bin,blastn_bin,muscle_bin,seaview_bin,phyml_bin = _init_binaries(system)
	annot_f = pathlib.Path(paramsdf.loc["SNP_annotation_file","Value"])
	threads = paramsdf.loc["num_threads", "Value"]

	# Grab the root logger instance and use it for logging
	logfile = outdir / pathlib.Path(".logfile.log")
	logger = logging.getLogger()
	fhandler = logging.FileHandler(filename=logfile, encoding="UTF-8")
	formatter = logging.Formatter("%(asctime)s\t%(message)s")
	fhandler.setFormatter(formatter)
	logger.addHandler(fhandler)
	logger.setLevel(logging.DEBUG)
	logging.info("#### DO NOT DELETE THIS FILE! ###")
	logging.info(F"Environment set successfully in {outdir} \u2705")
	logging.info("Parameters used in this run: ")
	
	indeces = paramsdf.index
	for index in indeces:
		val = paramsdf.loc[index,"Value"]
		logging.info(F"param: {index},{val}")

	# Environment setup
	env_setup.file_exists(query_f_path)
	env_setup.input_file_format(query_f_path)
	env_setup.dir_priviledges(outdir)


	# Filter input file
	workflow = "HPV16 filter"
	hpv16filt_results = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, exe=exe)
	hpv16_seqs = blast.identify_nonHPV16(hpv16filt_results)
	if len(hpv16_seqs) == 0:
		hpv16error = True
		return query_f_path, hpv16error
	query_f_filt = env_setup.filter_hpv16(query_f_path, hpv16_seqs, outdir)
	query_f = query_f_filt # Posix | Windows Path
	query_f_path = outdir / query_f
	logging.info(F"param: query,{str(query_f_path)}")
	query_f_index = SeqIO.index(str(query_f_path), "fasta")
	
	# BLASTn search
	workflow = "gene identification"
	geneid_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, exe=exe)
	geneid_blastn_res_path_xlsx, blastDF = blast.parse_GeneID_results(geneid_blastn_res_path,paramsdf, exe=exe)

	# Filter organisms that have many Ns inside genes
	aln_files = phylogeny.prepare_alns(outdir, query_f_path, geneid_blastn_res_path_xlsx, muscle_bin, exe_threads = threads, exe=exe)
	# records = [query_f_index[org] for org in query_f_index if org in kept_orgs]
	# with open(query_f_path, "w") as output_handle:
	# 	SeqIO.write(records, output_handle, "fasta")
	# output_handle.close()
	# query_f_index = SeqIO.index(str(query_f_path), "fasta")
	
	# Continue with renewed fasta_file
	workflow = "snp"
	snp_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, exe=exe)
	snp_blastn_res_path_xlsx, snpDF = blast.parse_SNP_results(snp_blastn_res_path, query_f_index, exe=exe)
	workflow = "cancer"
	C_snp_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, exe=exe)
	blast.parse_SNP_results(C_snp_blastn_res_path, query_f_index, exe=exe, cancer=True)
	annot.annotate_results(annot_f,snp_blastn_res_path_xlsx,exe=exe)
	
	blast.find_recombinants(blastDF, snpDF, outdir)

	# Gene alignments and trees
	profiledb_dir = paramsdf.loc["GenesProfile_directory","Value"]
	aln_files = phylogeny.profile_aln(outdir, aln_files, profiledb_dir, muscle_bin, exe = exe)
	
	organisms = len(query_f_index.keys()) # final seqs to analyze
	if organisms > 10: 
		method = "BioNJ" 
	else: 
		method = "PhyML"
	phylogeny.build_trees(outdir, aln_files, phyml_bin, seaview_bin, threads = threads, exe=exe, method = method)

	return query_f_path