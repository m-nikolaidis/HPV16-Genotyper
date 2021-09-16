import re
import os
import sys
import shlex
import pathlib
import logging
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
from multiprocessing.pool import ThreadPool

def _callMultiThreadProcc(cmd):
			system = sys.platform
			if "win" in system:
				p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			else:
				p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			p.communicate()

def prepare_alns(outdir: pathlib.Path, query_f_path: pathlib.Path, 
	geneid_blast_res_xls: pathlib.Path, muscle_bin: str, profiledb_dir: pathlib.Path, exe_threads: int,
	exe: bool = True,
) -> list:
	if not exe:
		batch_genes_dir = outdir / pathlib.Path(".tmp") /pathlib.Path("Gene_seqs")
		aln_files = batch_genes_dir.glob("**/*")
		aln_files = [f for f in aln_files if f.is_file()]
		return aln_files
	
	def _cut_genes(genomes_f: str, coords_df: pd.DataFrame, 
	) -> dict:
		"""
		Extract the gene sequences from the input fasta file
		"""
		final_seqs = {}
		seqs = SeqIO.index(str(genomes_f), "fasta")
		orgs = list(seqs.keys())
		for org in orgs:
			tmpdf = coords_df[coords_df["Query sequence"] == org]
			indeces = tmpdf.index
			for idx in indeces:
				gene = tmpdf.loc[idx,"Gene"]
				if gene not in final_seqs:
					final_seqs[gene] = {}
				if org not in final_seqs[gene]:
					final_seqs[gene][org] = ""
				start = tmpdf.loc[idx,"Query start"]
				end = tmpdf.loc[idx,"Query end"]
				seq = str(seqs[org].seq[int(start)-1:int(end)]).upper()
				final_seqs[gene][org] = seq
		return final_seqs		

	def _align_sequences(gene_files : list, muscle_bin : str,
		threads: int, outdir : pathlib.Path , profiledb : str, method: str = "muscle"
	) -> list :
		"""
		Align gene sequences to create profiles
		input: gene_files, list of pathlib paths
		return: list of pathlib paths pointing to each alignment file
		"""
		aln_out_dir = outdir / pathlib.Path("Alignments")
		aln_files = []
		pool = ThreadPool(threads)
		if method != "muscle":
			raise Exception(" No other alignment method is implemented yet ")
		gene_regex = re.compile(r'^\S+_([E|L]\d+).fa$')
		for f in gene_files:
			m = re.match(gene_regex,f.name)
			gene = m.group(1)
			fout = aln_out_dir / (f.stem + "_aln.fa")
			db_gene_f = profiledb / pathlib.Path(gene + "_profile.fa")
			arguments = " -profile -in1 " + str(f) + \
				" -in2 " + str(db_gene_f) + " -out " + str(fout)
			pool.apply_async(_callMultiThreadProcc, (muscle_bin + arguments,))
			aln_files.append(fout)
		pool.close()
		pool.join()		
		return aln_files

	logging.info(F"Creating the alignment files for each gene")
	coords_df = pd.read_excel(geneid_blast_res_xls,engine="openpyxl") 
	# Check for efficiency https://pandas.pydata.org/pandas-docs/stable/user_guide/scale.html
	batch_genes_dir = outdir / pathlib.Path(".tmp") / pathlib.Path("Gene_seqs")
	gene_files = []
	gene_seqs = _cut_genes(query_f_path, coords_df)
	for gene in gene_seqs:
		for org in gene_seqs[gene]:
			fout = batch_genes_dir / pathlib.Path(org + "_" + gene + ".fa")
			gene_files.append(fout)
			fhandle = open(fout,"w")
			# Maybe create SeqRecord objects to write? Faster IO? Less memory need?
			# TODO: Check
			# Try in python notebook to keep a record of all of these stuff
			str_to_write = ">" + org + "\n" + str(gene_seqs[gene][org]) + "\n"
			fhandle.write(str_to_write)
			fhandle.close()
	aln_files = _align_sequences(gene_files, muscle_bin, profiledb = profiledb_dir, threads = exe_threads, outdir = outdir)
	logging.info(F"Finished aligning files")
	return aln_files

def build_trees(outdir: pathlib.Path, aln_files: list, seaview_bin: str, 
	threads: int, method: str = "BioNJ", dist: str = "Kimura", 
	nj_bootstrap_repl: int = 1000, exe: bool = True
) -> pathlib.Path:
	"""
	Compute the Neighbour joining phylogenetic trees with 
	predefined parameters
	"""
	if not exe:
		trees_dir = outdir / pathlib.Path("Phylogenetic_Trees")
		return trees_dir
	
	def _clean_nj_output(tree_f: pathlib.Path, nj_regex: re.compile) -> None:
		"""
		Removes the leading info generated from SeaView in CLI mode
		----> [NJ \d+ sites Kimura \d+ repl.] <----
		"""
		text = tree_f.read_text()
		m = re.match(nj_regex,text)
		new_text = m.group(1)
		tree_f.write_text(new_text)

	trees_dir = outdir / pathlib.Path("Phylogenetic_Trees")
	
	logging.info(F"Initiating {method} tree calculation")
	if method == "BioNJ":
		pool = ThreadPool(threads)
		nj_regex = re.compile(r"^\[NJ \d+ sites Kimura.+\] (\S+)")
		results = []
		for aln_f in aln_files:
			tree_file_out = trees_dir / (aln_f.stem + "_NJ_tree.nwk")
			arguments = " -build_tree -NJ -distance " + dist \
				+ " -replicates " + str(nj_bootstrap_repl) + " -o " \
				+ str(tree_file_out) + " " + str(aln_f)
			results.append(pool.apply_async(_callMultiThreadProcc, (seaview_bin + arguments,)))
		pool.close()
		pool.join()
		for aln_f in aln_files:
			tree_file_out = trees_dir / (aln_f.stem + "_NJ_tree.nwk")
			_clean_nj_output(tree_file_out, nj_regex)
	logging.info(F"Finished")
	return trees_dir