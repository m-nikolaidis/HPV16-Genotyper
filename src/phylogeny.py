import re
import os
import sys
import pathlib
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, TreeFace, TextFace




def prepare_batch_alns(outdir, query_f_path, geneid_blast_res_xls, muscle_bin, batch_size=1):
	"""Split the input fasta file in batches of size N.
	Each batch will be utilized for 
	1) the extraction of genes
	2) creation of gene profiles
	3) computation of phylogenetic trees
	input: #TODO: Add this
	return: None
	"""
	def _split_batches(query_f_path, outdir, batch_size):
		"""Split the batches using a generator function to save memory.
		input: #TODO: Add this
		output: Batch files in tmp directory, File with sequence allocations for easy tracking
		"""
		def _batch_iterator(iterator, batch_size):
			"""Returns lists of length batch_size.
		
			This can be used on any iterator, for example to batch up
			SeqRecord objects from Bio.SeqIO.parse(), or simply
			lines from a file handle.
		
			This is a generator function, and it returns lists of the
			entries from the supplied iterator.  Each list will have
			batch_size entries, although the final list may be shorter.
			"""
			entry = True  # Make sure we loop once
			while entry:
				batch = []
				while len(batch) < batch_size:
					try:
						entry = next(iterator)
					except StopIteration:
						entry = None
					if entry is None: # End of file
						break
					batch.append(entry)
				if batch:
					yield batch
		
		record_iter = SeqIO.parse(open(query_f_path), "fasta")
		for i, batch in enumerate(_batch_iterator(record_iter, batch_size)):
			tmpfout = outdir / pathlib.Path("Batch_%i.fa" % (i + 1)) #TODO: Use pathlib
			#TODO, write a file with which sequences are contained in each batch file
			with open(tmpfout, "w") as handle:
				count = SeqIO.write(batch, handle, "fasta")


	def _cut_genes(genomes_f, coords_df):
		"""Extract the gene sequences for each batch fasta file produced by _split_batches()
		input: batch file (complete genome), pandas dataframe with gene coordinates
		return: dictionary of gene sequences
		"""
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
		
	def _align_sequences(gene_files,muscle_bin,method="muscle"):
		"""Align batch gene sequences to create profiles
		input: gene_files, list of pathlib paths
		return: list of pathlib paths pointing to each alignment file
		"""
		aln_out_dir = outdir / pathlib.Path("tmp_dir/Alignments")
		aln_files = []
		for f in gene_files:
			fout = aln_out_dir / (f.stem + "_aln.fa")
			if method == "muscle":
				aln_cline = MuscleCommandline(cmd = muscle_bin, input=f, out=fout)
			else:
				aln_cline = None
				raise Exception(" No other alignment method is implemented yet ")
			aln_cline()
			aln_files.append(fout)
		return aln_files

	batch_fa_dir = outdir / pathlib.Path("tmp_dir/Batch_fasta_files")
	_split_batches(query_f_path,batch_fa_dir,batch_size)

	coords_df = pd.read_excel(geneid_blast_res_xls,engine="openpyxl") # Check for efficiency https://pandas.pydata.org/pandas-docs/stable/user_guide/scale.html
	batch_genes_dir = outdir / pathlib.Path("tmp_dir/Batch_genes")
	batch_fa_files = os.listdir(batch_fa_dir)
	batch_fa_files = [batch_fa_dir / pathlib.Path(f) for f in batch_fa_files]
	gene_files = []
	for batch_f in batch_fa_files:
		batchname = batch_f.stem 
		gene_seqs = _cut_genes(batch_f, coords_df)
		for gene in gene_seqs:
			fout = batch_genes_dir / pathlib.Path(batchname + "_" + gene + ".fa")
			gene_files.append(fout)
			fhandle = open(fout,"w")
			for k in gene_seqs[gene]:
				str_to_write = ">" + k + "\n" + str(gene_seqs[gene][k]) + "\n"
				fhandle.write(str_to_write)
			fhandle.close()
	
	aln_files = _align_sequences(gene_files, muscle_bin)
	return aln_files

def profile_aln(outdir, aln_files, profiledb, muscle_bin, method="muscle"):
	"""Align each batch gene profile to the corresponding reference gene profile.
	Profile alignment input files MUST be in Fasta format
	input: outdir pathlib path pointing to output directoy
	aln_files, list of pathlib paths pointing to the batch gene alignment files
	return: list of pathlib paths pointing to the final alignments
	"""

	aln_out_dir = outdir / pathlib.Path("tmp_dir/Alignments")
	profile_aln_dir = outdir / pathlib.Path("Profile_alns")
	profile_aln_files = []
	gene_regex = re.compile(r'^Batch_\d+_(\S+)_aln.fa$')
	for f in aln_files:
		fout = profile_aln_dir / f.name
		profile_aln_files.append(fout)
		m = re.match(gene_regex,f.name)
		gene = m.group(1)
		db_gene_f = pathlib.Path(profiledb) / pathlib.Path(gene + "_profile.fa")
		if method == "muscle":
			aln_cline = MuscleCommandline(cmd = muscle_bin, in1=f, in2=db_gene_f, out=fout, profile=True)
		else:
			aln_cline = None
			raise Exception(" No other alignment method is implemented yet ")
		aln_cline()
	return profile_aln_files

def build_trees(outdir, aln_files,phyml_bin, seaview_bin, method="PhyML", dist="Kimura", nj_bootstrap_repl=1000):
	"""
	
	"""

	def _create_phylip_input(aln_f, phylip_aln_f):
		"""Transform fasta files to phylip format
		PhyML needs phylip input format to run
		input: 
		aln_f: pathlib path of the alignment file to convert
		phylip_aln_f: pathlib path of the phylip file
		"""
		# The relaxed phylip format needs for the sequence headers to not have spaces
		# CARE! Might create incompatabilities later?
		
		alignment = AlignIO.read(aln_f, "fasta")
		for i in range(len(alignment)):
			full_name = alignment[i].id
			filt_name = full_name.replace(" ","_")
			alignment[i].id = filt_name
		AlignIO.write([alignment], phylip_aln_f, "phylip-relaxed")


	def _clean_phyml_output(phylip_aln_f, trees_dir):
		"""Move _phyml_tree.txt files to the appropriate directory and remove _phyml_stats.txt
		"""
		phyml_tree_f = phylip_aln_f.parent / (phylip_aln_f.name + "_phyml_tree.txt")
		phyml_stats_f = phylip_aln_f.parent / (phylip_aln_f.name + "_phyml_stats.txt")
		phyml_tree_f_target = trees_dir / phyml_tree_f.name
		try:
			phyml_tree_f.replace(phyml_tree_f_target)
			phyml_stats_f.unlink()
		except:
			pass
	trees_dir = outdir / pathlib.Path("Phylogenetic_trees")
	alns_dir = aln_files[0].parent
	if method == "PhyML":
		phylip_dir = alns_dir / pathlib.Path("phylip")
		# TODO: Uncomment after testing
		# phylip_dir.mkdir() # Write error in logfile if folder exists
	for aln_f in aln_files:
		if method == "PhyML":
				phylip_aln_f = phylip_dir / aln_f.with_suffix(".phy").name
				_create_phylip_input(aln_f, phylip_aln_f)
				cmd = phyml_bin + " --quiet -o tl -s SPR -v estimated -m GTR -d nt -b -4 -f m -i " + str(phylip_aln_f)
				# os.system(cmd)
				_clean_phyml_output(phylip_aln_f,trees_dir)
		if method == "BioNJ":
			tree_file_out = trees_dir / (aln_f.stem + "_NJ_tree.nwk")
			cmd = seaview_bin + " -build_tree -NJ -distance " \
			+ dist + " -replicates " + str(nj_bootstrap_repl) + " -o " + str(tree_file_out) + " " + str(aln_f)
			# os.system(cmd)

	return trees_dir

def visualize_trees(trees_dir, batch):
	"""
	TODO: Finish this
	"""
	pass
def render_trees(trees_dir, batch):
	"""
	TODO: Finish this
	"""
	pass