import re
import os
import pathlib
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline

"""
This functions works in three steps
1) the extraction of genes
2) creation of gene profiles
3) computation of phylogenetic trees
The exe keyword argument dictates if the function is executed or just returns the desired output
"""

def prepare_alns(outdir: pathlib.Path, query_f_path: pathlib.Path, 
	geneid_blast_res_xls: pathlib.Path, muscle_bin: str, exe: bool = True
) -> list:

	if not exe:
		batch_genes_dir = outdir / pathlib.Path(".tmp") /pathlib.Path("Gene_seqs")
		aln_files = batch_genes_dir.glob("**/*")
		aln_files = [f for f in aln_files if f.is_file()]
		return aln_files

	def _filter_manyN(sequence: str, N_perc: float) -> bool:
		"""
		Check if the provided gene sequence has many Ns and filter the parent organism from the analysis
		Default N_perc=50%
		input: Sequence as uppercase string
		return: True to keep, False to discard
		"""
		seq_Nperc = sequence.count("N")/len(sequence)
		if seq_Nperc >= N_perc:
			return False
		return True

	def _cut_genes(genomes_f:str, coords_df: pd.DataFrame, 
		N_perc:float =0.5
	) -> dict:
		"""
		Extract the gene sequences from the input fasta file
		Default N_perc=50% for _filter_manyN function
		input: batch file (complete genome), 
		pandas dataframe with gene coordinates
		return: dictionary of gene sequences
		"""
		final_seqs = {}
		seqs = SeqIO.to_dict(SeqIO.parse(genomes_f, "fasta"))
		orgs = list(seqs.keys())
		skipped_orgs = {}
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
				seq = str(seqs[org].seq[int(start)-1:int(end)])
				keep = _filter_manyN(seq, N_perc)
				if not keep:
					skipped_orgs[org] = ""
				final_seqs[gene][org] = seq
		final_seqs_filt = final_seqs.copy()
		for gene in final_seqs:
			for org in final_seqs[gene]:
				if org in skipped_orgs:
					del final_seqs_filt[gene][org]
		
		return final_seqs_filt
		
	def _align_sequences(gene_files : list, muscle_bin : pathlib.Path,
		method: str = "muscle"
	) -> list:
		"""
		Align gene sequences to create profiles
		input: gene_files, list of pathlib paths
		return: list of pathlib paths pointing to each alignment file
		"""
		aln_out_dir = outdir / pathlib.Path(".tmp") / pathlib.Path("Alignments")
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

	coords_df = pd.read_excel(geneid_blast_res_xls,engine="openpyxl") # Check for efficiency https://pandas.pydata.org/pandas-docs/stable/user_guide/scale.html
	batch_genes_dir = outdir / pathlib.Path(".tmp") / pathlib.Path("Gene_seqs")
	gene_files = []
	query_f_name = query_f_path.stem 
	gene_seqs = _cut_genes(query_f_path, coords_df)
	for gene in gene_seqs:
		fout = batch_genes_dir / pathlib.Path(query_f_name + "_" + gene + ".fa")
		gene_files.append(fout)
		fhandle = open(fout,"w")
		for k in gene_seqs[gene]:
			str_to_write = ">" + k + "\n" + str(gene_seqs[gene][k]) + "\n"
			fhandle.write(str_to_write)
		fhandle.close()
	aln_files = _align_sequences(gene_files, muscle_bin)
	return aln_files

def profile_aln(outdir : pathlib.Path, aln_files : list, 
	profiledb : str, muscle_bin : str, method : str ="muscle", 
	exe : bool =True
) -> list :
	"""
	Align each batch gene profile to the corresponding reference gene profile.
	Profile alignment input files MUST be in Fasta format
	input: outdir pathlib path pointing to output directoy
	aln_files, list of pathlib paths pointing to the batch gene alignment files
	return: list of pathlib paths pointing to the final alignments
	"""
	if not exe:
		profile_aln_dir = outdir / pathlib.Path("Profile_alns")
		profile_aln_files = profile_aln_dir.glob("**/*")
		profile_aln_files = [f for f in profile_aln_files if f.is_file()]
		return profile_aln_files
	
	profile_aln_dir = outdir / pathlib.Path("Profile_alns")
	profile_aln_files = []
	gene_regex = re.compile(r'^\S+_(\S+)_aln.fa$')
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

def build_trees(outdir: pathlib.Path, aln_files: list, phyml_bin: str, 
	seaview_bin: str, method: str = "PhyML", dist: str = "Kimura", 
	nj_bootstrap_repl: int = 1000, exe: bool = True
) -> pathlib.Path:
	"""
	If the organisms in the analysis are more than 20, the NJ trees will be computed
	"""
	if not exe:
		trees_dir = outdir / pathlib.Path("Phylogenetic_Trees")
		return trees_dir
	
	def _create_phylip_input(aln_f : pathlib.Path, 
		phylip_aln_f: pathlib.Path
	) -> None:
		"""
		Transform fasta files to phylip format
		PhyML needs phylip input format to run
		input: 
		aln_f: alignment file to convert
		phylip_aln_f: the phylip file
		"""
		alignment = AlignIO.read(aln_f, "fasta")
		for i in range(len(alignment)):
			full_name = alignment[i].id
			filt_name = full_name.replace(" ","_")
			alignment[i].id = filt_name
		AlignIO.write([alignment], phylip_aln_f, "phylip-relaxed")
		return 

	def _clean_phyml_output(phylip_aln_f: pathlib.Path, 
		trees_dir: pathlib.Path
	) -> None:
		"""
		Move *_phyml_tree.txt files to the appropriate directory and delete *_phyml_stats.txt
		"""
		phyml_tree_f = phylip_aln_f.parent / (phylip_aln_f.name + "_phyml_tree.txt")
		phyml_stats_f = phylip_aln_f.parent / (phylip_aln_f.name + "_phyml_stats.txt")
		phyml_tree_f_target = trees_dir / phyml_tree_f.name
		try:
			phyml_tree_f.replace(phyml_tree_f_target)
			phyml_stats_f.unlink()
		except:
			pass
		return 
	
	trees_dir = outdir / pathlib.Path("Phylogenetic_Trees")
	alns_dir = aln_files[0].parent
	if method == "PhyML":
		phylip_dir = alns_dir / pathlib.Path("phylip")
		phylip_dir.mkdir(exist_ok=True) # TODO: Think if exist_ok is fine
	for aln_f in aln_files:
		if method == "PhyML":
				phylip_aln_f = phylip_dir / aln_f.with_suffix(".phy").name
				_create_phylip_input(aln_f, phylip_aln_f)
				cmd = phyml_bin + " --quiet -o tl -s SPR -v estimated -m GTR -d nt -b -4 -f m -i " + str(phylip_aln_f)
				os.system(cmd)
				_clean_phyml_output(phylip_aln_f,trees_dir)

		if method == "BioNJ":
			tree_file_out = trees_dir / (aln_f.stem + "_NJ_tree.nwk")
			cmd = seaview_bin + " -build_tree -NJ -distance " \
			+ dist + " -replicates " + str(nj_bootstrap_repl) + " -o " + str(tree_file_out) + " " + str(aln_f)
			os.system(cmd)

	return trees_dir
