import re
import os
import pathlib
import logging
from Bio import SeqIO

def create_dirs(outdir, exist_ok = False):
		"""
		Perform the necessary actions to set up the script working environment
		"""
		script_out = outdir
		tmp_out = script_out / pathlib.Path(".tmp")
		tmp_genes_out = tmp_out / pathlib.Path("Gene_seqs") # Gene files of batch size
		tmp_alignments_out = tmp_out / pathlib.Path("Alignments") # Alignment of batch genes
		blast_results = script_out /  pathlib.Path("BlastResults")
		profiles_out = script_out / pathlib.Path("Profile_alns") # Final profile alns for each gene
		trees_out = script_out / pathlib.Path("Phylogenetic_Trees")
		graphics_out = script_out / pathlib.Path("Graphics")
		gene_sim_graphics = graphics_out / pathlib.Path("GeneIdentification")
		snp_graphics = graphics_out / pathlib.Path("LineageSpecificSNPs")
		
		script_out.mkdir(exist_ok = exist_ok, parents=True)
		tmp_out.mkdir(exist_ok = exist_ok)
		tmp_genes_out.mkdir(exist_ok = exist_ok)
		tmp_alignments_out.mkdir(exist_ok = exist_ok)
		blast_results.mkdir(exist_ok = exist_ok)
		profiles_out.mkdir(exist_ok = exist_ok)
		trees_out.mkdir(exist_ok = exist_ok)
		graphics_out.mkdir(exist_ok = exist_ok)
		gene_sim_graphics.mkdir(exist_ok = exist_ok)
		snp_graphics.mkdir(exist_ok = exist_ok)

def file_exists(f: pathlib.Path) -> str or None:
	exist = os.path.exists(f)
	if not exist:
		logging.critical(F"{f} does not exist")
		msg = F" {f} does not exist"
		return msg
	logging.info(F"{f} exists \u2705")

def input_file_format(fasta_file_path: pathlib.Path) -> str:
	# TODO: Add the error message to the GUI
	"""
	Check if the input query file is in supported format
	Currently supporting formats: Fasta
	"""
	check_handle = open(fasta_file_path, "r")

	if fasta_file_path.stat().st_size == 0:
		logging.critical("Input file is empty, check input file ")
		msg = " Input file is empty, check input file "
		return msg
	try:
		first_line = check_handle.readline().rstrip()
		if re.match(r'^>', first_line): # Is fasta
			pass
		else:
			logging.critical("The input file is not fasta ")
			msg = " The input file is not fasta "
			return msg
	except UnicodeError:
		logging.critical("The input file is not in plain text file format ")
		msg = " The input file is not in plain text file format "
		return msg
	logging.info(F"Finished input file integrity check successfully \u2705")

def dir_priviledges(path):
	"""
	Check if the user has write priviledges in the specified directory
	Input: Posix | Window path from pathlib library
	"""
	write_priviledges = os.access(path, os.W_OK)
	if not write_priviledges:
		logging.critical(F"Program has no write permission on {path} \u274C")
		raise PermissionError(" You dont have the required priviledges in the specified output directory ")
	logging.info(F"Have permission to write on {path} \u2705")


def filter_hpv16(fasta_file_path: pathlib.Path, seqs_to_include: list, outdir: pathlib.Path) -> str:
	"""
	Write the filtered output file in the output directroy with name: Filtered_" + fasta_file_path.name
	"""

	def _remove_problematic_seqs(fasta_file_path: pathlib.Path, seqs_to_include: list) -> dict:
		"""
		Clean the HPV16 sequences for non ATGC nucleotides
		All the non canonical nucleotides will be turned to N
		"""
		sequences_tmp = SeqIO.index(str(fasta_file_path), "fasta")
		sequences = {}
		for included_seq in seqs_to_include:
			sequences[included_seq] = str(sequences_tmp[included_seq].seq).upper()
		for seqname in sequences:
			sequence = sequences[seqname]
			for x in range(len(sequence)):
				char = sequence[x]
				if char != "A" and char != "T" and char != "G" and char != "C" and char != "N":
					sequence = sequence[:x] + "N" + sequence[x+1:]
			sequences[seqname] =  sequence
		return sequences
	
	sequences = _remove_problematic_seqs(fasta_file_path, seqs_to_include)
	# Write the clean sequences
	fout = outdir / ("Filtered_" + fasta_file_path.name)
	fname = fout.name
	with open(fout, "w+") as output_file:
		for seqname in sequences:
			output_file.write(">" + seqname + "\n" + sequences[seqname] + "\n")
	output_file.close()
	logging.info("Filtered non HPV-16 sequences")
	return fname