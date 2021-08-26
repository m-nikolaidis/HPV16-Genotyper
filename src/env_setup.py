import re
import os
import pathlib
import logging
from Bio import SeqIO

def create_dirs(outdir, exist_ok=False):
		"""
		Perform the necessary actions to set up the script working environment
		"""
		script_out = outdir
		tmp_out = script_out / pathlib.Path("tmp_dir")
		tmp_genes_out = tmp_out / pathlib.Path("Gene_seqs") # Gene files of batch size
		tmp_alignments_out = tmp_out / pathlib.Path("Alignments") # Alignment of batch genes
		profiles_out = script_out / pathlib.Path("Profile_alns") # Final profile alns for each gene
		trees_out = script_out / pathlib.Path("Phylogenetic_Trees")
		graphics_out = script_out / pathlib.Path("Graphics")
		gene_sim_graphics = graphics_out / pathlib.Path("GeneSim")
		snp_graphics = graphics_out / pathlib.Path("SNP_detection")
		tree_render = graphics_out / pathlib.Path("Trees_images")
		
		script_out.mkdir(exist_ok=True)
		tmp_out.mkdir(exist_ok=exist_ok)
		tmp_genes_out.mkdir(exist_ok=exist_ok)
		tmp_alignments_out.mkdir(exist_ok=exist_ok)
		profiles_out.mkdir(exist_ok=exist_ok)
		trees_out.mkdir(exist_ok=exist_ok)
		graphics_out.mkdir(exist_ok=exist_ok)
		gene_sim_graphics.mkdir(exist_ok=exist_ok)
		snp_graphics.mkdir(exist_ok=exist_ok)
		tree_render.mkdir(exist_ok=exist_ok)

def file_exists(f):
	"""Empty docstring
	"""
	exist = os.path.exists(f)
	if not exist:
		logging.critical(F"{f} does not exist")
		raise FileNotFoundError(F" {f} does not exist")
	logging.info(F"{f} exists \u2705")

def input_file_format(fasta_file_path):
	"""Check if the input query file is in supported format
	Currently supporting formats: Fasta
	"""
	check_handle = open(fasta_file_path, "r")

	if fasta_file_path.stat().st_size == 0:
		logging.critical("Input file is empty, check input file ")
		raise Exception(" Input file is empty, check input file ")
	try:
		first_line = check_handle.readline().rstrip()
		if re.match(r'^>', first_line): # Is fasta
			pass
		else:
			logging.critical("The input file is not fasta ")
			raise Exception(" The input file is not fasta ")
	except UnicodeError:
		logging.critical("The input file is not in plain text file format ")
		raise Exception(" The input file is not in plain text file format ")
	
	logging.info(F"Finished input file integrity check successfully \u2705")

def filter_input_file(fasta_file_path,por_n=100,min_length=0):
	"""Need to perform a rename in the input file or create a file with renamed values
	Keep only the first \\S+ from each input sequence, to avoid conflicts
	Also filter for N like TRECS? TODO: Implement extra filters
	Input: Complete path of the fasta file (Posix | Windows Path)
	"""

	def _sequence_cleaner(fasta_file_path):
		"""Clean the sequences for non ATGC nucleotides
		All the non canonical nucleotides will be turned to N
		"""
		sequences = SeqIO.to_dict(SeqIO.parse(fasta_file_path, "fasta"))
		for seq_record in sequences:
			sequence = str(sequences[seq_record].seq).upper()
			for x in range(len(sequence)):
				char = sequence[x]
				if char != "A" and char != "T" and char != "G" and char != "C" and char != "N":
					sequence = sequence[:x] + "N" + sequence[x+1:]
			sequences[seq_record].seq =  sequence
		return sequences
	def _remove_problematic_seqs(sequences,por_n, min_length):
		"""min_length Default value 0, it means you don’t have to care about the minimum length
		 por_n the user defines the % of N is allowed. Default value 100, all sequences with ’N’ will be in your output, 
		 set value to 0 if you want no sequences with ”N” in your output
		"""
		# Create our hash table to add the sequences
		final_sequences = {}
	
		# Using the Biopython fasta parse we can read our fasta input
		for seqname in sequences:
			# Take the current sequence
			sequence = str(sequences[seqname].seq).upper()
			# Check if the current sequence is according to the user parameters
			if (
				len(sequence) >= min_length
				and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n
			):
				# If the sequence passed in the test "is it clean?" and it isn't in the
				# hash table, the sequence and its id are going to be in the hash
				if sequence not in final_sequences:
					final_sequences[sequence] = seqname
				# If it is already in the hash table, we're just gonna concatenate the ID
				# of the current sequence to another one that is already in the hash table
				else:
					final_sequences[sequence] += "_" + seqname
	
		# Write the clean sequences
	
		# Create a file in the same directory where you ran this script
		fout = fasta_file_path.with_name("Filt_" + fasta_file_path.name)
		fname = fout.name
		with open(fout, "w+") as output_file:
			for sequence in final_sequences:
				output_file.write(">" + final_sequences[sequence] + "\n" + sequence + "\n")
		return fname
	sequences = _sequence_cleaner(fasta_file_path)
	fname = _remove_problematic_seqs(sequences,por_n, min_length)
	logging.info("Cleaned the fasta file %s at N percentage %d "%(fasta_file_path.name,por_n))
	return fname

def dir_priviledges(path):
	"""Check if the user has write priviledges in the specified directory
	Input: Posix | Window path from pathlib library
	"""
	write_priviledges = os.access(path, os.W_OK)
	if not write_priviledges:
		logging.critical(F"Program has no write permission on {path} \u274C")
		raise PermissionError(" You dont have the required priviledges in the specified output directory ")
	logging.info(F"Have permission to write on {path} \u2705")

def clean_tmp_dir(path):
	pass
	logging.info("Cleaning - temporary directories ")