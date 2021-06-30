import re
import os
import pathlib
import logging
from Bio import SeqIO

def file_exists(f):
	"""Empty docstring
	"""
	exist = os.path.exists(f)
	if not exist:
		logging.critical("%s does not exist" %(f))
		raise FileNotFoundError("%s does not exist" %(f))
	logging.info("%s exists" %(f))

def input_file_format(fasta_file_path):
	"""Check if the input query file is in supported format
	Currently supporting formats: Fasta
	"""
	check_handle = open(fasta_file_path, "r")

	if fasta_file_path.stat().st_size == 0:
		logging.critical(" Input file is empty, check input file ")
		raise Exception(" Input file is empty, check input file ")
	try:
		first_line = check_handle.readline().rstrip()
		if re.match(r'^>', first_line): # Is fasta
			pass
		else:
			logging.critical(" The input file is not fasta ")
			raise Exception(" The input file is not fasta ")
	except UnicodeError:
		logging.critical(" The input file is not in plain text file format ")
		raise Exception(" The input file is not in plain text file format ")
	
	logging.info( " Finished input file integrity check successfully ")

def filter_input_file(fasta_file_path,por_n=100,min_length=0):
	"""Need to perform a rename in the input file or create a file with renamed values
	Keep only the first \\S+ from each input sequence, to avoid conflicts
	Also filter for N like TRECS? TODO: Implement extra filters
	Input: Complete path of the fasta file (Posix | Windows Path)
	"""

	def _sequence_cleaner(fasta_file_path,por_n, min_length):
		"""min_length Default value 0, it means you don’t have to care about the minimum length
		 por_n the user defines the % of N is allowed. Default value 100, all sequences with ’N’ will be in your output, 
		 set value to 0 if you want no sequences with ”N” in your output
		"""
		# Create our hash table to add the sequences
		sequences = {}
	
		# Using the Biopython fasta parse we can read our fasta input
		for seq_record in SeqIO.parse(fasta_file_path, "fasta"):
			# Take the current sequence
			sequence = str(seq_record.seq).upper()
			# Check if the current sequence is according to the user parameters
			if (
				len(sequence) >= min_length
				and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n
			):
				# If the sequence passed in the test "is it clean?" and it isn't in the
				# hash table, the sequence and its id are going to be in the hash
				if sequence not in sequences:
					sequences[sequence] = seq_record.id
				# If it is already in the hash table, we're just gonna concatenate the ID
				# of the current sequence to another one that is already in the hash table
				else:
					sequences[sequence] += "_" + seq_record.id
	
		# Write the clean sequences
	
		# Create a file in the same directory where you ran this script
		fout = fasta_file_path.with_name("Filt_" + fasta_file_path.name)
		fname = fout.name
		with open(fout, "w+") as output_file:
			for sequence in sequences:
				output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
		return fname
	fname = _sequence_cleaner(fasta_file_path,por_n, min_length)
	logging.info(" Cleaned the fasta file %s at N percentage %d "%(fasta_file_path.name,por_n))
	return fname

def dir_priviledges(path):
	"""Check if the user has write priviledges in the specified directory
	Input: Posix | Window path from pathlib library
	"""
	write_priviledges = os.access(path, os.W_OK)
	if not write_priviledges:
		logging.critical(" Program has no write permission on %s " %(path))
		raise PermissionError(" You dont have the required priviledges in the specified output directory ")
	logging.info(" Permission to write on %s " %(path))

# TODO: Verify integrity of params file?
# Check if all the prerequired fields are present ? This means that it will make a large excel file

def clean_tmp_dir(path):
	pass
	logging.info(" Cleaning temporary directories ")