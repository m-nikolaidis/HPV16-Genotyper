import os
from Bio import SeqIO

def _check_file_exists(f,outdir):
	"""
	Filelist is imported from blast.py return
	"""
	exist = os.path.exists(f)
	if exist == True:
		return True
	logging.debug("Blastn result file does not exist in %s" %(outdir))
	raise Exception("Blastn result file does not exist in %s" %(outdir))

def _check_input_fa_format():
	pass

def sequence_cleaner(fasta_file, outdir, min_length=0, por_n=100):
	# min_length Default value 0, it means you don’t have to care about the minimum length
	# por_n the user defines the % of N is allowed. Default value 100, all sequences with ’N’ will be in your ouput, set value to 0 if you want no sequences with ”N” in your output
		# Create our hash table to add the sequences
	sequences = {}

	# Using the Biopython fasta parse we can read our fasta input
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
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
	with open("clear_" + fasta_file, "w+") as output_file:
		# Just read the hash table and write on the file as a fasta format
		for sequence in sequences:
			output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

	print("CLEAN!!!\nPlease check clear_" + fasta_file)


def split_file(fin, indir, outdir, batch_size=1):
	# Split the fasta file to batches of size N for the phylogeny
	os.mkdir(outdir + "/tmp_dir/" + "Batch_fasta_files") # TODO: use pathlib for windows compatibility
	def batch_iterator(iterator, batch_size):
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
				if entry is None:
					# End of file
					break
				batch.append(entry)
			if batch:
				yield batch
	
	input_fa = indir + fin
	record_iter = SeqIO.parse(open(input_fa), "fasta")
	for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
		filename = outdir + "/tmp_dir/Batch_fasta_files/" + "Batch_%i.fasta" % (i + 1) #TODO: Use pathlib
		#TODO, write a file with which sequences are contained in each batch file
		with open(filename, "w") as handle:
			count = SeqIO.write(batch, handle, "fasta")
		print("Wrote %i records to %s" % (count, filename))