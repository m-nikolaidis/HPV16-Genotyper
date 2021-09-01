import re
import logging
import pathlib
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

def blastn_search(exe_params,query_f_path, workflow,makeblastdb_bin,blastn_bin,exe=True):
	"""
	Execute blast search
	Input query is ref vector fasta file, input database is genomes fasta for genotyping
	The params dict will be used for the parameters of blast
	Exe is a boolen to execute the command or not. Added it so i can just return the path in case the script fails after the search
	"""
	threads = exe_params.loc['num_threads',"Value"]
	outdir = pathlib.Path(exe_params.loc['out',"Value"])


	if workflow == "HPV16 filter":
		logging.info("Starting - BLASTn search to identify non HPV16 sequences")
		db_path = pathlib.Path(exe_params.loc['HPV_filter_db',"Value"])
		db_file = db_path.name
		blastres_f = query_f_path.stem + "_HPV16_Blastn_res.xml"
		evalue = exe_params.loc['Gene_identification_Evalue',"Value"]
		word_size = exe_params.loc['Gene_identification_Word_size',"Value"]
	
	if workflow == "gene identification":
		logging.info("Starting - BLASTn search for gene identification")
		db_path = pathlib.Path(exe_params.loc['GenesRef_database',"Value"])
		db_file = db_path.name
		blastres_f = "GeneIdentification_Blastn_results.xml"
		evalue = exe_params.loc['Gene_identification_Evalue',"Value"]
		word_size = exe_params.loc['Gene_identification_Word_size',"Value"]
	
	if workflow == "snp":
		logging.info("Starting - BLASTn search for Lineage specific SNPs")
		db_path = pathlib.Path(exe_params.loc['SNP_db_path',"Value"])
		db_file = db_path.name
		blastres_f = "lineageSpecificProbes_Blastn_results.xml"
		evalue = exe_params.loc["SNP_identification_Evalue","Value"]
		word_size = exe_params.loc["SNP_identification_Word_Size","Value"]

	if workflow == "cancer":
		logging.info("Starting - BLASTn search for increased cancer risk SNPs")
		db_path = pathlib.Path(exe_params.loc['cSNP_db_path',"Value"])
		db_file = db_path.name
		blastres_f = "cancerProbes_Blastn_results.xml"
		evalue = exe_params.loc["SNP_identification_Evalue","Value"]
		word_size = exe_params.loc["SNP_identification_Word_Size","Value"]
	
	tmp = pathlib.Path("tmp_dir")
	dbout = outdir / tmp / db_file
	blastres_f_path = outdir / blastres_f

	makedb = NcbimakeblastdbCommandline(cmd = makeblastdb_bin,
										dbtype ="nucl",
										input_file = db_path,
										out = dbout
	)
	blastn_search = NcbiblastnCommandline(cmd = blastn_bin,
											query = query_f_path,
											db = dbout,
											outfmt = 5,
											max_target_seqs = 100,
											out = blastres_f_path,
											evalue = evalue,
											num_threads = threads,
											word_size = word_size,
											dust="no"
	)
	if exe:
		makedb()
		blastn_search()
		logging.info("Finished - BLASTn search ")
	return blastres_f_path

def _sort_alphanumeric(iteratable):
	# Helper func Sort the given motif list alphanumerically :return: sorted list 
	int_convert = lambda text: int(text) if text.isdigit() else text 
	sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(iteratable, key = sorting_key)


def parse_GeneID_results(blastres_f_path: pathlib.Path, 
		params: pd.DataFrame, 
		exe:bool =True
	) -> pathlib.Path:
	"""
	Writes excel file
	"""
	if not exe:
		geneID_results_path = blastres_f_path.with_name(blastres_f_path.stem + "_GeneID.xlsx")
		return geneID_results_path
	
	def _get_db_lengths(db_path):
		dblengths = {}
		seqrecords = SeqIO.parse(db_path,"fasta")
		for seqrecord in seqrecords:
			name = seqrecord.id
			seqlen = len(str(seqrecord.seq))
			dblengths[name] = seqlen
		return dblengths
	logging.info("Parsing - gene identification results ")
	
	db_path = pathlib.Path(params.loc["GenesRef_database","Value"])
	dblengths = _get_db_lengths(db_path)
	blast_records= NCBIXML.parse(open(blastres_f_path))
	blast_res = {}
	for blast_record in blast_records: # Each blast record is the total res for a specific query
		qorg = blast_record.query
		if qorg not in blast_res:
			blast_res[qorg] = {}

		for alignment in blast_record.alignments: # One aln for each query - subject hit
			subjct = alignment.hit_def # alignment.hit_def is the subject name
			subjct_name = subjct[:-3]
			gene = subjct[-2:]
			if gene not in blast_res[qorg]:
				blast_res[qorg][gene] = {
				"qaccver":qorg,
				"saccver":None,
				"qstart":0,
				"qend":0,
				"sstart":0,
				"send":0,
				"length":0,
				"pident":0,
				"evalue":1,
				"subject_covhsp":0,
				# "subject_covtotal":0,
				}
			# alignment.length is the total aln length?
			for hsp in alignment.hsps: # Each query - subject hit has multiple hsps (blast hits)
				qstart =  hsp.query_start
				qend =  hsp.query_end
				sstart =  hsp.sbjct_start
				send =  hsp.sbjct_end
				evalue = hsp.expect 
				hsp_len = hsp.align_length
				pident =  hsp.identities / hsp.align_length # As calculated in NCBI blast2seq
				if evalue < blast_res[qorg][gene]["evalue"]: # Keep only the best hsp result
					blast_res[qorg][gene]["saccver"] = subjct
					blast_res[qorg][gene]["qstart"] = qstart
					blast_res[qorg][gene]["qend"] = qend
					blast_res[qorg][gene]["sstart"] = sstart 
					blast_res[qorg][gene]["send"] = send
					blast_res[qorg][gene]["evalue"] = evalue
					blast_res[qorg][gene]["length"] = hsp_len
					blast_res[qorg][gene]["pident"] = round((pident*100),3)
					blast_res[qorg][gene]["subject_covhsp"] = round((hsp_len/dblengths[subjct])*100,3)

	df_list = []
	for k in blast_res:
		tmpdf = pd.DataFrame.from_dict(blast_res[k],orient='index')
		df_list.append(tmpdf)

	blastdf = pd.concat(df_list)
	blastdf.index = range(len(blastdf))
	blastdf["Sequences"] = blastdf["qaccver"]
	indeces = blastdf.index
	for idx in indeces:
		saccver = blastdf.loc[idx,"saccver"]
		db, gene = saccver.split("_")
		blastdf.loc[idx,"Database"] = db
		blastdf.loc[idx,"Gene"] = gene
		blastdf.loc[idx,"Lineage"] = db[0] 
	geneID_results_path = blastres_f_path.parent / "GeneIdentification_results.xlsx"
	blastdf.to_excel(geneID_results_path,index=False)
	return geneID_results_path


def parse_SNP_results(blastres_f_path: pathlib.Path, 
	exe:bool = True, cancer:bool = False
	) -> pathlib.Path:
	"""
	Reads the blast output files and ...
	TODO: Write me
	cancer parameter is used to differentiate from Lineage specific SNPs and cancer SNPs
	Return None
	"""
	
	if cancer:
		probe_results_path = blastres_f_path.parent / "cancerSNP_results.xlsx"
		logging.info("Parsing - increased cancer risk SNP results ")
	else:
		probe_results_path = blastres_f_path.parent / "LineageSpecificSNPs.xlsx"
		logging.info("Parsing - lineage specific SNP results ")
	if not exe:
		return probe_results_path
	
	def _get_aln_nucl(query_aln,subjct_aln,sstart,probe_len):
		# I have used this approach for the case that gaps are inserted in the subject seq
		# If not return fixed position
		# TODO CHECK IF IT IS POSSIBLE
		sbjct_pos = 0
		middle_pos = (round(probe_len/2) - sstart) + 1 # Round so in case of odd number, i get int
		# This is not right. Because i dont know which part of the
		# probe has the hit
		# Need to fix it.
		# Should check the probe in the total length (subject probe fasta) to figure out
		# which part is aligned. And then get the original probe middle nucleotide
		# TODO: Fix me
		for x in range(len(subjct_aln)):
			if subjct_aln[x] != "-":
				sbjct_pos += 1
				if sbjct_pos == middle_pos:
					return (query_aln[x], sbjct_pos)
		return ["X",1]

	probe_len = 31 # TODO; should be calculated automatically or provided
	probe_table_res = {}
	blast_records= NCBIXML.parse(open(blastres_f_path))

	for blast_record in blast_records: # Each blast record is the total res for a specific query
		qorg = blast_record.query
		if qorg not in probe_table_res:
			probe_table_res[qorg] = {}
		for alignment in blast_record.alignments: # One aln for each query - subject hit
			subjct_probe = alignment.hit_def # alignment.hit_def is the subject name
			if subjct_probe not in probe_table_res[qorg]:
				if cancer:
					probe_table_res[qorg][subjct_probe] = {
					"Query nucleotide": "X",
					"Query position": 1,
					"E-value": 1
					} 
				else:
					probe_table_res[qorg][subjct_probe] = ["X",1] 
				# First value is going to be qnucl aligned at center of probe
				# Second value is going to be the evalue of hsp
			
			# alignment.length is the total aln length?
			for hsp in alignment.hsps: # Each query - subject hit has multiple hsps (blast hits)
				qstart =  hsp.query_start
				sstart =  hsp.sbjct_start
				evalue = hsp.expect 
				hsp_len = hsp.align_length
				if cancer:
					if evalue < probe_table_res[qorg][subjct_probe]["E-value"]:
						tmp = _get_aln_nucl(hsp.query,hsp.sbjct,sstart,probe_len)
						query_pos = qstart + tmp[1]  # Probe middle nucl pos in query
						probe_table_res[qorg][subjct_probe]["Query nucleotide"] = tmp[0]
						probe_table_res[qorg][subjct_probe]["Query position"] = query_pos
						probe_table_res[qorg][subjct_probe]["E-value"] = float('{:0.3e}'.format(evalue))
				else:
					if evalue < probe_table_res[qorg][subjct_probe][1]:
						tmp = _get_aln_nucl(hsp.query,hsp.sbjct,hsp.sbjct_start,probe_len)
						probe_table_res[qorg][subjct_probe][0] = tmp[0]
						probe_table_res[qorg][subjct_probe][1] = evalue
	if cancer:
		tmpdict = {}
		for k, item in probe_table_res.items():
			v = pd.DataFrame.from_dict(item,orient="index")
			v.index.name = "SNP"
			tmpdict[k] = v
		probedf = pd.concat(tmpdict,axis=0)
		probedf.index.names = ["Sequences", "SNP"]
		probedf = probedf.reset_index()
		probedf.to_excel(probe_results_path,na_rep="X")
	else:
		for qorg in probe_table_res:
			for subjct_probe in probe_table_res[qorg]:
				probe_table_res[qorg][subjct_probe] = probe_table_res[qorg][subjct_probe][0]
		probedf = pd.DataFrame.from_dict(probe_table_res,orient='index')
		probedf_cols = list(probedf.columns)
		probedf_cols = _sort_alphanumeric(probedf_cols)
		probedf = probedf[probedf_cols]
		probedf.index.name = "Sequences"
		probedf.to_excel(probe_results_path,na_rep="X")
	return probe_results_path

def find_recombinants(blastdf, exe=True):
	recombinants = []
	orgs = np.unique(blastdf["qaccver"])
	blastdf["Lineage"] = blastdf["saccver"]
	blastdf["Lineage"] = blastdf["Lineage"].apply(lambda x: x.split("_")[0][0])
	for org in orgs:
		tmpl = np.unique(blastdf[blastdf["qaccver"] == org]["Lineage"].values)
		if len(tmpl) > 1:
			# If more than one lineages exist it is recombinant
			recombinants.append(org)
	# TODO: Implement recombination status
	# I have to use results from the lineageSnps too
	return recombinants

def filter_nonHPV16(blastres_f_path, params, exe=True, filter="HPV16"):
	"""
	This function was created based on the minireview recommended on ICTV 2018 for the different types
	of alphapapillomavirus 9 species
	https://www.microbiologyresearch.org/content/journal/jgv/10.1099/jgv.0.001105
	"""
	# if not exe:
	# 	geneID_results_path = blastres_f_path.with_name(blastres_f_path.stem + "_GeneID.xlsx")
	# 	return geneID_results_path
	
	# def _get_db_lengths(db_path):
	# 	dblengths = {}
	# 	seqrecords = SeqIO.parse(db_path,"fasta")
	# 	for seqrecord in seqrecords:
	# 		name = seqrecord.id
	# 		seqlen = len(str(seqrecord.seq))
	# 		dblengths[name] = seqlen
	# 	return dblengths
	# logging.info("Identifying HPV16 sequences ")
	# db_path = pathlib.Path(params.loc["GenesRef_database","Value"])
	# dblengths = _get_db_lengths(db_path)

	blast_records= NCBIXML.parse(open(blastres_f_path))
	blast_res = {}
	for blast_record in blast_records: # Each blast record is the total res for a specific query
		qorg = blast_record.query
		if qorg not in blast_res:
			blast_res[qorg] = {
				"qaccver":qorg,
				"saccver":None,
				"qstart":0,
				"qend":0,
				"sstart":0,
				"send":0,
				"length":0,
				"pident":0,
				"evalue":1,
				"subject_covhsp":0,
			}
		for alignment in blast_record.alignments: # One aln for each query - subject hit
			subjct = alignment.hit_def # alignment.hit_def is the subject name
			subjct_name = subjct.split("_")[0]

			# alignment.length is the total aln length?
			for hsp in alignment.hsps: # Each query - subject hit has multiple hsps (blast hits)
				if evalue < blast_res[qorg]["evalue"]: # Keep only the best hsp result
					blast_res[qorg]["saccver"] = subjct
					blast_res[qorg]["qstart"] = hsp.query_start
					blast_res[qorg]["qend"] = hsp.query_end
					blast_res[qorg]["sstart"] = hsp.sbjct_start 
					blast_res[qorg]["send"] = hsp.sbjct_end
					blast_res[qorg]["evalue"] = hsp.expect
					blast_res[qorg]["length"] = hsp.align_length
					blast_res[qorg]["pident"] = hsp.identities / hsp.align_length # As calculated in NCBI blast2seq
					# blast_res[qorg]["subject_covhsp"] = round((hsp_len/dblengths[subjct])*100,4)

	df_list = []
	for k in blast_res:
		tmpdf = pd.DataFrame.from_dict(blast_res[k],orient='index')
		df_list.append(tmpdf)

	blastdf = pd.concat(df_list)
	hpvID_results_path = blastres_f_path.with_name(blastres_f_path.stem + "_HPV16_ID.xlsx")
	blastdf.to_excel(hpvID_results_path,index=False)
	hpv16_seqs = list(blastdf[blastdf["saccver"] == filter].index)
	print(hpv16_seqs)
	#TODO: Finish this
	return hpvID_results_path
