from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import os
import logging
import re
import pathlib

def _sort_alphanumeric(iteratable):
	# Helper func Sort the given motif list alphanumerically :return: sorted list 
	int_convert = lambda text: int(text) if text.isdigit() else text 
	sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(iteratable, key = sorting_key)



def parse_probe_blast_results(blastres_f,outdir, query_f):
	"""
	Reads the blast output files and ...
	Return None
	"""
	def _get_aln_nucl(query_aln,subjct_aln,sstart,probe_len):
		# I have used this approach for the case that gaps are inserted in the subject seq
		# If not return fixed position
		# TODO CHECK IF IT IS POSSIBLE
		sbjct_pos = 0
		middle_pos = (round(probe_len/2) - sstart) + 1 # Round so in case of odd number, i get int
		# print(middle_pos,probe_len,sstart)		
		for x in range(len(subjct_aln)):
			if subjct_aln[x] != "-":
				sbjct_pos += 1
				if sbjct_pos == middle_pos:
					return (query_aln[x], sbjct_pos)
		return ["X",1]

	probe_len = 31 # TODO; should be calculated automatically or provided
	blast_res = {}
	probe_table_res = {}
	blast_records= NCBIXML.parse(open(blastres_f))

	for blast_record in blast_records: # Each blast record is the total res for a specific query
		qorg = blast_record.query
		if qorg not in blast_res:
			blast_res[qorg] = {"qaccver":qorg,
				"saccver":"",
				"qstart":0,
				"qend":0,
				"sstart":0,
				"send":0,
				"length":0,
				"pident":0,
				"mismatch":0,
				"evalue":1,
				"alignment":"",
				"Query_nucl":"X",
				"Query_genomic_pos":0
			}
		if qorg not in probe_table_res:
			probe_table_res[qorg] = {}

		for alignment in blast_record.alignments: # One aln for each query - subject hit
			subjct_probe = alignment.hit_def # alignment.hit_def is the subject name
			if subjct_probe not in probe_table_res[qorg]:
				probe_table_res[qorg][subjct_probe] = ["X",1] 
				# First value is going to be qnucl aligned at center of probe
				# Second value is going to be the evalue of hsp
			
			# alignment.length is the total aln length?
			for hsp in alignment.hsps: # Each query - subject hit has multiple hsps (blast hits)
				qstart =  hsp.query_start
				qend =  hsp.query_end
				sstart =  hsp.sbjct_start
				send =  hsp.sbjct_end
				evalue = hsp.expect 
				hsp_len = hsp.align_length
				# pident =  hsp.identities / hsp. Need to find a way to calcute percent identity correctly (blast2seq)
				hsp_aln_str =  hsp.query + "\n" + hsp.match + "\n" + hsp.sbjct
				if evalue < probe_table_res[qorg][subjct_probe][1]:
					tmp = _get_aln_nucl(hsp.query,hsp.sbjct,hsp.sbjct_start,probe_len)
					probe_table_res[qorg][subjct_probe][0] = tmp[0]
					probe_table_res[qorg][subjct_probe][1] = evalue

				if evalue < blast_res[qorg]["evalue"]:
					blast_res[qorg]["saccver"] = subjct_probe
					blast_res[qorg]["qstart"] = qstart
					blast_res[qorg]["qend"] = qend
					blast_res[qorg]["sstart"] = sstart 
					blast_res[qorg]["send"] = send
					blast_res[qorg]["evalue"] = evalue
					blast_res[qorg]["length"] = hsp_len
					# blast_res[qorg]["pident"] = 
					# blast_res[qorg]["mismatch"] = 
					blast_res[qorg]["alignment"] = hsp_aln_str
					tmp = _get_aln_nucl(hsp.query,hsp.sbjct,hsp.sbjct_start,probe_len)
					blast_res[qorg]["Query_nucl"] = tmp[0]
					blast_res[qorg]["Query_genomic_pos"] = qstart + tmp[1] - 1 # 1 because it is indexed in the genome (1 base) 
	
	blastdf = pd.DataFrame.from_dict(blast_res,orient='index')
	blastdf.to_excel(outdir + query_f + "Blastnres.xlsx",index=False)
	for qorg in probe_table_res:
		for subjct_probe in probe_table_res[qorg]:
			probe_table_res[qorg][subjct_probe] = probe_table_res[qorg][subjct_probe][0]
	probedf = pd.DataFrame.from_dict(probe_table_res,orient='index')
	return probedf

# TODO: Make evrything using path and os modules for compatibility
paramsdf = pd.read_excel("../input/params.xlsx",index_col=0)
params = paramsdf.to_dict()['Value']
outdir = params['out']
query_f = params['query']
blastn_f = outdir + query_f + "_Blastn_res.xml"
# _check_file_exists(blastn_f, outdir) # Import from _checks
# probedf = parse_probe_blast_results(blastn_f, outdir, query_f)
# probedf_cols = list(probedf.columns)
# probedf_cols = _sort_alphanumeric(probedf_cols)
# probedf = probedf[probedf_cols]
# probedf.to_excel(outdir + query_f + "_Probes.xlsx",na_rep="X")

def _get_db_lengths():
	# TODO: add the database path as input to function
	db_path = os.path.join(str(pathlib.Path(__file__)), str(pathlib.Path('../resources/sequences/16refs_gene_db.fa'))).replace(os.path.basename(__file__) + "/", "")
	# TODO: This should be given from blast script or should be updated for windows / linux
	dblengths = {}
	seqrecords = SeqIO.parse(db_path,"fasta")
	for seqrecord in seqrecords:
		name = seqrecord.id
		seqlen = len(str(seqrecord.seq))
		dblengths[name] = seqlen
	return dblengths

def parse_gene_identification_blast_results(blastres_f,outdir, query_f, dblengths):
	blast_records= NCBIXML.parse(open(blastres_f))
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
				"subject_covtotal":0,
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
					blast_res[qorg][gene]["pident"] = pident
					blast_res[qorg][gene]["subject_covhsp"] = (hsp_len/dblengths[subjct])*100
					# print(subjct,(hsp_len/dblengths[subjct])*100, qorg,hsp)
					# input()
	df_list = []
	for k in blast_res:
		tmpdf = pd.DataFrame.from_dict(blast_res[k],orient='index')
		df_list.append(tmpdf)

	blastdf = pd.concat(df_list)
	blastdf.to_excel(outdir + query_f + "Blastnres.xlsx",index=False)


dblengths = _get_db_lengths()
parse_gene_identification_blast_results(blastn_f, outdir, query_f, dblengths)