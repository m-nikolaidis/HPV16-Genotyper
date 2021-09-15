import re
import os
import logging
import pathlib
import numpy as np
import pandas as pd
from collections import Counter
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

def blastn_search(exe_params: pd.DataFrame, query_f_path: pathlib.Path, 
	workflow: str, makeblastdb_bin: str, blastn_bin: str, exe: bool = True
	) -> pathlib.Path:
	"""
	Execute blast search for the various workflows
	"""
	threads = exe_params.loc['num_threads',"Value"]
	outdir = pathlib.Path(exe_params.loc['out',"Value"])

	if workflow == "HPV16 filter":
		logging.info("Starting - BLASTn search to identify non HPV16 sequences")
		db_path = pathlib.Path(exe_params.loc['HPV_filter_db',"Value"])
		db_file = db_path.name
		blastres_f = query_f_path.stem + "_HPV16_Blastn_res.txt"
		evalue = exe_params.loc['Gene_identification_Evalue',"Value"]
		word_size = exe_params.loc['Gene_identification_Word_size',"Value"]
	if workflow == "gene identification":
		logging.info("Starting - BLASTn search for gene identification")
		db_path = pathlib.Path(exe_params.loc['GenesRef_database',"Value"])
		db_file = db_path.name
		blastres_f = "GeneIdentification_Blastn_results.txt"
		evalue = exe_params.loc['Gene_identification_Evalue',"Value"]
		word_size = exe_params.loc['Gene_identification_Word_size',"Value"]
	if workflow == "snp":
		logging.info("Starting - BLASTn search for Lineage specific SNPs")
		db_path = pathlib.Path(exe_params.loc['SNP_db_path',"Value"])
		db_file = db_path.name
		blastres_f = "lineageSpecificProbes_Blastn_results.txt"
		evalue = exe_params.loc["SNP_identification_Evalue","Value"]
		word_size = exe_params.loc["SNP_identification_Word_Size","Value"]

	if workflow == "cancer":
		logging.info("Starting - BLASTn search for increased cancer risk SNPs")
		db_path = pathlib.Path(exe_params.loc['cSNP_db_path',"Value"])
		db_file = db_path.name
		blastres_f = "cancerProbes_Blastn_results.txt"
		evalue = exe_params.loc["SNP_identification_Evalue","Value"]
		word_size = exe_params.loc["SNP_identification_Word_Size","Value"]
	
	tmp = pathlib.Path(".tmp")
	dbout = outdir / tmp / db_file
	blastres_f_path = outdir / "BlastResults" / blastres_f
	outfmt = "6 qaccver saccver qstart qend sstart send evalue bitscore length pident qcovs qcovhsp slen"
	makedb = NcbimakeblastdbCommandline(cmd = makeblastdb_bin,
										dbtype ="nucl",
										input_file = db_path,
										out = dbout
	)
	blastn_search = NcbiblastnCommandline(cmd = blastn_bin,
											query = query_f_path,
											db = dbout,
											outfmt = outfmt,
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

def _sort_alphanumeric(iteratable: list) -> list:
	int_convert = lambda text: int(text) if text.isdigit() else text 
	sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(iteratable, key = sorting_key)

def parse_GeneID_results(blastres_f_path: pathlib.Path, 
		params: pd.DataFrame, 
		exe:bool =True
	) -> pathlib.Path:
	"""
	In order to avoid the small gene fragments and to filter the genes that
	have many Ns (more than one hsp usually), a subject coverage filter is applied per hsp (>=90%)
	Writes excel file
	"""
	
	geneID_results_path = blastres_f_path.parent / ".." / "GeneIdentification_results.xlsx"
	if not exe:
		return geneID_results_path
	
	logging.info("Parsing - gene identification results ")
	
	col_names = ["Query sequence", "Subject sequence", "Query start", "Query end",
			 "Subject start", "Subject end", "E-value", "Bitscore", "Aln length", 
			 "Perc. identity", "Query coverage", "Query coverage hsp", "Subject length"
	]
	df = pd.read_csv(blastres_f_path, sep="\t",names = col_names)
	df["Subject coverage"] = ((df["Subject end"] - df["Subject start"] + 1)/ df["Subject length"] ) * 100
	df = df[df["Subject coverage"] >= 90.0]
	df["Gene"] = df["Subject sequence"].apply(lambda x: x.split("_")[2])
	df["Lineage"] = df["Subject sequence"].apply(lambda x: x.split("_")[0][0])
	
	def _filterLineage(x):
		tmpLineageList = list(np.unique(x.values))
		concatLineages = "/".join(tmpLineageList)
		return concatLineages

	def _filterSublineage(x):
		tmpSublineageList = list(np.unique(x.values))
		tmpSublineageList = [org.split("_")[0] for org in tmpSublineageList]
		concatSublineages = "/".join(tmpSublineageList)
		return concatSublineages
	
	t = df.groupby(["Query sequence","Gene"])
	idx = t["Bitscore"].transform(max) == df["Bitscore"]
	df = df[idx]
	t = df.groupby(["Query sequence","Gene"])
	df["Lineage"] = t["Lineage"].transform(lambda x: _filterLineage(x))
	df["Sublineage"] = t["Subject sequence"].transform(lambda x: _filterSublineage(x))
	df["Subject sequence"] = t["Subject sequence"].transform(lambda x: _filterLineage(x))
	df = df.drop_duplicates()	
	
	df.index = range(len(df))
	df["Perc. identity"] = df["Perc. identity"].apply(lambda x: "{:0.2f}".format(x))
	df["E-value"] = df["E-value"].apply(lambda x: F"{x:.2e}")
	
	df = df[[
		"Query sequence",
		"Gene",
		"Lineage",
		"Sublineage",
		"Query start",
		"Query end",
		"Subject sequence",
		"Subject start",
		"Subject end",
		"Aln length",
		"Perc. identity",
		"E-value"
	]]
	df.to_excel(geneID_results_path, index = False)
	return geneID_results_path, df

def parse_SNP_results(blastres_f_path: pathlib.Path, seqindex: dict,
	exe:bool = True, cancer:bool = False
	) -> pathlib.Path:
	"""
	Cancer parameter is used to differentiate from Lineage specific SNPs and cancer SNPs
	"""
	if cancer:
		probe_results_path = blastres_f_path.parent / ".." / "cancerSNP_results.xlsx"
		logging.info("Parsing - increased cancer risk SNP results ")
	else:
		probe_results_path = blastres_f_path.parent / ".." / "LineageSpecificSNPs.xlsx"
		logging.info("Parsing - lineage specific SNP results ")
	if not exe:
		return probe_results_path
	
	def _get_aln_nucl(qseq, qstart, sstart, probe_len) -> list:
		"""
		Return the query nucleotide that has hit the center of the probe
		"""
		# Round so in case of odd number, i get int
		middle_pos = (round(probe_len/2) - sstart)
		return 	[qseq[qstart + middle_pos - 1]  , middle_pos] #  - 1 because the sequence is 0 based

	probe_len = 31
	probe_table_res = {}
	col_names = ["Query sequence", "Subject sequence", "Query start", "Query end",
			 "Subject start", "Subject end", "E-value", "Bitscore", "Aln length", 
			 "Perc. identity", "Query coverage", "Query coverage hsp", "Subject length"
	]
	df = pd.read_csv(blastres_f_path, sep = "\t", names = col_names)
	
	indeces = df.index
	for idx in indeces:
		qorg = df.loc[idx, "Query sequence"]
		subjct_start = df.loc[idx, "Subject start"]
		subjct_probe = df.loc[idx, "Subject sequence"]
		qstart = df.loc[idx, "Query start"]
		if qorg not in probe_table_res:
			probe_table_res[qorg] = {}
		if subjct_probe not in probe_table_res[qorg]:
			if cancer:
				probe_table_res[qorg][subjct_probe] = {
					"Query nucleotide": "X",
					"Query position": 1,
					"E-value": 1
				}
			else:
				probe_table_res[qorg][subjct_probe] = ["X",1]
		evalue = df.loc[idx, "E-value"]
		qseq = str(seqindex[qorg].seq).upper()
		if cancer:
			if evalue < probe_table_res[qorg][subjct_probe]["E-value"]:
				tmp = _get_aln_nucl(qseq, qstart, subjct_start, probe_len)
				query_pos = qstart + tmp[1]  # Probe middle nucl pos in query
				probe_table_res[qorg][subjct_probe]["Query nucleotide"] = tmp[0]
				probe_table_res[qorg][subjct_probe]["Query position"] = query_pos
				probe_table_res[qorg][subjct_probe]["E-value"] = float('{:0.3e}'.format(evalue))
		else:
			if evalue < probe_table_res[qorg][subjct_probe][1]:
				tmp = _get_aln_nucl(qseq, qstart, subjct_start, probe_len)
				probe_table_res[qorg][subjct_probe][0] = tmp[0]
				probe_table_res[qorg][subjct_probe][1] = evalue
	if cancer:
		tmpdict = {}
		for k, item in probe_table_res.items():
			v = pd.DataFrame.from_dict(item,orient="index")
			v.index.name = "SNP"
			tmpdict[k] = v
		probedf = pd.concat(tmpdict,axis=0)
		probedf.index.names = ["Query sequence", "SNP"]
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
	return probe_results_path, probedf

def find_recombinants(blastdf: pd.DataFrame, lineageSnpDf: pd.DataFrame, 
	outdir: pathlib.Path
) -> dict:
	"""
	Use the cut-offs to identify putative recombinants / artifacts in the analysis
	Gene identificatio cut-off: Atleast 1 gene with different main lineage
	Lineage specific SNPs: 3 or more SNPs other than the dominant
	"""
	fout = outdir / "putativeRecombinants.txt"
	recombStatus = {}
	orgs = np.unique(blastdf["Query sequence"])
	for org in orgs:
		recombStatus[org] = {"geneID":0,"lSNP":0}
		# Recombination by gene identification
		tmpdf = blastdf[blastdf["Query sequence"] == org]
		tmplineages = Counter(tmpdf["Lineage"])
		dom_lineage = max(tmplineages,key=lambda key: tmplineages[key])
		for tmplineage in tmplineages:
			if dom_lineage not in tmplineage:
				recombStatus[org]["geneID"] = 1
		tmpdf = lineageSnpDf[lineageSnpDf.index == org].tail(1).T.head(67) # lineageSnpDf["Index"] will become lineageSnpDf.index
		tmplineages = dict(Counter(tmpdf[tmpdf.columns[0]].values))
		dom_lineage = max(tmplineages,key=lambda key: tmplineages[key])
		indeces = tmpdf.index
		non_dom_consec = 0
		for idx in indeces:
			val = tmpdf.loc[idx].values[0]
			if val == dom_lineage:
				non_dom_consec = 0
			if dom_lineage == "Lin_A":
				if val != dom_lineage and val != "Other":
					non_dom_consec += 1
					if non_dom_consec >= 3:
						recombStatus[org]["lSNP"] = 1	
			else:
				if dom_lineage == "Other":
					continue
				else:
					if val == "Lin_A" and val != "Other":
						non_dom_consec += 1
						if non_dom_consec >= 3:
							recombStatus[org]["lSNP"] = 1
					
	recombinants = [org for org in recombStatus if recombStatus[org]["geneID"] != 0 and recombStatus[org]["lSNP"] != 0]
	fout.write_text(os.linesep.join(recombinants))
	

def identify_nonHPV16(blastres_f_path: pathlib.Path) -> list:
	"""
	In Mirabello et al., 2018 doi:10.3390/v10020080 it is reviewed that HPV types in the same species 
	(HPV16 vs 31; which are both Alphapapillomavirus-9 'are defined by difference of atleast 10% in L1 nucleotide sequence'
	Also
	within each of these HPV types there are variant lineages and sublineages with intratypic genome sequence differences of 1.0–10% and 0.5–1.0%,
	respectively
	"""

	logging.info("Identifying HPV16 sequences")

	col_names = ["Query sequence", "Subject sequence", "Query start", "Query end", "Subject start", "Subject end",
			"E-value", "Bitscore", "Aln length", "Perc. identity", "Query coverage", "Query coverage hsp", "Subject length"
	]
	df = pd.read_csv(blastres_f_path, sep="\t",names = col_names)
	hpv16filt_results_path = blastres_f_path.parent / ".." / "HPV16filt_BlastN_results.xlsx"
	kept_indeces = []
	seqs = np.unique(df["Query sequence"])
	for seq in seqs:
		tmpdf = df.loc[(df["Query sequence"] == seq ) & (df["Perc. identity"] >= 90.00) & (df["Query coverage"] >= 78) & (df["Subject sequence"] == "HPV16_cg")] # 78% was used in our relaxed analysis (45 recombination events), to allow many Ns
		if not tmpdf.empty:
			evalue = 1
			final_idx = None
			indeces = tmpdf.index
			for idx in indeces:
				if tmpdf.loc[idx, "E-value"] < evalue:
					evalue = tmpdf.loc[idx, "E-value"]
					final_idx = idx
			kept_indeces.append(final_idx)
	df = df.loc[kept_indeces]
	hpv16_sequences = list(df["Query sequence"])
	df.to_excel(hpv16filt_results_path)

	logging.info(F"Finished, identified {len(hpv16_sequences)} HPV16 sequences ")
	return hpv16_sequences