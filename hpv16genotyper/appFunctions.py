# from copy import error
import re
import os
import sys
import shlex
import pathlib
import logging
import subprocess
import multiprocessing
import pandas as pd
import numpy as np
from multiprocessing.pool import ThreadPool
from collections import Counter
from PyQt5.QtCore import pyqtSignal, QObject
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

pd.options.mode.chained_assignment = None

class MainFunctions(QObject):
	countSignal = pyqtSignal(int)
	pBarIdxSignal = pyqtSignal(int)
	errorSignal = pyqtSignal(str)
	startTimerSignal = pyqtSignal(int)
	stopTimerSignal = pyqtSignal(int)
	def __init__(self):
		super().__init__()
	
	def emitSignal(self, updateValue: int = None, barIdx: int = None, errorMsg: str = None ) -> None:
		"""
		Use this signal to update the progress bars in main Window ui
		The update step is going to be the current value of the progress bar
		"""
		if barIdx != None:
			self.pBarIdxSignal.emit(barIdx)
		if updateValue != None:
			self.countSignal.emit(updateValue)
		if errorMsg != None:
			self.errorSignal.emit(errorMsg)

	def startTimer(self) -> None:
		self.startTimerSignal.emit(1)
	
	def stopTimer(self) -> None:
		self.stopTimerSignal.emit(1)
	
	##### General functions
	def _sort_alphanumeric(self, iteratable: list) -> list:
		int_convert = lambda text: int(text) if text.isdigit() else text 
		sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
		return sorted(iteratable, key = sorting_key)

	##### Set-up environment
	def create_dirs(self, outdir : pathlib.Path, exist_ok: bool = True) -> pathlib.Path:
		"""
		Perform the necessary actions to set up the script working environment
		"""
		if re.match(r".+_run\d+$",outdir.name):
			tmpname = outdir.name.split("_run")[0]
			parentDir = outdir.parent
			outdirs = os.listdir(parentDir)
			if len(outdirs) != 0:
				outdirs = [newdir for newdir in outdirs if tmpname in newdir]
				outdirs = self._sort_alphanumeric(outdirs)
				outdir  = pathlib.Path(outdirs[-1])
				dirname = outdir.name
				tmp = dirname.split("_run")
				newrun = int(tmp[1]) + 1
				newname = "".join([tmp[0], "_run", str(newrun)])
				outdir = parentDir / newname

		script_out = outdir
		tmp_out = script_out / pathlib.Path(".tmp")
		tmp_genes_out = tmp_out / pathlib.Path("Gene_seqs") # Gene files of batch size
		blast_results = script_out /  pathlib.Path("BlastResults")
		profiles_out = script_out / pathlib.Path("Alignments") # Final profile alns for each gene
		trees_out = script_out / pathlib.Path("Phylogenetic_Trees")
		graphics_out = script_out / pathlib.Path("Graphics")
		gene_sim_graphics = graphics_out / pathlib.Path("GeneIdentification")
		snp_graphics = graphics_out / pathlib.Path("LineageSpecificSNPs")
		
		script_out.mkdir(exist_ok = exist_ok, parents=True)
		tmp_out.mkdir(exist_ok = exist_ok)
		tmp_genes_out.mkdir(exist_ok = exist_ok)
		blast_results.mkdir(exist_ok = exist_ok)
		profiles_out.mkdir(exist_ok = exist_ok)
		trees_out.mkdir(exist_ok = exist_ok)
		graphics_out.mkdir(exist_ok = exist_ok)
		gene_sim_graphics.mkdir(exist_ok = exist_ok)
		snp_graphics.mkdir(exist_ok = exist_ok)
		return script_out


	def file_exists(self, f: pathlib.Path) -> str or None:
		exist = f.exists()
		if exist == False:
			print(F"ERROR:CRITICAL: {f} does not exist")
			msg = F" Input file {f} does not exist "
			self.emitSignal(errorMsg = msg)
			raise Exception(msg)
		logging.info(F"{f} exists")
		return None
	
	def input_file_format(self, fasta_file_path: pathlib.Path) -> str or None:
		"""
		Check if the input query file is in supported format
		Currently supporting formats: Fasta
		"""
		check_handle = open(fasta_file_path, "r")

		if fasta_file_path.stat().st_size == 0:
			print("ERROR:CRITICAL: Input file is empty, check input file ")
			msg = " Input file is empty, please check input file and restart the application"
			self.emitSignal(errorMsg = msg)
			return msg
		try:
			first_line = check_handle.readline().rstrip()
			if re.match(r'^>', first_line): # Is fasta
				pass
			else:
				print("ERROR:CRITICAL: The input file is not fasta ")
				msg = " The input file is not fasta "
				self.emitSignal(errorMsg = msg)
				return msg
		except UnicodeError:
			print("ERROR:CRITICAL: The input file is not in plain text file format ")
			msg = " The input file is not in plain text file format, please check input file and restart the application"
			self.emitSignal(errorMsg = msg)
			return msg
		logging.info(F"Finished input file integrity check successfully")
		return None

	def dir_priviledges(self, path: pathlib.Path) -> str or None:
		"""
		Check if the user has write priviledges in the specified directory
		Input: Posix | Window path from pathlib library
		"""
		write_priviledges = os.access(path.parent, os.W_OK)
		if not write_priviledges:
			print(F"ERROR:CRITICAL: Program has no write permission on {path.parent}")
			msg = "You cannot write in the specified output directory, please check output directory and restart the application"
			self.emitSignal(errorMsg = msg)
			return msg
		logging.info(F"Have permission to write on {path}")
		return None

	def filter_hpv16(self, fasta_file_path: pathlib.Path, 
		seqs_to_include: list, outdir: pathlib.Path
	) -> str:
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
	
	##### BLAST search
	def blastProgress(self):
		column_names = ["qaccver","saccver","qstart","qend","sstart","send","evalue","bitscore","length","pident","qcovs","qcovhsp","slen"]
		def _helper(f):
			done = 0
			if f.exists():
				df_readCSVPartial = pd.read_csv(f,sep="\t", names = column_names, usecols = ["qaccver"]) # Read partial data for lower memory consuption
				df_readCSVPartial["qaccver"] = df_readCSVPartial["qaccver"].astype("category") # Lower memory consuption for large datasets
				done = len(np.unique(df_readCSVPartial["qaccver"]))
			return done
		for i in self.blastPool.imap(_helper, [self.blastres_f_path]):
			done = i
			done_perc = int((done/self.total)*100)
			self.emitSignal(done_perc, self.barIdx)	

	def blastn_search(self, exe_params: pd.DataFrame, query_f_path: pathlib.Path, 
			workflow: str, makeblastdb_bin: str, blastn_bin: str, total: int
		) -> pathlib.Path:
		"""
		Execute blast search for the various workflows
		"""
		threads = exe_params.loc['num_threads',"Value"]
		outdir = pathlib.Path(exe_params.loc['out',"Value"])

		if workflow == "HPV16 filter":
			barIdx = 0
			logging.info("Starting - BLASTn search to identify non HPV16 sequences")
			db_path = pathlib.Path(exe_params.loc['HPV_filter_db',"Value"])
			db_file = db_path.name
			blastres_f = query_f_path.stem + "_HPV16_Blastn_res.txt"
			evalue = exe_params.loc['Gene_identification_Evalue',"Value"]
			word_size = exe_params.loc['Gene_identification_Word_size',"Value"]
		if workflow == "gene identification":
			barIdx = 1
			logging.info("Starting - BLASTn search for gene identification")
			db_path = pathlib.Path(exe_params.loc['GenesRef_database',"Value"])
			db_file = db_path.name
			blastres_f = "GeneIdentification_Blastn_results.txt"
			evalue = exe_params.loc['Gene_identification_Evalue',"Value"]
			word_size = exe_params.loc['Gene_identification_Word_size',"Value"]
		if workflow == "snp":
			barIdx = 2
			logging.info("Starting - BLASTn search for Lineage specific SNPs")
			db_path = pathlib.Path(exe_params.loc['SNP_db_path',"Value"])
			db_file = db_path.name
			blastres_f = "lineageSpecificProbes_Blastn_results.txt"
			evalue = exe_params.loc["SNP_identification_Evalue","Value"]
			word_size = exe_params.loc["SNP_identification_Word_Size","Value"]
		if workflow == "cancer":
			barIdx = 3
			logging.info("Starting - BLASTn search for increased cancer risk SNPs")
			db_path = pathlib.Path(exe_params.loc['cSNP_db_path',"Value"])
			db_file = db_path.name
			blastres_f = "cancerProbes_Blastn_results.txt"
			evalue = exe_params.loc["SNP_identification_Evalue","Value"]
			word_size = exe_params.loc["SNP_identification_Word_Size","Value"]
		tmp = pathlib.Path(".tmp")
		dbout = outdir / tmp / db_file
		self.blastres_f_path = outdir / "BlastResults" / blastres_f
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
												out = self.blastres_f_path,
												evalue = evalue,
												num_threads = threads,
												word_size = word_size,
												dust="no"
		)
	
		makedb()
		self.blastPool = ThreadPool(1)
		self.total = total
		self.barIdx = barIdx
		self.startTimer() # Send signal to start the blastTimer
		blastn_search()
		self.blastPool.close()
		self.blastPool.join()
		self.stopTimer() # Send signal to main app to stop the timer
		self.emitSignal(100, barIdx)
		logging.info("Finished - BLASTn search ")
		return self.blastres_f_path

	def parse_GeneID_results(self, blastres_f_path: pathlib.Path, 
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

	def parse_SNP_results(self, blastres_f_path: pathlib.Path, seqindex: dict,
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

		def _cancerSNP(snpName):
			"""
			Return the cancer mutation depending on the queried reference SNP
			"""
			if snpName == "NC001526_NuclPos_145_G":
				return "T"
			if snpName == "NC001526_NuclPos_335_C":
				return "T"
			if snpName == "NC001526_NuclPos_350_T":
				return "G"
			if snpName == "NC001526_NuclPos_647_A":
				return "G"
			if snpName == "NC001526_NuclPos_749_C":
				return "T"
			if snpName == "NC001526_NuclPos_2860_C":
				return "A"
			if snpName == "NC001526_NuclPos_3410_C":
				return "T"
			if snpName == "NC001526_NuclPos_3684_C":
				return "A"
			if snpName == "NC001526_NuclPos_4042_A":
				return "G/C/T"

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
			probedf.index.names = ["Query sequence", "Reference genome position"]
			probedf = probedf.reset_index()
			probedf["Cancer SNP"] = probedf["Reference genome position"].apply(lambda x: _cancerSNP(x))
			probedf["Reference genome position"] =  probedf["Reference genome position"].apply(lambda x: "_".join(x.split("_")[:3]))
			probedf = probedf[["Query sequence", "Reference genome position", 
				"Cancer SNP", "Query nucleotide", "Query position", "E-value"]
			]
			probedf.to_excel(probe_results_path,na_rep="X")
		else:
			for qorg in probe_table_res:
				for subjct_probe in probe_table_res[qorg]:
					probe_table_res[qorg][subjct_probe] = probe_table_res[qorg][subjct_probe][0]
			probedf = pd.DataFrame.from_dict(probe_table_res,orient='index')
			probedf_cols = list(probedf.columns)
			probedf_cols = self._sort_alphanumeric(probedf_cols)
			probedf = probedf[probedf_cols]
			probedf.index.name = "Sequences"
			probedf.to_excel(probe_results_path,na_rep="X")
		return probe_results_path

	def identify_nonHPV16(self, blastres_f_path: pathlib.Path) -> list:
		"""
		In Mirabello et al., 2018 doi:10.3390/v10020080 it is reviewed that HPV types in the same species 
		(HPV16 vs 35; which are both Alphapapillomavirus-9 'are defined by difference of atleast 10% in L1 nucleotide sequence'
		Also
		within each of these HPV types there are variant lineages and sublineages with intratypic genome sequence differences of 1.0–10% and 0.5–1.0%,
		respectively
		HPV35 is the closest relative to HPV16 based on fig1 of 10.1371/journal.pone.0020183
		"""

		logging.info("Identifying HPV16 sequences")

		col_names = ["Query sequence", "Subject sequence", "Query start", "Query end", "Subject start", "Subject end",
				"E-value", "Bitscore", "Aln length", "Perc. identity", "Query coverage", "Query coverage hsp", "Subject length"
		]
		df = pd.read_csv(blastres_f_path, sep="\t",names = col_names)
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

		nonHPV16Seqs = list(set(seqs).difference(hpv16_sequences))
		nonHPV16fout = blastres_f_path.parent.parent / "nonHPV16.txt"
		nonHPV16fout.write_text(os.linesep.join(nonHPV16Seqs))
		
		logging.info(F"Finished, identified {len(hpv16_sequences)} HPV16 and {len(nonHPV16Seqs)} non HPV16 sequences ")
		return hpv16_sequences

	def annotate_results(self, annot_f: pathlib.Path,probe_res: pathlib.Path,
		empty: str ="Other",exe: bool=True
		) -> int:
		"""
		Probe results is the file that has mapped each nucleotide results of the probe blast
		For the non ATGC character value will be automatically Other
		"""
		xl_engine = "openpyxl"
		annotdf = pd.read_excel(annot_f,index_col=0, engine=xl_engine)
		probedf = pd.read_excel(probe_res,index_col=0, engine=xl_engine)
		to_replace = {"Lin_ABD":"Other", "Lin_ABC":"Other", 
			"Lin_ABD":"Other", "Lin_ACD":"Other"
		}
		if not exe:
			num_snps = len(probedf.columns)
			return num_snps
		# probedf_cp = probedf.copy()

		ref_prefix_regex = re.compile(r'(\S+)_\d+$') # Take this based on the probedf
		probedf_first_col = probedf.columns[0]
		m = re.match(ref_prefix_regex,probedf_first_col)
		ref_prefix = m.group(1)
		tmpl = annotdf.index.to_list()
		tmpl = [ref_prefix + "_" + str(l) for l in tmpl]
		annotdf.index = tmpl # So i can match the columns of probe results file

		annotdict = annotdf.to_dict(orient='index')

		possible_nucl_chars=list(np.unique(probedf.values)) # Non ATGC characters will be automatically value to empty

		for char in possible_nucl_chars:
			if char == "A" or char == "T" or char == "G" or char == "C":
				pass
			else:
				for probe in annotdict:
					annotdict[probe][char] = empty
		probedf = probedf.replace(annotdict)
		probedf = probedf.replace(to_replace)
		num_snps = len(probedf.columns)
		counters = {}
		indeces = probedf.index
		for idx in indeces:
			 counters[idx] = dict(Counter(probedf.loc[idx]))
		t = pd.DataFrame.from_dict(counters,orient='index')
		t.fillna(0.0,inplace=True)
		t = t.apply(lambda x:round((x/num_snps)*100,2))
		cols = t.columns.to_list()
		cols.remove(empty)
		indeces = t.index
		for idx in indeces:
			maxval = t.loc[idx,cols].astype(float).idxmax()
			t.loc[idx,"Dominant lineage"] = maxval
		probedf = probedf.join(t)
		annotdf = annotdf.T
		extra_row = {"Nucleotides":{}}
		for col in annotdf.columns:
			for idx in annotdf.index:
				if col not in extra_row["Nucleotides"]:
					extra_row["Nucleotides"][col] = idx + ":" + annotdf.loc[idx,col]
				else:
					extra_row["Nucleotides"][col] += "," + idx + ":" + annotdf.loc[idx,col]
		tmpdf = pd.DataFrame.from_dict(extra_row,orient='index')
		# tmpdf = pd.concat([tmpdf,probedf_cp])
		tmpdf = pd.concat([tmpdf,probedf])
		tmpdf.index.name = "Query sequence"
		tmpdf.fillna("X",inplace=True)
		tmpdf.to_excel(probe_res, engine=xl_engine)
		return num_snps, tmpdf

	def find_recombinants(self, blastdf: pd.DataFrame, lineageSnpDf: pd.DataFrame, 
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
			# Recombination by lineage SNPs 
			tmpdf = lineageSnpDf[lineageSnpDf.index == org].T.head(67)
			tmplineages = dict(Counter(tmpdf[tmpdf.columns[0]].values))
			dom_lineage = max(tmplineages,key=lambda key: tmplineages[key])
			if dom_lineage == "Other":
				continue
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
					if val == "Lin_A" and val != "Other":
						non_dom_consec += 1
						if non_dom_consec >= 3:
							recombStatus[org]["lSNP"] = 1
		recombinants = [org for org in recombStatus if recombStatus[org]["geneID"] != 0 and recombStatus[org]["lSNP"] != 0]
		fout.write_text(os.linesep.join(recombinants))

	def _callMultiThreadProcc(self, cmd):
		system = sys.platform
		if "win" in system:
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		else:
			p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		p.communicate()

	def prepare_alns(self, outdir: pathlib.Path, query_f_path: pathlib.Path, 
		geneid_blast_res_xls: pathlib.Path, muscle_bin: str, 
		profiledb_dir: pathlib.Path, threads: int
	) -> list:
		"""
		Align gene sequences to create profiles
		input: gene_files, list of pathlib paths
		return: list of pathlib paths pointing to each alignment file
		"""

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

		
		logging.info(F"Creating the alignment files for each gene")
		coords_df = pd.read_excel(geneid_blast_res_xls,engine="openpyxl") 
		# TODO: Check for efficiency https://pandas.pydata.org/pandas-docs/stable/user_guide/scale.html
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

		# Alignment
		aln_files = []
		total_processes = len(gene_files)
		done_processes = 0
		aln_out_dir = outdir / pathlib.Path("Alignments")
		pool = ThreadPool(threads)
		gene_regex = re.compile(r'^\S+_([E|L]\d+).fa$')
		def _alnFunc(f):
				"""
				Initialize the command to run and return 1 once the subprocess is done
				"""
				m = re.match(gene_regex,f.name)
				gene = m.group(1)
				fout = aln_out_dir / (f.stem + "_aln.fa")
				aln_files.append(fout)
				db_gene_f = profiledb_dir / pathlib.Path(gene + "_profile.fa")
				arguments = " -profile -in1 " + str(f) + \
					" -in2 " + str(db_gene_f) + " -quiet -out " + str(fout)
				
				self._callMultiThreadProcc(muscle_bin + arguments)
				return 1
		


		for i in pool.imap_unordered(_alnFunc, gene_files):
			done_processes += i
			self.emitSignal(int((done_processes/total_processes)*100), 4)	
		pool.close()
		pool.join()	
		logging.info(F"Finished aligning files")
		return aln_files

	def build_trees(self, outdir: pathlib.Path, aln_files: list, seaview_bin: str, 
		threads: int, method: str = "BioNJ", dist: str = "Kimura", 
		nj_bootstrap_repl: int = 100, exe: bool = True
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
			def _treeFunc(aln_f):
				tree_file_out = trees_dir / (aln_f.stem + "_NJ_tree.nwk")
				arguments = " -build_tree -NJ -distance " + dist \
					+ " -replicates " + str(nj_bootstrap_repl) + " -o " \
					+ str(tree_file_out) + " " + str(aln_f)
				self._callMultiThreadProcc(seaview_bin + arguments)
				return 1

			pool = ThreadPool(threads)
			nj_regex = re.compile(r"^\[NJ \d+ sites Kimura.+\] (\S+)")
			self.finishedProcesses = 0
			total_processes = len(aln_files)
			done_processes = 0
			for i in pool.imap_unordered(_treeFunc, aln_files):
				done_processes += i
				self.emitSignal(int((done_processes/total_processes)*100), 5)	
			pool.close()
			pool.join()
			for aln_f in aln_files:
				tree_file_out = trees_dir / (aln_f.stem + "_NJ_tree.nwk")
				_clean_nj_output(tree_file_out, nj_regex)
		logging.info(F"Finished")
		return trees_dir

	def main(self, paramsdf: pd.DataFrame, exe: bool = True) -> pathlib.Path:
		# Basic variables
		hpv16error = False # To stop execution if no hpv16 sequences are found
		outdir = pathlib.Path(paramsdf.loc["out","Value"])
		system = paramsdf.loc["system","Value"]
		indir = pathlib.Path(paramsdf.loc["in","Value"])
		query_f = pathlib.Path(paramsdf.loc["query","Value"])
		query_f_path = indir / query_f
		makeblastdb_bin, blastn_bin, muscle_bin, seaview_bin = _init_binaries(system)
		annot_f = pathlib.Path(paramsdf.loc["SNP_annotation_file","Value"])
		threads = paramsdf.loc["num_threads", "Value"]
		
		# Critical checks
		priviledgeMsg = self.dir_priviledges(outdir)
		existMsg = self.file_exists(query_f_path)
		fmtMsg = self.input_file_format(query_f_path)
		if priviledgeMsg != None or existMsg != None or fmtMsg != None:
			return None
		outdir = self.create_dirs(outdir)
		paramsdf.loc["out","Value"] = outdir
		
		# Grab the root logger instance and use it for logging
		logfile = outdir / pathlib.Path(".logfile.log")
		# logging.config.fileConfig(str(logfile))
		logging.getLogger().handlers.clear()
		logger = logging.getLogger()
		fhandler = logging.FileHandler(filename=logfile, encoding="UTF-8")
		formatter = logging.Formatter("%(asctime)s\t%(message)s")
		fhandler.setFormatter(formatter)
		logger.addHandler(fhandler)
		logger.setLevel(logging.DEBUG)
		
		logging.info("#### DO NOT DELETE THIS FILE! ###")
		logging.info(F"Environment set successfully in {outdir}")
		logging.info("Parameters used in this run: ")
		indeces = paramsdf.index
		for index in indeces:
			val = paramsdf.loc[index,"Value"]
			logging.info(F"param: {index},{val}")
		
		### Basic workflow
		# Filter input file
		workflow = "HPV16 filter"

		query_f_index = SeqIO.index(str(query_f_path), "fasta")
		query_orgs_num = len(query_f_index.keys())

		hpv16filt_results = self.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin, total = query_orgs_num)
		hpv16_seqs = self.identify_nonHPV16(hpv16filt_results)
		if len(hpv16_seqs) == 0:
			hpv16error = True
			return query_f_path, hpv16error
		# query_f_filt = env_setup.filter_hpv16(query_f_path, hpv16_seqs, outdir)
		query_f_filt = self.filter_hpv16(query_f_path, hpv16_seqs, outdir)
		query_f = query_f_filt
		query_f_path = outdir / query_f
		logging.info(F"param: query,{str(query_f_path)}")

		query_f_index = SeqIO.index(str(query_f_path), "fasta")
		query_orgs_num = len(query_f_index.keys())
		# BLASTn search
		workflow = "gene identification"
		geneid_blastn_res_path = self.blastn_search(paramsdf, query_f_path, workflow, 
			makeblastdb_bin, blastn_bin, query_orgs_num
		)
		geneid_blastn_res_path_xlsx, blastDF = self.parse_GeneID_results(geneid_blastn_res_path,paramsdf, exe=exe)
		
		workflow = "snp"
		snp_blastn_res_path = self.blastn_search(paramsdf, query_f_path, workflow, 
			makeblastdb_bin, blastn_bin, query_orgs_num
		)
		snp_blastn_res_path_xlsx = self.parse_SNP_results(snp_blastn_res_path, query_f_index, exe=exe)
		_, snpDFannot = self.annotate_results(annot_f,snp_blastn_res_path_xlsx, exe=exe)
		
		workflow = "cancer"
		C_snp_blastn_res_path = self.blastn_search(paramsdf, query_f_path, workflow, 
			makeblastdb_bin, blastn_bin, query_orgs_num
		)
		self.parse_SNP_results(C_snp_blastn_res_path, query_f_index, exe=exe, cancer=True)
		self.find_recombinants(blastDF, snpDFannot, outdir)

		# Gene alignments and trees
		profiledb_dir = paramsdf.loc["GenesProfile_directory","Value"]
		aln_files = self.prepare_alns(
			outdir, query_f_path, geneid_blastn_res_path_xlsx, 
			muscle_bin, profiledb_dir, threads = threads
		)
		self.build_trees(outdir, aln_files, seaview_bin, threads = threads, method = "BioNJ")
		return query_f_path, outdir, hpv16error

def _init_binaries(system: sys.platform) -> list:
	"""
	Initialize the paths for the needed binaries depending on the system
	input: System info
	return: List of string
	"""
	if "win" in system:
		makeblastdb_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "makeblastdb.exe"
		blastn_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "blastn.exe"
		muscle_bin= pathlib.Path(__file__).parent / "resources" / "Win_bin" / "muscle3.8.31.exe"
		seaview_bin = pathlib.Path(__file__).parent / "resources" / "Win_bin" / "seaview.exe"
	if "linux" in system:
		makeblastdb_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "makeblastdb"
		blastn_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "blastn"
		muscle_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "muscle3.8.31_i86linux"
		seaview_bin = pathlib.Path(__file__).parent / "resources" / "Linux_bin" / "seaview4"
	return [str(makeblastdb_bin), str(blastn_bin), str(muscle_bin), str(seaview_bin)]

def _defaultparams() -> dict:
	system = sys.platform
	input_dir = ""
	query_f = ""
	outdir = pathlib.Path(__file__).absolute().parent / pathlib.Path("Results_directory")
	threads_to_use =  multiprocessing.cpu_count() - 2
	params = {
		"system": system,
		"in": input_dir,
		"out": outdir,
		"num_threads": threads_to_use,
		"query": query_f,
		"SNP_identification_Evalue": 0.005,
		"SNP_identification_Word_Size": 4,
		"SNP_annotation_file": pathlib.Path(__file__).parent / 'resources'/ 'SNPs_annotation.xlsx',
		"SNP_db_path": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'NC_001526_probes.fa',
		"cSNP_db_path": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'NC_001526_cancer_probes.fa',
		"Gene_identification_Evalue": 1e-5, 
		"Gene_identification_Word_size": 7,
		"GenesProfile_directory": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'profiles',
		"GenesRef_database":pathlib.Path(__file__).parent / 'resources' / 'sequences' / '16refs_gene_db.fa',
		"SimplotRef_database":pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'Reference_genomes_profile_mafft_GINSI.fa',
		"HPV_filter_db": pathlib.Path(__file__).parent / 'resources' / 'sequences' / 'alphapapillomavirus9.fa'
	}
	return params
