import re
import pandas as pd
import pathlib
import numpy as np
from collections import Counter

def create_annot_file(f: pathlib.Path, annot_f: pathlib.Path, 
	prefix: str="Lin", empty: str="Other"
	) -> None:

	results = {}
	df = pd.read_excel(f,index_col=0, engine="openpyxl")
	lineages = list(df.columns)
	indeces = df.index
	for idx in indeces:
		if idx not in results:
			results[idx] = {"A":"","T":"","G":"","C":""}
		for lin in lineages:
			if type(df.loc[idx,lin]) == str:
					nucl = df.loc[idx,lin]
			else:
					nucl = df.loc[idx,lin].head(1).values[0]
			if lin not in results[idx][nucl]:
				results[idx][nucl] += lin
	for k in results:
		for nucl in results[k]:
			if results[k][nucl] == "":
				results[k][nucl] = empty
			else:
				results[k][nucl] = prefix + "_" + str(results[k][nucl])
	resdf = pd.DataFrame.from_dict(results,orient='index')
	resdf.to_excel(annot_f)

def annotate_results(annot_f: pathlib.Path,probe_res: pathlib.Path,
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
	probedf_cp = probedf.copy()

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
	tmpdf = pd.concat([tmpdf,probedf_cp])
	tmpdf = pd.concat([tmpdf,probedf])
	tmpdf.index.name = "Index"
	tmpdf.fillna("X",inplace=True)
	tmpdf.to_excel(probe_res, engine=xl_engine)
	return num_snps