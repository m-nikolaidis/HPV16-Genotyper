import re
import pandas as pd
import pathlib
import numpy as np
from collections import Counter

def create_annot_file(f, annot_f, prefix="Lin", empty="Other"):
	"""
	TODO: Write me
	"""
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
	
def annotate_results(annot_f,probe_res,empty="Other",exe=True):
	"""
	Probe results is the file that has mapped each nucleotide results of the probe blast
	For the non ATGC character value will be automatically Other
	"""
	xl_engine = "openpyxl"
	annotdf = pd.read_excel(annot_f,index_col=0, engine=xl_engine)
	probedf = pd.read_excel(probe_res,index_col=0, engine=xl_engine)

	if not exe:
		num_snps = len(probedf.columns)
		print(num_snps)
		return num_snps

	probedf_cp = probedf.copy()
	# TODO: Explain what you are doing and why

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
	outdir = probe_res.parent
	outf = outdir / "SNPs_results.xlsx"
	probedf.to_excel(outf)

	# Create new probedf
	# with the Lineages instead of the Nucleotides
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
	tmpdf.index.name = "Index"
	tmpdf.to_excel(probe_res, engine=xl_engine)
	
	return num_snps

#TODO: Delete after testing
if __name__ == "__main__":
	annot_f = pathlib.Path("/media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/resources/SNPs_annotation.xlsx")
	probe_res = pathlib.Path("/media/marios/Data/Lab/HPV16/HPV16_genotyping_tool/src/mock_2021_08_24/Filt_mock_Probe_Blastn_res_SNP_nucl.xlsx")
	annotate_results(annot_f=annot_f,probe_res=probe_res)