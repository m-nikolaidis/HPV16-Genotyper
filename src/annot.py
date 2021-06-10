import pandas as pd
import numpy as np
from collections import Counter

def create_annot_file(f, annot_f, prefix="Lin", empty="Other"):
	"""
	"""
	results = {}
	lineages = ["A","B","C","D"] # TODO Should take this from the excel, so we can do sublineage identification?
	# Write in MD that the excel must have the lineage names in the columns (e.x. A, B, C or A1,A2,B2 etc) depending
	# on the data
	# The first column should always be 'NC_001526 genomic position'
	df = pd.read_excel(f,index_col=0)
	columns = ['NC_001526 genomic position', 'A', 'B', 'C','D',] # Should have the option to have other names
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
	# Prefix and empty should be provided from user, default Lin_, Other

def annotate_results(annot_f,probe_res,empty="Other"):
	#For the non ATGC character value will be automatically Other
	"""
	Probe results is the file that has mapped each nucleotide results of the probe blast
	"""
	annotdf = pd.read_excel(annot_f,index_col=0)
	probedf = pd.read_excel(probe_res,index_col=0)
	tmpl = annotdf.index.to_list()
	tmpl = ["NC_001526_" + str(l) for l in tmpl]
	annotdf.index = tmpl # So i can match the columns of probe results file
	annotdict = annotdf.to_dict(orient='index')
	
	posible_nucl_chars=list(np.unique(probedf.values)) # Non ATGC characters will be automatically value to empty
	for char in posible_nucl_chars:
		if char == "A" or char == "T" or char == "G" or char == "C":
			pass
		else:
			for probe in annotdict:
				annotdict[probe][char] = "Other"
	probedf = probedf.replace(annotdict)
	num_probes = len(probedf.columns)
	counters = {}
	indeces = probedf.index
	for idx in indeces:
		 counters[idx] = dict(Counter(probedf.loc[idx]))
	t = pd.DataFrame.from_dict(counters,orient='index')
	t.fillna(0,inplace=True)
	t = t.apply(lambda x:round((x/num_probes)*100,2))
	probedf = probedf.join(t)
	probedf.to_excel("ProbeF_lineage.xlsx") # TODO fix correctly the output file outputdir + query_f + "_ProbeF_lineage"


fin = "../input/signature_nucl_data_20210602_for_params.xlsx" # TODO: Fix to Take it from params
annot_f = "../input/annot.xlsx" # TODO: Fix to Take it from params
# create_annot_file(fin,annot_f)

probe_res = "../Script_out/10359_HPV16_genbank.fasta_Probes.xlsx"
# TODO: Fix to Take it from params
annotate_results(annot_f,probe_res)
