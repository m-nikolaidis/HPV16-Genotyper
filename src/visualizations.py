import re
import pathlib
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import AlignIO
from ete3_custom import Tree, TreeStyle, NodeStyle, faces, AttrFace, TreeFace, TextFace

def GeneID_fig(blastresdf, output_dir, profile_dir, color_dict, recombinants):
	"""
	"""
	columns = ["saccver","qstart","qend","pident"]
	blastresdf = blastresdf.sort_values(by="qstart",axis=0)
	blastresdf["Sequences"] = blastresdf["qaccver"]
	indeces = blastresdf.index
	for idx in indeces:
		saccver = blastresdf.loc[idx,"saccver"]
		db, gene = saccver.split("_") # Need to care for this. Either the input must have a special character or i need to find another way
		blastresdf.loc[idx,"Database"] = db
		blastresdf.loc[idx,"Gene"] = gene
		blastresdf.loc[idx,"Lineage"] = db[0] 
	blastresdf["pident"] = blastresdf["pident"].round(3) * 100
	fig = px.scatter(blastresdf, x="Gene", y="Sequences", size="pident", color="Lineage",text="Database",
		hover_data=["Database","qstart","qend","pident"], title="Gene identification all sequences", color_discrete_map = color_dict)
	html_outf = str(output_dir / pathlib.Path("Graphics/GeneSim") / pathlib.Path("All_sequences.html"))
	fig.write_html(html_outf,include_plotlyjs="cdn")
	# Recombinants only figure
	rec_only_df = blastresdf[blastresdf["qaccver"].isin(recombinants)]
	recfig = px.scatter(rec_only_df, x="Gene", y="Sequences", size="pident", color="Lineage",text="Database",
			hover_data=["Database","qstart","qend","pident"], title="Gene identification recombinant sequences", color_discrete_map = color_dict)
	html_outf = str(output_dir / pathlib.Path("Graphics/GeneSim") / pathlib.Path("Recombinant_sequences.html"))
	recfig.write_html(html_outf,include_plotlyjs="cdn")
	rec_only_df.to_excel(output_dir/"RecOnlyGeneID_to_usedf.xlsx")

def SNP_fig(snp_resultsdf_annot, output_dir, color_dict, snp_resultsdf, recombinants):
	"""
	"""
	c_dict_k = list(color_dict.keys())
	def _rename(x):
		if x in c_dict_k:
			return x
		else:
			return "Other"
	def _get_pos(x):
		m = re.match(r"^\S+_(\d+)$",x)
		return m.group(1)
	snp_resultsdf["SNP_pos"] = snp_resultsdf.index
	snp_resultsdf["SNP_pos"] = snp_resultsdf["SNP_pos"].apply(lambda x: _get_pos(x))
	snp_resultsdf.index = snp_resultsdf["SNP_pos"]
	snp_resultsdf_annot_tr = snp_resultsdf_annot.applymap(lambda x:_rename(x))
	snp_resultsdf_annot_tr["SNP_pos"] = snp_resultsdf_annot_tr.index
	snp_resultsdf_annot_tr["SNP_pos"] = snp_resultsdf_annot_tr["SNP_pos"].apply(lambda x: _get_pos(x))
	snp_resultsdf_annot_tr.index = snp_resultsdf_annot_tr["SNP_pos"]
	snp_resultsdf_annot_tr = snp_resultsdf_annot_tr.drop("SNP_pos",axis=1)
	rows, col = snp_resultsdf_annot_tr.shape
	indeces = range(rows*col)
	cols = ["Sequences","SNP reference pos","Nucleotide","Lineage"]
	df = pd.DataFrame(index=indeces,columns=cols)
	snp_pos_l = snp_resultsdf_annot_tr.index
	cols = list(snp_resultsdf_annot_tr.columns)
	count = 0
	for snp_pos in snp_pos_l:
		for seq in cols:
			lineage = snp_resultsdf_annot_tr.loc[snp_pos,seq]
			nucl = snp_resultsdf.loc[snp_pos,seq]
			df.loc[count,"Sequences"] = seq
			df.loc[count,"SNP reference pos"] = snp_pos 
			df.loc[count,"Lineage"] = lineage
			df.loc[count,"Nucleotide"] = nucl
			count += 1
	fig = px.scatter(df, x="SNP reference pos", y="Sequences", color="Lineage",
		hover_data=["SNP reference pos","Lineage","Nucleotide"], title="SNP identification all sequences", color_discrete_map=color_dict)
	fig.update_layout(xaxis_type = 'linear')
	html_outf = str(output_dir / pathlib.Path("Graphics/SNP_detection") / pathlib.Path("All_sequences.html"))
	fig.write_html(html_outf,include_plotlyjs="cdn")
	
	# Recombinants only figure
	recombinants.extend(["Nucleotides"])
	rec_only_df = df[df["Sequences"].isin(recombinants)]
	recfig = px.scatter(rec_only_df, x="SNP reference pos", y="Sequences", color="Lineage",
		hover_data=["SNP reference pos","Lineage","Nucleotide"], title="SNP identification recombinant sequences", color_discrete_map=color_dict
	)
	recfig.update_layout(xaxis_type = 'linear')
	html_outf = str(output_dir / pathlib.Path("Graphics/SNP_detection") / pathlib.Path("Recombinant_sequences.html"))
	recfig.write_html(html_outf,include_plotlyjs="cdn")
	rec_only_df.to_excel(output_dir / "RecOnlySNP_to_usedf.xlsx")

def render_trees(trees_dir, outdir, mode="jpeg"):
	"""
	Type of file to render the phylogenetic trees (Default jpg) \n available types: pdf, jpeg, png
	"""
	gene_name_regex = re.compile(r"^\S+_(\S+)_aln.+$")
	
	# Load trees
	trees = trees_dir.glob("**/*")
	
	mode = "." + mode
	trees = [f for f in trees if f.is_file()]

	for tree_f in trees:
		tree_f_str = str(tree_f)
		t = Tree(tree_f_str)
		# Beautify tree
		t.ladderize(direction=1)
		ts = TreeStyle()
		m = re.match(gene_name_regex, tree_f.name)
		gene = m.group(1)
		ts.title.add_face(TextFace(gene),column=1)
		ts.show_branch_support = True
		# Outputs
		tree_fout = tree_f.with_suffix(mode)
		tree_fout = str(outdir / pathlib.Path("Graphics/Trees_images") / tree_fout.name)
		t.render(tree_fout, tree_style=ts)
