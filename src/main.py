import os
import sys
import pathlib
import logging
import datetime
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
# Tool modules
import annot
import blast
import params
import _env_setup
import phylogeny
import plotly_graphs

def _init_binaries(system):
	"""Initialize the paths for the needed binaries depending on the system
	input: System info
	return: List of pathlib paths
	"""
	if "win" in system:
		logging.critical(" Windows BLAST not implemented yet ")
		raise NotImplementedError(" Windows BLAST not implemented yet ")
		muscle_bin= os.path.join(str(pathlib.Path(__file__)), 
			str(pathlib.Path('../resources/Win_bin/muscle3.8.31.exe'))).replace(os.path.basename(__file__) + "\\", "")
		raise NotImplementedError(" Tree creation on windows is not implemented yet")

	if "linux" in system:
		makeblastdb_bin = os.path.join(str(pathlib.Path(__file__)), 
			str(pathlib.Path('../resources/Linux_bin/makeblastdb'))).replace(os.path.basename(__file__) + "/", "")
		blastn_bin = os.path.join(str(pathlib.Path(__file__)), 
			str(pathlib.Path('../resources/Linux_bin/blastn'))).replace(os.path.basename(__file__) + "/", "")
		muscle_bin = os.path.join(str(pathlib.Path(__file__)), 
			str(pathlib.Path('../resources/Linux_bin/muscle3.8.31_i86linux'))).replace(os.path.basename(__file__) + "/", "")
		seaview_bin = os.path.join(str(pathlib.Path(__file__)), 
			str(pathlib.Path('../resources/Linux_bin/seaview4'))).replace(os.path.basename(__file__) + "/", "")
		phyml_bin = os.path.join(str(pathlib.Path(__file__)), 
			str(pathlib.Path('../resources/Linux_bin/phyml-mpi'))).replace(os.path.basename(__file__) + "/", "")
	return [makeblastdb_bin,blastn_bin,muscle_bin,seaview_bin,phyml_bin]

def _create_dirs(outdir):
		"""
		Perform the necessary actions to set up the script working environment
		"""
		script_out = outdir
		tmp_out = script_out / pathlib.Path("tmp_dir")
		tmp_batch_fasta_out = tmp_out / pathlib.Path("Batch_fasta_files") # Batch complete genomes
		tmp_batch_genes_out = tmp_out / pathlib.Path("Batch_genes") # Gene files of batch size
		tmp_alignments_out = tmp_out / pathlib.Path("Alignments") # Alignment of batch genes
		profiles_out = script_out / pathlib.Path("Profile_alns") # Final profile alns for each gene
		trees_out = script_out / pathlib.Path("Phylogenetic_Trees")
		graphics_out = script_out / pathlib.Path("Graphics")
		gene_sim_graphics = graphics_out / pathlib.Path("GeneSim")
		snp_graphics = graphics_out / pathlib.Path("SNP_detection")

		script_out.mkdir(exist_ok=False)
		tmp_out.mkdir(exist_ok=False)
		tmp_batch_fasta_out.mkdir(exist_ok=False)
		tmp_batch_genes_out.mkdir(exist_ok=False)
		tmp_alignments_out.mkdir(exist_ok=False)
		profiles_out.mkdir(exist_ok=False)
		trees_out.mkdir(exist_ok=False)
		graphics_out.mkdir(exist_ok=False)
		gene_sim_graphics.mkdir(exist_ok=False)
		snp_graphics.mkdir(exist_ok=False)

		# TODO: Warn the user before trying to initialize. Throw message in case the app fails

# Create the needed arguments
parser = argparse.ArgumentParser(description="HPV16 genotyping tool CLI")
parser.add_argument('--p_in', metavar="F",type=str,help='Absolute or relative path to the existing parameters file')
parser.add_argument('--p_out', metavar="F",type=str,help='Directory to write the default parameters file')
parser.add_argument('--clean', metavar="bool",type=bool,help='Remove the tmp directories created during execution (True / False)')
# TODO: Add a dual switch to these. Either read or write parameters file, and make them positional
args = parser.parse_args()

params_o = args.p_out
params_i = args.p_in
params_clean = args.clean

now = datetime.datetime.now()
start_dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

if params_o:
	params.default_param_file(params_o)
	sys.exit(0) # Kill the CLI when outputing parameters file
create_annot = False
if create_annot:
	fin = "../input/signature_nucl_data_20210602_for_params.xlsx" # TODO: Fix to Take it from params
	annot_f = pathlib.Path("../input/annot.xlsx") # TODO: Fix to Take it from params
	annot.create_annot_file(fin,annot_f)
	# Output this in the input directory?
paramsdf = pd.read_excel(params_i,index_col=0)
outdir = pathlib.Path(paramsdf.loc["out","Value"])
# _create_dirs(outdir)
logfile = outdir / pathlib.Path("logfile.log")
logging.basicConfig(filename=logfile, format='%(asctime)s %(message)s', level=logging.DEBUG)
logging.info(" Loaded parameters file successfully ")
logging.info(" Environment set successfully in dir %s" %(outdir))
system = paramsdf.loc["system","Value"]
indir = pathlib.Path(paramsdf.loc["in","Value"])
query_f = pathlib.Path(paramsdf.loc["query","Value"])
query_f_path = indir / query_f
_env_setup.file_exists(query_f_path)
_env_setup.input_file_format(query_f_path)
query_f_filt = _env_setup.filter_input_file(query_f_path)
query_f = query_f_filt # Posix | Windows Path
query_f_path = indir / query_f
_env_setup.dir_priviledges(outdir)
makeblastdb_bin,blastn_bin,muscle_bin,seaview_bin,phyml_bin = _init_binaries(system)
workflow = "gene identification"
geneid_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin)
geneid_blastn_res_path_xlsx = blast.parse_GeneID_results(geneid_blastn_res_path,paramsdf)
workflow = "snp"
snp_blastn_res_path = blast.blastn_search(paramsdf, query_f_path, workflow, makeblastdb_bin, blastn_bin)
snp_blastn_res_path_xlsx = blast.parse_SNP_results(snp_blastn_res_path)
aln_files = phylogeny.prepare_batch_alns(outdir, query_f_path,geneid_blastn_res_path_xlsx, muscle_bin ,batch_size=1)
profiledb_dir = paramsdf.loc["GenesProfile_directory","Value"]
aln_files = phylogeny.profile_aln(outdir,aln_files,profiledb_dir, muscle_bin)
trees_dir = phylogeny.build_trees(outdir, aln_files, phyml_bin, seaview_bin)

visualize_trees = False
if visualize_trees:
	phylogeny.visualize_trees() # TODO: Needs implementation

geneid_res_df = pd.read_excel(geneid_blastn_res_path_xlsx)
organisms = list(np.unique(geneid_res_df["qaccver"]))
for seqname in organisms:
	plotly_graphs.GeneID_fig(geneid_res_df,seqname, outdir, profiledb_dir)

annot_f = pathlib.Path("../input/annot.xlsx") # TODO: Need to make this automatic
num_snps = annot.annotate_results(annot_f,snp_blastn_res_path_xlsx)

# if params_clean:
# 	_env_setup.clean_tmp()
# 	# TODO: Implement this
snp_results = outdir / "SNPs_results.xlsx"
snp_resultsdf = pd.read_excel(snp_results,index_col=0)
snp_resultsdf = snp_resultsdf.T
snp_resultsdf = snp_resultsdf.head(num_snps)
snp_resultsdf["Color"] = ""
color_dict = {
    "Lin_A":"#00ff00",
    "Lin_B":"#0080ff",
    "Lin_C":"#ff8000",
    "Lin_D":"#ff0000",
    "Lin_BCD":"#ebeb34",
    "Other":"#ffffff"
}
for seqname in organisms:
	plotly_graphs.SNP_fig(snp_resultsdf, seqname, outdir, color_dict)
# TODO: Implement function to choose colors?. On gui for sure

end_dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print("Application started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))
# logging.info("The scripts started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))
