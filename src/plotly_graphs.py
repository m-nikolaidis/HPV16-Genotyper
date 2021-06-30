import re
import pathlib
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import AlignIO

def _sort_alphanumeric(iteratable):
    # Helper func Sort the given motif list alphanumerically :return: sorted list 
    int_convert = lambda text: int(text) if text.isdigit() else text 
    sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(iteratable, key = sorting_key)

def _add_reference_organisms(profile_dir,df):
    """In order to plot the different reference organisms, i need to add zero values to the first gene
    Use the first gene in the profile files
    Dataframe is the sequence specific dataframe, sorted by query start
    I am going to add at the start of the dataframe the 0 values for each reference, sorted alphanumerically
    """
    df_first_row = df.head(1)
    qstart  = df_first_row["qstart"].values[0]
    saccver = df_first_row["saccver"].values[0]
    _, gene = saccver.split("_") # Need to care for this. Either the input must have a special character or i need to find another way

    profile_f = pathlib.Path(gene + "_profile.fa") # Need to care for this, might create compatibility problems
    profile_f = profile_dir / profile_f
    alignment = AlignIO.read(profile_f,"fasta")
    headers = [alignment[x].id for x in range(len(alignment))]
    headers = _sort_alphanumeric(headers)
    for x in range(len(headers)):
        h = headers[x]
        m = re.match(r"HPV31.+",h) # Outgroup regex, should generalise this in a way
        if m:
            break
    del headers[x]

    tmp = []
    # headers.reverse() # In order to appear phylogenetically in HPV16 paradigm
    for org in headers:
        row = df_first_row.copy()
        row.loc[row.index,"saccver"] = org
        row.loc[row.index,"pident"] = 0
        tmp.append(row)
    tmpdf = pd.concat(tmp)
    final_df = pd.concat([tmpdf,df])
    return final_df

def GeneID_fig(blastresdf, seqname, output_dir, profile_dir):
    """
    """
    tmpdf = blastresdf[blastresdf["qaccver"] == seqname]
    columns = ["saccver","qstart","qend","pident"]
    tmpdf = tmpdf[columns]
    tmpdf = tmpdf.sort_values(by="qstart",axis=0)
    
    tmpdf = _add_reference_organisms(profile_dir, tmpdf)
    tmpdf.index = range(len(tmpdf))
    indeces = tmpdf.index
    for idx in indeces:
        saccver = tmpdf.loc[idx,"saccver"]
        db, gene = saccver.split("_") # Need to care for this. Either the input must have a special character or i need to find another way
        tmpdf.loc[idx,"Database"] = db
        tmpdf.loc[idx,"Gene"] = gene
        tmpdf.loc[idx,"Lineage"] = db[0] 
    tmpdf["pident"] = tmpdf["pident"].round(3) * 100
    discrete_colors = ["#00ff00","#0080ff","#ff8000","#ff0000"]
    fig = px.scatter(tmpdf, x="Gene", y="Database", size="pident", color="Lineage",
        hover_data=["qstart","qend","pident"], title=seqname, color_discrete_sequence = discrete_colors)
    html_outf = str(output_dir / pathlib.Path("Graphics/GeneSim") / pathlib.Path(seqname + ".html"))
    fig.write_html(html_outf,include_plotlyjs="cdn")



def SNP_fig(snp_resultsdf, seqname, output_dir, color_dict):
    """
    """
    tmpdf = snp_resultsdf[[seqname,"Color"]]
    indeces = tmpdf.index
    for i in indeces:
        m = re.match(r'^\S+_(\d+)$',i)
        lineage = tmpdf.loc[i,seqname]
        if lineage in color_dict:
            tmpdf.loc[i,seqname] = lineage
        else:
            tmpdf.loc[i,seqname] = "Other"
        tmpdf.loc[i,"SNP_pos"] = m.group(1)
    tmpdf["y"] = 1
    color_list = list(tmpdf["Color"])
    # Write in MD. This plot renames all the positions who are not specified in color_dict
    # as Other
    fig = px.scatter(tmpdf, x="SNP_pos", y="y", color=seqname,
        hover_data=[seqname], title=seqname, color_discrete_map=color_dict)
    html_outf = str(output_dir / pathlib.Path("Graphics/SNP_detection") / pathlib.Path(seqname + ".html"))
    fig.write_html(html_outf,include_plotlyjs="cdn")
