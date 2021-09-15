# The tool
The tool in CLI has been tested to run in Linux, WSL and Windows (W10).
CHeck in new install of windows if the binaries run correctly
----------------------------------------------------------------------
# Use DataFrame.to_markdown for the parameters excel

# Write in MD, that the user should have the required priviledges in the chosen output directory

From params.py
""" To write in .md
	systems supported at the moment: Linux, WSL
	The script can be installed in whatever folder that has user privileges
	input directory for files is script/../input by default
	Evalue notation accepts both scientific and floating point (due to pandas)
	Word_size should be int >= 4
	num_threads is system threads - 2 by default
	Currently using blast+ 2.11.0 #TODO: Update, 2.12
	All the files should be put in the desired input directory, except for the params excel. which can be wherever
	The database and query should be only the file names. Not the paths
	Do not delete the first line of the excel file that has the word Value. It is needed by the script
	TODO: the annotation should be given only as the file name, not the whole directory
"""


## Blastn
"""Write in md:
oufmt columns can be altered based on user needs directly from the excel params file.
The addition of the text should be in the same cell and separated by spaces
The raw blast_results are automatically turned into an excel file for ease of interpretation (The txt file is also kept)
The script will read either the csv or the excel file, whichever is available. If none is, it will throw an error
TODO: Add list of all possibilities for non-cli users

Add something about the workflow on the MD file?
"""

### Recombinant sequences
Recombinant sequences are considered sequences with conflicting lineaeges in the gene identification results

# Phylogeny
## Profile alignment
We are using muscle profile alignment to save computational time and storage space.
Muscle and seaview binaries are taken from wsl compiles. Work on WSL. 
Check versions to add to MD

## Trees
Trees are built with PhyML by default

# SNP Annotation
## Create SNP annotation file
Input excel must have the specific format
Ref genomic position	A	B	C	D
The first column should always be the reference genomic position (e.x. NC_001526 genomic position)
NC_001526 genomic position	A	B	C	D
286	T	A	A	A
289	A	G	G	G
335	C	T	T	T

make the above thing a table

create_annot_file function:
Prefix and empty should be provided from user, default Lin_, Other

# Graphical outputs
SNP graphical outputs are a bit buggy in the windows execution
TODO: Need to investigate
## Gene best hits
 requires an active internet connection to view the plotly graphs

## SNPs
This plot renames all the SNPs (Non A,B,C,D, BCD) who are not specified in color_dict as Other
requires an active internet connection to view the plotly graphs

## Tree visualizations (interactive)

## Tree rendering

# The time consuming functions have on-off switch
Blastn search and phylogentic tree construction
Going to use them as checkpoints



# Used parameters
Check if the provided gene sequence has many Ns and filter the parent organism from the analysis
Default N_perc=50%

If more than 20 organisms are in the filtered input, then the BioNJ method will be used

Also filtering the non HPV-16 genomes


# output 
If no output directory is specified the results are thrown into the script directory with the Fasta input

filename prefix (without the extension) and the current date

# lineage snp graph
This graph is showing the A,B,C,D specfic SNPs and the BCD specific
So we can divide the A from the non-A

# Recombinant identification
Use the cut-offs to identify putative recombinants / artifacts in the analysis
Gene identificatio cut-off: Atleast 1 gene with different main lineage
Lineage specific SNPs: 3 or more SNPs other than the dominant lineage