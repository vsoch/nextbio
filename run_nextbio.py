#!/usr/bin/python

# STEP 0 ---- import the module
import nextbio as nb


# STEP 1 ---- extract genes

# Specify our input files (can also be done by script)
files = ["/home/vanessa/Documents/Dropbox/Code/Python/nextbio/GSE465_b1.txt","/home/vanessa/Documents/Dropbox/Code/Python/nextbio/GSE574_b1.csv","/home/vanessa/Documents/Dropbox/Code/Python/nextbio/GSE1007_b1.csv"]

# Create "gene lists" for each input file - will be put into the present working directory
# This assumes the common format downloaded by nextbio.  if there are errors in formatting,
# you will get an error you need to troubleshoot.
for f in files:
  nb.extractGenes(f)

# STEP 2 --- PAIRWISE COMPARISONS
# When you finish the above, you will have files appended with *_genes.txt, each a list of genes
# specified in the file.  Now you should put the ones you want to pairwise compare in a folder,
# and run the following:
nb.pairwiseCompare("/home/vanessa/Documents/Dropbox/Code/Python/nextbio/lists")

# When you finish with the above, you will have the following output in the PWD:

# r66_overlap_[date].txt: A three column file with file1,file2, and gene overlap with route 66:
# This means we take the two file overlap, and then overlap that set with route 66 genes
file1,file2,r66overlap_genes
GSE574_b1_genes.txt,GSE574_b1_genes.txt,['ITGB1', 'RBMX', 'AGTPBP1', 'UBE2D3', 'ASAH1']


# r66_counts_[date}.txt: A three column file with file1,file2, overlap with route 66
file1,file2,r66overlap_count
GSE574_b1_genes.txt,GSE574_b1_genes.txt,5

# r66_matrix_[date].txt: A file with equivalent info, but in a square matrix format.  Labels are in:
# labels_[date].txt, a single column that describes the row/column label
GSE1007_b1_genes.txt
GSE465_b1_genes.txt
GSE574_b1_genes.txt

# pairwiseGeneCount_matrix_[date].txt:  This has the RAW overlap counts for pairwise comparisons BEFORE
# route 66 is involved, with same corresponding labels.



