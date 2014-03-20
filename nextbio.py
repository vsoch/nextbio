# nextbio.py  - Manipulation of nextbio files for Leticia

# This script includes a set of functions for preprocessing of nextbio
# downloaded files, for download.  To use, you should import it into
# a python script, eg:
#
# import nextbio as nb
#
# Extract genes from a file
# nb.extractGenes("/path/to/myfile.csv")
# 
# Do pairwise comparisons for all gene lists in a folder
# nb.pairwiseCompare('/path/to/folder')

# See ALGORITHM.txt for the original algorithm (from Leticia) 
# see run_nextbio.py for more details about running and output
# vsochat@stanford.edu
# Wall-Lab, March 2014

import numpy as nu
from copy import deepcopy
import os
from os import listdir
from os.path import isfile, join
import datetime

# -- PREPROCESSING FUNCTIONS --------------------------------------------------------------

# Extract just UNIQUE gene names from a .csv or .txt file
# The extension MUST be .csv or .txt, and the extraction
# is done differently depending on the type
# Files are output to the same directory as the input
def extractGenes(filey):
    print "Extracting genes from " + filey
    myfile = open(filey,'r')        
    raw = myfile.readlines()
    symbols = []
    for f in range(5, len(raw)):
      symbols.append(raw[f].split(",")[0].strip("\"")) 
    symbols = nu.unique(symbols)
    outfile, extension = os.path.splitext(filey)
    # Print to file
    outfile = open(outfile + "_genes.txt",'w')
    for s in symbols:
      outfile.writelines(s + "\n")
    outfile.close()

# -- ANALYSIS FUNCTIONS --------------------------------------------------------------

# Pairwise compare does pairwise comparisons of all gene lists in a folder 
# of interest
def pairwiseCompare(folder):
    # Get files in the folder
    if folder[-1] not in ["/","\\"]:
      folder = folder + "/"

    # Read in route 66 file
    r66 = []
    tmp = open('resource/Route66Genes.txt','r').readlines()
    for r in range(0,len(tmp)):
      r66.append(tmp[r].strip("\n").strip("\r"))

    files = [ f for f in listdir(folder) if isfile(join(folder,f)) ]
    # Read in each file to a dictionary
    genes = dict()
    for f in files:
      print "Reading genes from " + f
      filey = open(folder + f,'r').readlines()
      tmp = []
      for line in filey:
        tmp.append(line.strip("\n"))
      genes[f] = tmp  

    # This matrix will hold the raw intersection counts
    matrix = nu.zeros(shape=(len(files),len(files)))
    # This matrix will hold route66 overlap counts
    matrix66 = nu.zeros(shape=(len(files),len(files)))
    # This dictionary will hold pairwise overlaps with route 66 for each pair
    route66 = dict()
    # We will use the analysis time for the output file
    analysisTime = datetime.datetime.now().strftime('%b-%d-%Y-%I-%M-%S')

    # Iterate through pairs, this could be done more efficiently :)
    print "Calculating overlap between genes and route66..."
    for f1 in range(0,len(files)):
      for f2 in range(0,len(files)):
        pairid = files[f1] + "||" + files[f2]
        # Save overlap counts
        matrix[f1,f2] = len(set.intersection(set(genes[files[f1]]),set(genes[files[f2]])))
        geneoverlap = set.intersection(set(genes[files[f1]]),set(genes[files[f2]]))
        route66[pairid] =  set.intersection(set(r66),geneoverlap)
        matrix66[f1,f2] = len(set.intersection(set(r66),geneoverlap))

    # Now we have a matrix of raw overlap between genes, and route66 overlap   
    # Print to file a list of route 66 overlap,
    print "Printing output files to present working directory..."
    r66overlap = open("r66_overlap_" + analysisTime + ".txt",'w')
    r66count = open("r66_counts_" + analysisTime + ".txt",'w')
    r66overlap.writelines("file1,file2,r66overlap_genes\n")
    r66count.writelines("file1,file2,r66overlap_count\n")
    for key,val in route66.iteritems():  
      pairs = key.split('||')
      r66overlap.writelines(pairs[0] + "," + pairs[1] + "," + str(list(val)) + "\n")
      r66count.writelines(pairs[0] + "," + pairs[1] + "," + str(len(val)) + "\n")
    r66overlap.close()
    labels = open("labels_" + analysisTime + ".txt","w")       
    for f in files:
       labels.writelines(f + "\n")
    labels.close() 
    nu.savetxt("r66_matrix_" + analysisTime + ".txt", matrix66, delimiter=",",fmt="%g")
    nu.savetxt("pairwiseGeneCount_matrix_" + analysisTime + ".txt", matrix, delimiter=",",fmt="%g")
  
if __name__ == "__main__":
  print "Please import as a module"
