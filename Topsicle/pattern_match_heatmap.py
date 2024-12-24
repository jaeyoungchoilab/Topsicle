import sys
import os
# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import gzip
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO 

import ruptures as rpt
import re
import seaborn as sns
import matplotlib.pyplot as plt

# setting up 
colors = sns.color_palette("colorblind",n_colors=30)
sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

# need to find chr_all -- telo length of all reads from all pattern in that chromosome 

def read_pattern_matches(fastqz_loc, patterns, chr_all):
    '''
    fastqz_loc: string, full path location of sequence, but just 1 path at a time 
    chr_all: telo length of all reads from all pattern in that chromosome 
    return: a dataframe with pattern, match, type of base pair, length of telomere, read id 
    '''

    # patterns we are looking at 
    patterns = patterns  # List of patterns to search for

    # main 
    # Loop through reads in file, find interested read. Then find their matches and plot 
    match_heatmap = []
    with gzip.open(fastqz_loc, "rt") as handle:
        for seq in SeqIO.parse(handle, "fastq"):  # loop over each row in the file 
            #print(seq.id)
           
            for i in range(len(patterns)): 
                pattern = patterns[i]  # get pattern from the list of pattern 
                seq_r = str(seq.seq)
                regex = re.compile(f"{pattern}(.{{3}})")
                
                # Find all matches
                matches = regex.finditer(seq_r)
                # Append the pattern and matches to the results list
                for match in matches:
                    match_heatmap.append((pattern, match.group(1),"3bp",seq.id))    # get pattern and the match after it 
                    
            
            # pattern with 4 bp 
            for i in range(len(patterns_2)): 
                pattern = patterns_2[i]
                        
                seq_r = str(seq.seq)
                #regex = re.compile(f"{pattern}(.{{2}})")
                        # Find all matches
                #matches = re.finditer(regex,seq_r)

                matches = re.finditer(f"(?<={pattern}).{{2}}", seq_r)
                
                    # Append the pattern and matches to the results list
                for match in matches:
                    match_heatmap.append((pattern, match.group(),"4bp",seq.id))    # get pattern and the match after it 
                     
                    #print('next from 4bp')
    
    matches_all = pd.DataFrame(match_heatmap, columns=["Pattern", "Match","Type","read id"])

    return matches_all
