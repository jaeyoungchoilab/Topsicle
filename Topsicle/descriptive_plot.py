# Descritive plot to have an overview of the data
# input: fasta, fa.gz, fastq.gz, or fastq files
# step 1: calculate Telomere-like Repeat Count (TRC) to identify read start / end has telomere 
# Step 2: telomere length identified through mean values change of all telomere-like through window sliding 

import sys
import os
# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import gzip
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO 
import logging

import re
import seaborn as sns
import matplotlib.pyplot as plt

colors = sns.color_palette("colorblind",n_colors=30)

sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

version_number = "1.0.0"
BoundTeloNano_output_prefix = "BoundTeloNano"

# check file type before processing 
def check_file_type(filepath):
    open_func = gzip.open if filepath.endswith(".gz") else open
    try:
        with open_func(filepath, "rt", encoding='utf-8') as handle:
            first_line = handle.readline().strip()
            if first_line.startswith("@"):
                return "fastq"
            elif first_line.startswith(">"):
                return "fasta"
            else:
                logging.warning('Format cannot be identified. Check the input.')
                return 0
    except Exception as e:
        logging.error(f"Error checking file type: {e}")
        return 0


def unzip_file(filepath):
    """Read sequences from a FASTQ or FASTA file (compressed or uncompressed)."""
    if not isinstance(filepath, str):
        logging.error("Input must be a string representing the file path.")
        return None

    file_type = check_file_type(filepath)
    if not file_type:
        logging.error("File type could not be determined or is unsupported.")
        return None

    open_func = gzip.open if filepath.endswith(".gz") else open
    try:
        with open_func(filepath, "rt", encoding='utf-8') as handle:
            for sequence in SeqIO.parse(handle, file_type):
                yield sequence
    except Exception as e:
        logging.error(f"Error parsing file: {e}")

        #sys.exit("there is problem with input in", filepath, " Terminate.")

# descriptive plot
def descriptive_plot(filepath, pattern, minSeqLength):
    '''
    Plot to view the distribution of telomere-like pattern in the sequence
    filepath: the path of input file
    pattern: pattern to plot - will plot the input pattern here
    output: a plot showing where are telomere-like pattern distribute
    '''
    base_name = filepath.split('/')[-1]
    file_name = base_name.split('.')[0]
    fig, ax = plt.subplots(figsize=(10, 15))
    sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})
    markers = ['|','|']
    trans_table = str.maketrans('ACGT', 'TGCA')
    
    # Apply the translation to the sequence
    patterns = [pattern.upper(),pattern.translate(trans_table).upper()]
    # For viz purpose 
    labels_pre = [pattern.upper(),pattern.translate(trans_table).upper()[::-1]]
    k_line = 0
    read_ids=[]
    added_labels=set()
    iin=0
    if unzip_file(filepath) == None: 
        print('problem in filepath, can not have descriptive plot')
        return None 
    for read in unzip_file(filepath):
        if len(read.seq) > minSeqLength:   # length of sequence > minSeqLength, which default - 9kbp 
            iin+=1
            
            seq = str(read.seq[:minSeqLength]).upper()   # no need to plot all, just to 9kbp is enough to observe
            read_ids.append(read.name)
            seq_2 = str(read.seq[::-1][:minSeqLength]).upper() # comp strand (reverse of reverse comp)
            
            # Loop through each pattern
            for i, pattern in enumerate(patterns):
                matches = [m.start() for m in re.finditer(re.compile(pattern), seq)]
                if pattern not in added_labels:
                    ax.scatter(matches, [k_line] * len(matches), color=colors[i], marker='|', label=pattern, zorder=2)
                    added_labels.add(pattern)  # Mark this pattern as added
                else:
                    ax.scatter(matches, [k_line] * len(matches), color=colors[i], marker='|', zorder=2)
           
                matches_2 = [m.start() for m in re.finditer(re.compile(pattern), seq_2)]
                # tab in to indicate the pattern found at the end 
                ax.scatter([x for x in matches_2], [k_line] * len(matches_2), color=colors[i], marker='|', zorder=2)
            # Increment the line position for the next read
            k_line += 2
            
        else:
            pass
        # Add labels and title
        ax.set_title(f'Location of telomere patterns in {file_name}')
        ax.set_xlabel('Position')

        # Add legend to differentiate patterns, but make it reverse 
        handles, labels = ax.get_legend_handles_labels()
    
        # handles.reverse()
        # labels.reverse()
        # Set the legend with reversed order
        labels=labels_pre
        ax.legend(handles, labels, title="Pattern")  

        # Set y-axis ticks and labels
        ax.set_yticks([i*2 for i in range(len(read_ids))])
        ax.set_yticklabels([read_ids[i] for i in range(len(read_ids))])

        # get the grid 
        ax.xaxis.grid(True)  # Enable horizontal grid lines
        ax.yaxis.grid(True)  # Disable vertical grid lines
        plt.tight_layout()

        if iin > 40:
            print("file has more than 40 reads, but it is not reccomended to have plot with that many reads")
            print("so the output plot will have 40 reads only")
            # show plot
            return "plotted"
    # show plot
    return "plotted"

def pattern_scramble_telo (pattern, cut_length):
    '''
    pattern: string, input desired patterns to scramble
    cut_length: list, length wanna cut, recommend having 3 and 4 cut 
    '''
    # double the pattern since telomere sequences are repetitive sequences of pattern
    unique_cuts = set()
    pattern = pattern + pattern 
    pattern = pattern.upper()

    # transform to list for further work 
    if type(cut_length)!= list:
        cut_length = [cut_length]
       
    # Loop over each cut length (3 and 4)
    for length in cut_length:
        for i in range(len(pattern) - length + 1):
            # Extract the substring of the given length
            substring = pattern[i:i + length]
            # Add to set (to ensure uniqueness)
            unique_cuts.add(substring)
    
    # Convert set to sorted list to keep order
    return sorted(list(unique_cuts))

def patterns_to_search (telopattern, cut_length):
    '''
    telopattern: str or list of str, pattern of telomere
    cut_length: 4-bp from 6-bp pattern 
    '''
    if "|" in telopattern and type(telopattern) == str: 
        ## origin pattern: 'AACC|ACCG',...
        origin_pattern = telopattern
        origin_pattern_rev = origin_pattern[::-1]
        
        ## comp_pattern, the "TTTGGC" one 
        ### make trans_table 
        trans_table = str.maketrans('ACGT', 'TGCA')
    
        # Apply the translation to the sequence
        comp_pattern = origin_pattern.translate(trans_table) 
        comp_pattern_rev = comp_pattern[::-1] 
        pattern_all = origin_pattern.upper() + comp_pattern.upper()
        
    if "|" not in telopattern and type(telopattern) == str: 
        # pattern scramble and work   
        ## origin pattern directly from the telopattern, "AAACCG"
        origin_pattern = pattern_scramble_telo(telopattern, cut_length=[cut_length]) 
        origin_pattern_rev = [seq[::-1] for seq in origin_pattern]
        
        ## comp_pattern, the "TTTGGC" one 
        ### make trans_table 
        trans_table = str.maketrans('ACGT', 'TGCA')
    
        # Apply the translation to the sequence
        comp_pattern = [seq.translate(trans_table) for seq in origin_pattern]
        comp_pattern_rev = [seq[::-1] for seq in comp_pattern]

        # all patterns 
        pattern_all = origin_pattern + comp_pattern
        pattern_all=[element.upper() for element in pattern_all]

    if type(telopattern) == list:
        pattern_all=[element.upper() for element in telopattern]

    return pattern_all

def patterns_vs_match_heatmap (filepath, telopattern,telophrase,minSeqLength):
    '''
    fastqz_loc: string, full path location of sequence, but just 1 path at a time 
    '''
    # Split by '/' to get the last component
    base_name = filepath.split('/')[-1]
    file_name = base_name.split('.')[0]
    sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})
    trans_table = str.maketrans('ACGT', 'TGCA')
    
    # Apply the translation to the sequence
    pattern_all = pattern_scramble_telo (telopattern, cut_length=telophrase)
    # double check patterns we look for
    # print("Patterns topsicle is finding:")
    # print(pattern_all)

    matches=[]
    matches_2_list=[]
    if unzip_file(filepath) == None: 
        print('problem in filepath, can not have heatmap')
        return None 
    
    for read in unzip_file(filepath):
        if len(read.seq) > minSeqLength:
            read_ids=[]
            seq = str(read.seq[:minSeqLength]).upper()   # no need to plot all, just to 9kbp is enough to observe
            read_ids.append(read.name)
            seq_2 = str(read.seq[::-1][:minSeqLength]).upper() # reverse comp)
            trans_table = str.maketrans('ACGT', 'TGCA')
    
            # Apply the translation to the sequence
            seq_2 = seq_2.translate(trans_table)

            # Loop through each pattern
            for i, pattern in enumerate(pattern_all):
                finding=int(len(telopattern)-telophrase)
                dots='.'*finding  # number of match after pattern 
                regex = re.compile(rf"{pattern}({dots})")
                # Find all matches
                # print("regex",regex)
                matches_1 = regex.finditer(seq)
                matches_2 = regex.finditer(seq_2)
                
                y_axis_order=set()
                # Append the pattern and matches to the results list
                for match in matches_1:
                    y_axis_order.add(match.group(1))
                    matches.append((pattern, match.group(1),read_ids))    # get pattern and the match after it
                for match in matches_2:
                    y_axis_order.add(match.group(1))
                    matches_2_list.append((pattern, match.group(1),read_ids))    # get pattern and the match after it 
    
    matches_df = pd.DataFrame(matches, columns=["Pattern", "Match","read id"]) 
    matches_2_df = pd.DataFrame(matches_2_list, columns=["Pattern", "Match","read id"]) 
    y_axis_order=list(y_axis_order)
    matches_df['Match'] = pd.Categorical(matches_df['Match'],y_axis_order)
    matches_2_df['Match'] = pd.Categorical(matches_2_df['Match'],y_axis_order)

    fig, ax = plt.subplots(1, 2, figsize=(10, 8))  # 1 row, 2 columns
    y_axis_order=list(y_axis_order)
    # print('matches_df',matches_df)

    # matches_2_df = matches_2_df.loc[y_axis_order]
    sns.histplot(data=matches_df, ax=ax[0],
                x="Pattern", y="Match",cbar=True, cbar_kws=dict(shrink=.75)).set_title("forward strand")
    ax[0].tick_params(axis='x', rotation=45)
    sns.histplot(data=matches_2_df,ax=ax[1],
                x="Pattern", y="Match",cbar=True, cbar_kws=dict(shrink=.75)).set_title("reverse strand")
    ax[1].tick_params(axis='x', rotation=45)
    plt.suptitle(f"{telophrase}-bp patterns and matches from reads in \n {file_name}")
    plt.tight_layout()

    return matches_df
