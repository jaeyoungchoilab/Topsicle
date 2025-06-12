# Topsicle
# input: fasta, fa.gz, fastq.gz, or fastq files
# step 1: calculate Telomere-like Repeat Count (TRC) to identify read start / end has telomere 
# Step 2: telomere length identified through mean values change of all telomere-like through window sliding 
# Also have visualizations (density plot of mean window change)

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import gzip
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO 
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# algo packages 
import ruptures as rpt
import re

# logging 
import logging
logging.basicConfig(level=logging.ERROR)
import warnings

# plotting
import seaborn as sns
import matplotlib.pyplot as plt
colors = sns.color_palette("colorblind",n_colors=30)
sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

version_number = "1.0.0"

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

# step 1: will that read have telomere?

# step 0: pre processing inputs
## cut out the telomere pattern from 6-bp to cut_length basepair
## such as 3bp or 4bp patterns
def pattern_scramble_telo (pattern, cut_length):
    '''
    pattern_scramble_telo: get the unique cuts from the initial sequence, along with pattern_to_search
    pattern: string, input desired patterns to scramble
    cut_length: list, length wanna cut, recommend having 3 and 4 cut 
    output: list of unique k-mers
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
    patterns_to_search: get the unique k-mers from the initial sequence, along with pattern_scramble_telo
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

def unzip_file(filepath):
    """
    Read sequences from a FASTQ or FASTA file (compressed or uncompressed)
    input: the path of the file
    yield sequence content 
    """
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

# step 1: count TRC 
def patternTRC_count(filepath, telopattern, read_length=0, kmer=4, no_bp=1000, cutoff=0.5):
    '''
    Find the telomere-like repeat count (TRC) in the first and last no_bp base pairs of the sequence
    filepath: str, location of file, can be either fastqz.gz or fasta, just 1 file at the time 
    no_bp: int, number of first bp to find 
    telopattern: str or list of str, pattern of telomere
    read_length: minimum length require for the sequence (filter out too short sequences)
    kmer: length of scanning pattern divided by chunk when we wanna have flexibility, default = 4
    cutoff: cutoff for mean_window value to be cutoff point for telomere boundary, default = 0.5
    yeild: list of reads with patterns and TRC value. Won't filter out reads at this step
    '''
    find_match = []  # empty list to store found pattern, ratio, start/end 
    raw_count = []

    # get the list of patterns to search for 
    pattern_all = patterns_to_search(telopattern, cut_length=kmer)
    compiled_patterns = [re.compile(pattern) for pattern in pattern_all]

    if isinstance(filepath, list):
        print("Can only process 1 file path at the time, please loop paths through the list")
        return None 

    for seq in unzip_file(filepath): 
        if len(seq.seq) > read_length:  # have to be longer than the minimum length required 
            seq_start = str(seq.seq[:no_bp].upper())  # get first no_bp base pairs 
            seq_end = str(seq.seq[-no_bp:][::-1].upper())  # get last no_bp base pairs   
            ratio_perfect_hit = no_bp / len(telopattern)
            itcs_indi = [['seq id', 'pattern', 'TRC forward strand', 'TRC reverse strand']]
            
            for pattern in compiled_patterns: 
                matches_start = len([m.start() for m in pattern.finditer(seq_start)])
                matches_end = len([m.start() for m in pattern.finditer(seq_end)])
                
                count_pattern_start = matches_start / ratio_perfect_hit
                count_pattern_end = matches_end / ratio_perfect_hit
                itcs_indi.append([seq.id, pattern.pattern, count_pattern_start, count_pattern_end])

            raw_count.append(itcs_indi)
            max_row_start = max(itcs_indi[1:], key=lambda row: row[2])
            max_row_end = max(itcs_indi[1:], key=lambda row: row[3])
            
            if max_row_start[2] > max_row_end[3]:
                if max_row_start[2] > cutoff:    
                    find_match.append([max_row_start[0], max_row_start[1], 'forward', max_row_start[2]]) 
            else:
                if max_row_end[3] > cutoff: 
                    find_match.append([max_row_end[0], max_row_end[1], 'reverse', max_row_end[3]])
                    
    if check_file_type(filepath) is None:
        print("can not read in file - can not run step 1")
        return None

    return find_match

# step 2: having telomere - how long is that telomere? 
def seq_cut_windows(s, window_size, step):
    """
    Generate overlapping windows of a given size from a sequence along with their start and end indices.
    
    seq: The input sequence (can be a list, string, etc.)
    window_size: The size of each window
    step: The number of elements to slide between windows (default is 1)

    return: A tuple containing the window start index, end index, and the window itself.
    """
    # Loop to create windows with overlap
    windows=[]
    for i in range(0, len(s) - window_size + 1, step):
        start = i  # Starting index of the window
        end = i + window_size - 1  # Ending index of the window
        if end > len(s):
            end = len(s)
        windows.append((start, s[start:end]))
    return windows  # Return indices and window
    
def bound_detect(filepath, read, pattern_telo, windowSize, slide, trimfirst, maxlengthtelo, cut_length, tail=None, plot_yes_no=None,plotcp_range=None):
    '''
    Find telomere-subtelomere boundary point in the sequence by using mean window value change and changepoint algo 
    filepath: location of fasta file
    read: read name, but just 1 read at the time 
    pattern: pattern wanna find, maybe ACCG. multiple patterns at each time
    windowSize: size of window
    slide: step of each window
    region: cut the sequence - in telomere, just check for first 5000 bp
    trimfirst: number of nucleotide trim off 
    maxlengthtelo: max possible length of telomere in this species

    return: telomere boundary point of each window 
    '''

    # Loop through reads in file, find interested read. Then find their matches and plot 
    read_ids = []
    boundary = []
    if not isinstance(read, str):
        print("can only read in 1 read at a time")
        return None

    patterns = patterns_to_search(pattern_telo, cut_length=cut_length)
    compiled_patterns = [re.compile(pattern) for pattern in patterns]

    sequences = unzip_file(filepath)
    if sequences is None:
        print("Problem in filepath, please double check")
        return ["didn't run", filepath, 0]

    for seq in sequences:
        if seq.id != read or windowSize is None:
            continue

        mean_s = []
        mean_e = []
        if len(seq) < maxlengthtelo:
            maxlengthtelo = len(seq)

        # check both ends of the read to find telo boundary 
        seq_start = str(seq.seq[trimfirst:maxlengthtelo]).upper()
        read_ids.append([seq.id, "forward"])

        seq_end = str(seq.seq[::-1]).upper()
        seq_end = seq_end[trimfirst:maxlengthtelo]  # flip chr end
        read_ids.append([seq.id, "reverse"])

        # cut sequence by windowSize and sliding
        windows_start = seq_cut_windows(seq_start, windowSize, slide)
        windows_end = seq_cut_windows(seq_end, windowSize, slide)

        # initiate matches 
        for start, seq_cut in windows_start:
            count_matches_read = [
                len([m.start() for m in pattern.finditer(seq_cut)]) or 1
                for pattern in compiled_patterns
            ]
            mean_s.append(("forward", start, sum(count_matches_read) / len(count_matches_read)))

        for start, seq_cut_2 in windows_end:
            count_matches_read = [
                len([m.start() for m in pattern.finditer(seq_cut_2)]) or 1
                for pattern in compiled_patterns
            ]
            mean_e.append(("reverse", start, sum(count_matches_read) / len(count_matches_read)))

        # specify the tail of read that have telomere 
        if tail == 'forward':
            mean_e = []
        elif tail == 'reverse':
            mean_s = []

        # check mean of end tail 
        def process_mean(mean, direction):
            if not mean:
                return

            x = [item[1] + trimfirst for item in mean if item[1] + trimfirst <= maxlengthtelo]
            y = [item[2] for item in mean if item[1] + trimfirst <= maxlengthtelo]

            if not x:
                return

            algo = rpt.Binseg(model="l2").fit(np.array(y))
            result = algo.predict(pen=4, n_bkps=1)
            all_cps = [x[cp] for cp in result[:-1]]

            if all_cps:
                telo_boundary_point = int(all_cps[0])
                if plot_yes_no:
                    plt.figure(figsize=(7.5, 3), dpi=300)
                    plt.plot(x, y, color='#000000', linestyle='-', linewidth=2)
                    plt.axvline(x=telo_boundary_point, color='#FF2C2C', linewidth=2, linestyle='--', label=f'x = boundary point: {telo_boundary_point}')
                    plt.title(f'mean window + boundary point of {seq.id}')
                    plt.xlabel('base pair (bp)')
                    plt.ylabel('mean window value')
                    if plotcp_range:
                        plt.xlim(0,plotcp_range)
                    else:
                        plt.xlim(0, maxlengthtelo)
                    plt.tight_layout()
                    plt.grid(True)

                if telo_boundary_point <= maxlengthtelo and telo_boundary_point != 0:
                    boundary.append([read, telo_boundary_point]) #save the boundary point
                else:
                    boundary.append([read, 0])

        process_mean(mean_e, "reverse")
        process_mean(mean_s, "forward")

    return boundary

def plot_patterns(seq, patterns, read_ids, added_labels, ax, direction):
    trans_table = str.maketrans('ACGT', 'TGCA')
    for patt in patterns:
        patterns.append(patt)
        comp_pattern = patt.translate(trans_table)
        patterns.append(comp_pattern)

    print("patterns:", patterns)

    for i, pattern in enumerate(patterns):
        pattern_name = pattern
        read_ids.append(pattern_name)
        matches = [m.start() for m in re.finditer(re.compile(pattern), seq)]
        if pattern not in added_labels:
            ax.scatter(matches, [(i + 1) * 2] * len(matches), color=colors[i], marker='|', label=pattern, zorder=2)
            added_labels.add(pattern)
        else:
            ax.scatter(matches, [(i + 1) * 2] * len(matches), color=colors[i], marker='|', zorder=2)

def rawCountPattern(filepath, read, pattern_telo, windowSize, slide, trimfirst, cut_length, minSeqLength, maxlengthtelo, tail=None, plot_raw=False):
    '''
    To get raw count of telomere-like repeat in the sequence
    filepath: location of fasta file
    read: read name, but just 1 read at the time 
    pattern: pattern wanna find, maybe ACCG. multiple patterns at each time
    windowSize: size of window
    slide: step of each window
    region: cut the sequence - in telomere, just check for first 5000 bp
    trimfirst: number of nucleotide trim off 

    return: raw count  
    '''

    # Loop through reads in file, find interested read. Then find their matches and plot 
    if not isinstance(read, str):
        print("can only read in 1 read at a time")
        return None

    patterns = patterns_to_search(pattern_telo, cut_length=cut_length)
    compiled_patterns = [re.compile(pattern) for pattern in patterns]

    sequences = unzip_file(filepath)
    if sequences is None:
        print("Problem in filepath, please double check")
        return ["didn't run", filepath, 0]

    rawcount_all = []

    for seq in sequences:
        if seq.id != read or windowSize is None:
            continue

        # print('working on:', seq.id)

        # check both ends of the read to find telo boundary 
        seq_start = str(seq.seq[trimfirst:maxlengthtelo]).upper()
        seq_end = str(seq.seq[::-1]).upper()[trimfirst:maxlengthtelo]

        windows_start = seq_cut_windows(seq_start, windowSize, slide)
        windows_end = seq_cut_windows(seq_end, windowSize, slide)

        rawpattern_s = [
            ("forward", start, pattern.pattern, len([m.start() for m in pattern.finditer(seq_cut)]) or 1)
            for start, seq_cut in windows_start
            for pattern in compiled_patterns
        ]

        rawpattern_e = [
            ("reverse", start, pattern.pattern, len([m.start() for m in pattern.finditer(seq_cut_2)]) or 1)
            for start, seq_cut_2 in windows_end
            for pattern in compiled_patterns
        ]

        if tail == 'forward':
            rawpattern_e = []
        elif tail == 'reverse':
            rawpattern_s = []

        rawcount_all.extend(rawpattern_s)
        rawcount_all.extend(rawpattern_e)

    plot_raw=False
    if (tail == 'forward' or tail is None) and plot_raw==True:
        sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})
        fig, ax = plt.subplots(figsize=(12, 8))
        read_ids = []
        added_labels = set()

        for seqraw in sequences:
            if len(seqraw.seq) > maxlengthtelo and seqraw.id == read:
                seq = str(seqraw.seq[:maxlengthtelo]).upper()
                plot_patterns(seq, patterns, read_ids, added_labels, ax, "forward")

        ax.set_title(f'Matches of telomere phrases in {read} forward')
        ax.set_xlabel('Position (bp)')
        ax.set_yticks([(i + 1) * 2 for i in range(len(read_ids))])
        ax.set_yticklabels(read_ids)
        plt.xlim(0, maxlengthtelo)
        plt.tight_layout()
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)
        plt.savefig(filename_s, format='png', dpi=300)

    if (tail == 'reverse' or tail is None) and plot_raw==True:
        sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})
        fig, ax = plt.subplots(figsize=(12, 8))
        read_ids = []
        added_labels = set()

        for seqraw in sequences:
            if len(seqraw.seq) > maxlengthtelo and seqraw.id == read:
                seq_2 = str(seqraw.seq[::-1][:maxlengthtelo]).upper()
                plot_patterns(seq_2, patterns, read_ids, added_labels, ax, "reverse")

        ax.set_title(f'Matches of telomere phrases in {read} reverse')
        ax.set_xlabel('Position (bp)')
        ax.set_yticks([(i + 1) * 2 for i in range(len(read_ids))])
        ax.set_yticklabels(read_ids)
        plt.xlim(0, maxlengthtelo)
        plt.tight_layout()
        ax.xaxis.grid(True)
        ax.yaxis.grid(True)
        plt.savefig(filename_e, format='png', dpi=300)

    return pd.DataFrame(rawcount_all, columns=['tail', 'position', 'pattern', 'count'])

# automate get TRC values and its corresponding telomere length
def fit_quadratic_and_find_vertex(trc_list, telo_length_list):
    """
    Fit a quadratic curve to TRC and telomere length data,
    and find the TRC value where the derivative is zero (vertex).
    """
    trc_arr = np.array(trc_list)
    telo_arr = np.array(telo_length_list)
    coeffs = np.polyfit(trc_arr, telo_arr, 2)
    a, b, c = coeffs
    vertex_x = -b / (2 * a)
    vertex_y = a * vertex_x**2 + b * vertex_x + c
    return vertex_x, vertex_y, coeffs