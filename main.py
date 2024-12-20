# This script is used to run the BoundTeloNano analysis on a set of fastq or fasta files.

import sys
import os
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO 
# os.environ['MPLCONFIGDIR'] = os.path.join(os.path.dirname(__file__), 'config')

#for parallel 
from multiprocessing import Pool
from multiprocessing import cpu_count

import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# for validate performance 
# import cProfile
# import pstats
# import io

# Add project root to sys.path 
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
colors = sns.color_palette("colorblind",n_colors=30)  
# import our package here 
from Pocky.allsteps import *

sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

verbose = False
def vprint(*args, **kwargs):
    if verbose:
        print(*args, **kwargs)

#check number of cores 
num_cores = cpu_count()
print(f"Number of cores: {num_cores}")

def process_file(args, seq_loc, telo_phrase, pattern):
    base_name = os.path.basename(seq_loc)
    file_name = os.path.splitext(base_name)[0]
    # print("Working on:", seq_loc)
    
    read_w_telo_mver = patternTRC_count(filepath=seq_loc, no_bp=1000, read_length=args.minSeqLength, telopattern=args.pattern, cutoff=args.cutoff, kmer=telo_phrase)
    
    read_ID_w_telo_mver = [item[0] for item in read_w_telo_mver]
    read_ID_w_telo_mver_tail = {item[0]: item[2] for item in read_w_telo_mver}

    bound_all_detected = []
    image_num = 1
    for read in read_ID_w_telo_mver:
        tail = read_ID_w_telo_mver_tail.get(read, None)
        print("step 2 on:", read)

        bound_res = bound_detect(filepath=seq_loc, read=read, pattern_telo=pattern, windowSize=args.windowSize, tail=tail, cut_length=telo_phrase,
                                 slide=args.slide, trimfirst=args.trimfirst, plot_yes_no=args.plot, maxlengthtelo=args.maxlengthtelo)
        bound_all_detected.append((file_name, telo_phrase, bound_res))
        
        if args.plot:
            plt.savefig(f"{args.outputDir}/plot_{telo_phrase}_{image_num}.png", format='png', dpi=300)
            plt.close()

        if args.rawcountpattern:
            allrawcount = rawCountPattern(filepath=seq_loc, read=read, pattern_telo=pattern, windowSize=args.windowSize, maxlengthtelo=args.maxlengthtelo,
                                          slide=args.slide, trimfirst=args.trimfirst, cut_length=telo_phrase, tail=tail, minSeqLength=args.minSeqLength, plot_raw=False)
            if args.rawcountpattern:
                allrawcount.to_csv(f"{args.outputDir}/rawcount_{telo_phrase}_{image_num}.csv")

        image_num += 1

    return bound_all_detected

def analysis_run(args):
    print("Begin running")
    telo_phrases = args.telophrase
    bound_all_detected = []

    for telo_phrase in telo_phrases:
        if telo_phrase > len(args.pattern):
            print(f"Cannot have length of subset larger than length of pattern")
            print(f"Cannot get {telo_phrase}-bp cut from {len(args.pattern)}-bp pattern")
            sys.exit()

        pattern = patterns_to_search(telopattern=args.pattern, cut_length=telo_phrase)
        print("patterns to search:", pattern)

        filenames = []
        if not os.path.exists(args.outputDir):
            os.makedirs(args.outputDir)

        if os.path.isdir(args.inputDir):
            for root, dirs, files in os.walk(args.inputDir):
                for filename in files:
                    filenames.append(os.path.join(root, filename))
        else:
            filenames.append(args.inputDir)

        print("begin processing reads")
        with Pool(processes=num_cores) as pool: 
            results = pool.starmap(process_file, [(args, seq_loc, telo_phrase, pattern) for seq_loc in filenames])

        print("start writing results")
        for result in results:
            bound_all_detected.extend(result)

    flat_data = []
    for record in bound_all_detected:
        file_name = record[0]
        phrases = record[1]
        nested_data = record[2]

        for sub_record in nested_data:
            flat_data.append([file_name] + [phrases] + sub_record)

    bound_all_detected_df = pd.DataFrame(flat_data, columns=['file_number', 'phrase', 'readID', 'telo_length'])
    bound_all_detected_df.to_csv(f"{args.outputDir}/telolengths_all.csv")
    print(f"output: {args.outputDir}/telolengths_all.csv")
    return bound_all_detected_df

version_number = "1.0.0"
BoundTeloNano_output_prefix = "Pocky"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Command line input handling for run_analysis function')
    parser.add_argument('--inputDir', type=str, help='Path to the input folder directory', required=True)
    parser.add_argument('--outputDir', type=str, help='Path to the output folder directory', required=True)
    parser.add_argument('--pattern', type=str, help='Telomere pattern, in Mver, AAACCG', required=True)
    parser.add_argument('--minSeqLength', type=int, help='Minimum of long read sequence, default = 9kbp', default=9000)
    parser.add_argument('--rawcountpattern', action='store_true', help='Print raw count of number of times see that pattern in each window')
    parser.add_argument('--telophrase', nargs='+', type=int, help=' Step 1 - Length of telomere cut, can be 4 or 5 or so on', default=4)
    parser.add_argument('--cutoff', type=float, help='Step 1 - Cutoff of TRC value to have telomere', default=0.4)
    parser.add_argument('--windowSize', type=int, help=' Step 2 - Window size for sliding', default=100)
    parser.add_argument('--slide', type=int, help=' Step 2 - Window sliding step', default=6)
    parser.add_argument('--trimfirst', type=int, help='Step 2 - Trimming off first number of base pair in case of adapter', default=100)
    parser.add_argument('--maxlengthtelo', type=int, help='Step 2 - Longest value can be for telomere or sequence', default=20000)
    parser.add_argument('--plot', action='store_true', help='Step 2 - Plot of changes in mean window and change point detected, boolean, presence=True')

    args = parser.parse_args()
    analysis_run(args)

    # for checking performance
    # pr = cProfile.Profile()
    # pr.enable()

    # # Code to profile
    # analysis_run(args)
    # pr.disable()

    # # Print profiling results
    # s = io.StringIO()
    # ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
    # ps.print_stats()
    # print(s.getvalue())