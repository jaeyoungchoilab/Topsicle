# This script is used to run the Topsicle analysis on a set of fastq or fasta files.
# it will return a lot of files - descriptive plot, heatmap and mean window plots 
# it's not fully developed yet so user needs to run as 
# python3 /PATH/overview_plot.py --commands

import sys
import os
import gzip
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO 
from itertools import chain
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import time
import datetime
from collections import defaultdict

# Add project root to sys.path - so will not get error in importing code
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
colors = sns.color_palette("colorblind",n_colors=30)

from Topsicle.descriptive_plot import *
sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

version_number = "1.0.0"
Topsicle_output_prefix = "Topsicle"

def tprint(*args, **kwargs):
    msg = " ".join(str(a) for a in args)
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}]", msg)

def plot_running(args):
    # begin plotting
    filenames = []
    # print("Data type: Multiple files in the directory")

    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir)

    if os.path.isdir(args.inputDir):
        for root, dirs, files in os.walk(args.inputDir):
            for filename in files:
                filenames.append(os.path.join(root, filename))
    else:
        filenames.append(args.inputDir)

    # in case user not provide telomere phrase
    if args.telophrase is None:
        telo_phrases= [len(args.pattern) - 2]
        tprint(f"No telophrase provided, use kmer: {telo_phrases}")
    else:
        telo_phrases = args.telophrase if isinstance(args.telophrase, list) else [args.telophrase]
    
            
    print("Loaded all data, start plotting")

    # descriptive plotting 
    for i, seq_loc in enumerate(filenames, start=1):
        print(f"Descriptive plot on: {seq_loc}")

        # descriptive plot
        descriptive_plot(seq_loc, pattern=args.pattern, minSeqLength=args.minSeqLength)
        plt.savefig(f"{args.outputDir}/descriptive_plot_{i}.png", format='png', dpi=300)
        plt.close()
    
    print(f"Descriptive plot is in here: {args.outputDir}")

    # heatmap plotting
    if args.recfindingpattern:
        for i, seq_loc in enumerate(filenames, start=1):
            for phrase in telo_phrases:
                print(f"Heatmap on {seq_loc}")
                heatmap = patterns_vs_match_heatmap(seq_loc, args.pattern, phrase, args.minSeqLength)
                plt.savefig(f"{args.outputDir}/heatmap_{i}.png", format='png', dpi=300)
                plt.close()

                if args.rawcount:  # in case user wants to save the raw count to replot the heatmap in another way
                    csv_path = f"{args.outputDir}/heatmap_rawcount_{i}.csv"
                    print(f"Saving raw count of heatmap to {csv_path}")
                    heatmap.to_csv(csv_path, index=False)

    
    print(f"Heatmap is in here: {args.outputDir}")

    return 'plotted the plot'

if __name__ == "__main__":
   # Create the argument parser
    parser = argparse.ArgumentParser(description='Command line input handling for run_analysis function')

    # Add the arguments
    parser.add_argument('--inputDir', type=str, help='Path to the input folder directory')
    parser.add_argument('--outputDir', type=str, help='Path to the output folder directory')
    parser.add_argument('--pattern', type=str, help='Telomere pattern, in Mver, AAACCG')
    parser.add_argument('--minSeqLength', type=int, help='Minimum of long read sequence, default = 9kbp', default=9000)
    parser.add_argument('--telophrase', nargs='+', type=int, help='Step 1 - Length of telomere cut, can be 4 or 5 or so on')
    parser.add_argument('--recfindingpattern', action='store_true', help='Boolean, use this to plot the heatmap of patterns vs match')
    parser.add_argument('--rawcount', action='store_true', help='Boolean, save raw count results to CSV for flexibility of plotting')


    # Parse the command line arguments
    args = parser.parse_args()
    plot_running(args)
