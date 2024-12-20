# This script is used to run the Pocky analysis on a set of fastq or fasta files.
# it will return a lot of files

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

# Add project root to sys.path - so will not get error in importing nov5_allcode.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
colors = sns.color_palette("colorblind",n_colors=30)

from Pocky.descriptive_plot import *
sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

version_number = "1.0.0"
Pocky_output_prefix = "Pocky"

def plot_running(args):
    # # checking inputs first 
    # if os.path.isdir(args.inputDir):
    #     print(f"Problem with input for {args.inputDir}")
    #     print("If you only have 1 file, please get the path of that file, not the directory.")
    #     print("If you have multiple files in a directory, please delete singlefilePath.")
    #     sys.exit()

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
            for phrase in args.telophrase:
                print(f"Heatmap on {seq_loc}")
                patterns_vs_match_heatmap(seq_loc, args.pattern, phrase, args.minSeqLength)
                plt.savefig(f"{args.outputDir}/heatmap_{i}.png", format='png', dpi=300)
                plt.close()
    
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
    parser.add_argument('--telophrase', nargs='+', type=int, help='Step 1 - Length of telomere cut, can be 4 or 5 or so on', default=4)
    parser.add_argument('--recfindingpattern', action='store_true', help='Boolean, use this to plot the heatmap of patterns vs match')


    # Parse the command line arguments
    args = parser.parse_args()
    plot_running(args)
