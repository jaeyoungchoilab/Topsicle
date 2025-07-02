# This script is used to run the Topsicle analysis on a set of fastq or fasta files, return telomere length. Main code.

import sys
import os
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO 
import time
import datetime
from collections import defaultdict

#for parallel 
import multiprocessing
from multiprocessing import Pool, Manager, Lock, cpu_count

import argparse

# result output
import seaborn as sns
import matplotlib.pyplot as plt
colors = sns.color_palette("colorblind",n_colors=30) 
sns.set_style("whitegrid", {'grid.color': 'grey', 'grid.linestyle': '--'})

import csv 

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# import our package here 
from Topsicle.allsteps import *

def get_log_path(args):
    # Use outputDir if available, else current directory
    log_dir = getattr(args, 'outputDir', '.')
    os.makedirs(log_dir, exist_ok=True)
    return os.path.join(log_dir, "topsicle_run.log")

def tprint(*args, **kwargs):
    msg = " ".join(str(a) for a in args)
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{now}] {msg}"
    print(line)
    # Write to log file if available
    if hasattr(tprint, "logfile"):
        with open(tprint.logfile, "a") as f:
            f.write(line + "\n")

verbose = False
def vprint(*args, **kwargs):
    if verbose:
        print(*args, **kwargs)

def process_file(args, seq_loc, telo_phrase, pattern,sliding_val,lock):
    tprint("subsetting raw dataset based on TRC cutoff")
    base_name = os.path.basename(seq_loc)
    file_name = os.path.splitext(base_name)[0]
    min_cutoff = min(args.cutoff) if isinstance(args.cutoff, (list, tuple)) else args.cutoff
    read_w_telo_mver = patternTRC_count(filepath=seq_loc, no_bp=1000, read_length=args.minSeqLength, telopattern=args.pattern, cutoff=min_cutoff, kmer=telo_phrase)
    
    read_ID_w_telo_mver = [item[0] for item in read_w_telo_mver]
    trc_read_w_telo_all = [item[3] for item in read_w_telo_mver]
    read_ID_w_telo_mver_tail = {item[0]: item[2] for item in read_w_telo_mver}

    # get reads that potentially have telomere
    fasta_temp = os.path.join(args.outputDir, f"{file_name}_trc_over_{min_cutoff}.fasta")
    if os.path.exists(fasta_temp):
        tprint(f"Temporary fasta file already exists: {fasta_temp}. Using existing file.")
    else:
        # Check if seq_loc is gzipped
        if seq_loc.endswith(".gz"):
            open_func = lambda f: gzip.open(f, "rt", encoding="utf-8")
            seq_format = "fastq" if seq_loc.endswith(".fastq.gz") or seq_loc.endswith(".fq.gz") else "fasta"
        else:
            open_func = lambda f: open(f, "rt")
            seq_format = "fastq" if seq_loc.endswith(".fastq") or seq_loc.endswith(".fq") else "fasta"

        # Decide output format and extension
        if seq_format == "fastq":
            out_format = "fastq"
            fasta_temp = os.path.join(args.outputDir, f"{file_name}_trc_over_{min_cutoff}.fastq")
        else:
            out_format = "fasta"
            fasta_temp = os.path.join(args.outputDir, f"{file_name}_trc_over_{min_cutoff}.fasta")

        with open_func(seq_loc) as in_handle, open(fasta_temp, "w") as out_handle:
            for record in SeqIO.parse(in_handle, seq_format):
                if record.id in read_ID_w_telo_mver:
                    SeqIO.write(record, out_handle, out_format)
        tprint(f"Temporary fasta file with TRC more than {min_cutoff}:", fasta_temp)

    bound_all_detected = []
    image_num = 1
    # in case wanting to check a specific read only 
    if args.read_check:
        tprint("checking specific read:", args.read_check)
        tail = read_ID_w_telo_mver_tail.get(args.read_check, None)
        tprint("step 2 on:", args.read_check)
        bound_res = bound_detect(filepath=fasta_temp, read=args.read_check, pattern_telo=pattern, windowSize=args.windowSize, tail=tail, cut_length=telo_phrase,
                                 slide=sliding_val, trimfirst=args.trimfirst, plot_yes_no=args.plot, maxlengthtelo=args.maxlengthtelo,plotcp_range=args.rangecp)
        
        readID, telolen = bound_res[0]

        try:
            idx = read_ID_w_telo_mver.index(args.read_check)
            trc_val = trc_read_w_telo_all[idx]
        except ValueError:
            trc_val = ""

        with open(f'{args.outputDir}/telolengths_all.csv', mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([file_name, telo_phrase, f"{trc_val:.3f}", readID,telolen])

        bound_all_detected.append((file_name, telo_phrase, bound_res, trc_val))
        
        if args.plot:
            plt.savefig(f"{args.outputDir}/plot_{telo_phrase}_{image_num}.png", format='png', dpi=300)
            plt.close()

        if args.rawcountpattern:
            allrawcount = rawCountPattern(filepath=fasta_temp, read=read, pattern_telo=pattern, windowSize=args.windowSize, maxlengthtelo=args.maxlengthtelo,
                                          slide=sliding_val, trimfirst=args.trimfirst, cut_length=telo_phrase, tail=tail, minSeqLength=args.minSeqLength, plot_raw=False)
            if args.rawcountpattern:
                allrawcount.to_csv(f"{args.outputDir}/rawcount_{telo_phrase}_{image_num}.csv")

    # check all filtered reads, no specifying read
    else:   
        for idx, read in enumerate(read_ID_w_telo_mver):
            tail = read_ID_w_telo_mver_tail.get(read, None)
            #tprint("step 2 on:", read)

            bound_res = bound_detect(filepath=fasta_temp, read=read, pattern_telo=pattern, windowSize=args.windowSize, tail=tail, cut_length=telo_phrase,
                                    slide=sliding_val, trimfirst=args.trimfirst, plot_yes_no=args.plot, maxlengthtelo=args.maxlengthtelo,plotcp_range=args.rangecp)


            readID, telolen = bound_res[0]
            trc_val = trc_read_w_telo_all[idx]
            with lock:
                with open(f'{args.outputDir}/telolengths_all.csv', mode='a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow([file_name, telo_phrase,f"{trc_val:.3f}", readID,telolen])

            bound_all_detected.append((file_name, telo_phrase, bound_res,trc_val))
            
            if args.plot:
                plt.savefig(f"{args.outputDir}/plot_{telo_phrase}_{image_num}.png", format='png', dpi=300)
                plt.close()

            if args.rawcountpattern:
                allrawcount = rawCountPattern(filepath=fasta_temp, read=read, pattern_telo=pattern, windowSize=args.windowSize, maxlengthtelo=args.maxlengthtelo,
                                            slide=sliding_val, trimfirst=args.trimfirst, cut_length=telo_phrase, tail=tail, minSeqLength=args.minSeqLength, plot_raw=False)
                if args.rawcountpattern:
                    allrawcount.to_csv(f"{args.outputDir}/rawcount_{telo_phrase}_{image_num}.csv")

            image_num += 1

    return bound_all_detected

def analysis_run(args):
    print("---- Topsicle run parameters ---")
    for k, v in vars(args).items():
        tprint(f"{k}: {v}")
    print("---------------------")

    tprint("Starting Topsicle analysis")

    manager = Manager()
    bound_all_detected = manager.list()
    os.makedirs(args.outputDir, exist_ok=True)  #make sure we have output directory

    if args.threads is not None:
        num_cores = args.threads
        tprint(f"Specified number of cores are/is: {num_cores}")
    else:
        try:
            num_cores = len(os.sched_getaffinity(0))
            tprint(f"By default, Topsicle allocates nmber of cores: {num_cores}")
        except AttributeError:
            num_cores = cpu_count()
            tprint(f"By default, Topsicle allocates number of cores: {num_cores}")

    output_csv = f'{args.outputDir}/telolengths_all.csv'
    tprint(f"Output will be here: {output_csv}")
    if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
        if args.override:
            tprint(f"Output file {output_csv} already exists and will be overridden becuz having --override flag.")
            os.remove(output_csv)
        else:
            tprint(f"Output file {output_csv} already exists and is not empty. Exiting to avoid overwrite. Use --override to force overwrite.")
            sys.exit(1)   

    if args.telophrase is None:
        telo_phrases= [len(args.pattern) - 2]
        tprint(f"No telophrase provided, use kmer: {telo_phrases}")
    else:
        telo_phrases = args.telophrase if isinstance(args.telophrase, list) else [args.telophrase]
    
        
    print("---------------------")

    with open(f'{args.outputDir}/telolengths_all.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['file_number', 'phrase','trc', "readID", 'telo_length'])
    
    bound_all_detected = []
    phrase_to_telo = defaultdict(list)
    phrase_to_trc = defaultdict(list)

    for telo_phrase in telo_phrases:
        if telo_phrase > len(args.pattern):
            tprint(f"Cannot have length of subset larger than length of pattern")
            tprint(f"Cannot get {telo_phrase}-bp cut from {len(args.pattern)}-bp pattern")
            sys.exit()
        
        if args.slide:
            sliding_val = args.slide
        else:
            sliding_val=len(args.pattern)
            
        pattern = patterns_to_search(telopattern=args.pattern, cut_length=telo_phrase)
        tprint("patterns to search:", pattern)

        filenames = []
        if not os.path.exists(args.outputDir):
            os.makedirs(args.outputDir)

        if os.path.isdir(args.inputDir):
            for root, dirs, files in os.walk(args.inputDir):
                for filename in files:
                    filenames.append(os.path.join(root, filename))
        else:
            filenames.append(args.inputDir)

        tprint("begin processing reads")
        with Manager() as manager:
            lock = manager.Lock()
            with Pool(processes=num_cores) as pool: 
                results = pool.starmap(process_file, [(args, seq_loc, telo_phrase, pattern,sliding_val,lock) for seq_loc in filenames])

        tprint("finished processing all reads")
        print("---------------------")

        for file_result in results:
            for entry in file_result:
                telophrase = entry[1]
                telolen = float(entry[2][0][1])
                trc_val = float(entry[3])
                phrase_to_telo[telophrase].append(telolen)
                phrase_to_trc[telophrase].append(trc_val)

        # Compute medians for each telo_phrase
    telo_phrases_sorted = sorted(phrase_to_telo)
    telo_phrase_medians = []
    median_telolen_list = []
    median_trc_list=[]

    if isinstance(args.cutoff, (list, tuple)):
        inputtrc = args.cutoff[0]
    else:
        inputtrc = args.cutoff

    for phrase in telo_phrases_sorted:
        median_telo = np.median(phrase_to_telo[phrase])
        median_trc = np.median(phrase_to_trc[phrase])
        telo_phrase_medians.append(phrase)
        median_telolen_list.append(median_telo)
        median_trc_list.append(median_trc)
        # telling TRC distribution and print obs 
        tprint(f"k-mer: {phrase}, with TRC >= {inputtrc}, median telomere length is {median_telo:.2f} bp")

        if len(phrase_to_telo[phrase]) >= 3:
            max_trc=max(phrase_to_trc[phrase])
            plot_path = os.path.join(args.outputDir, f"quadfit_{phrase}mer_{args.pattern}.png")
            vertex_x, vertex_y, coeffs = fit_quadratic_and_find_vertex(
                phrase_to_trc[phrase], phrase_to_telo[phrase], 
                inputtrc=inputtrc, median_trc=median_trc, save_path=plot_path)
            
            #vertex_x, vertex_y, coeffs = fit_quadratic_and_find_vertex(phrase_to_trc[phrase], phrase_to_telo[phrase])

            if vertex_x > max_trc:
                tprint(f"Asymptotic TRC {vertex_x:.3f} is greater than max TRC, which is not expected. See plot.")
                if median_trc < 1.0:
                    tprint(f"Using median TRC value ({median_trc:.3f}) as asymptotic TRC instead.")
                    vertex_x = median_trc
                else:
                    tprint(f"Using 0.9 as asymptotic TRC instead, since asymptotic is greater than 1.0.")
                    vertex_x = 0.9
            if vertex_x < 0.4:
                tprint("Quadratic fit suggests asymptotic TRC less than 0.4. See plot with fit line")
                if max_trc < 0.4:
                    tprint(f"Maximum TRC value in data is {max_trc:.3f}, which is less than 0.4, indicating low confidence in telomere detection.")
                if vertex_x < inputtrc:
                    tprint(f"Asymptotic TRC {vertex_x:.3f} is less than input cutoff {inputtrc:.3f}. Topsicle declares input TRC (={inputtrc}) as asymptotic TRC.")
                    vertex_x = inputtrc
            
            tprint(f"asymptotic TRC, or recommended cutoff: {vertex_x:.3f}")

            # Get telomere lengths where TRC > vertex_x
            filtered_telolen = [
                telo for trc, telo in zip(phrase_to_trc[phrase], phrase_to_telo[phrase])
                if trc >= vertex_x
            ]
            if filtered_telolen:
                median_filtered_telolen = np.median(filtered_telolen)
                tprint(f"Median telomere length for reads with TRC cutoff >= {vertex_x:.3f}: {median_filtered_telolen:.2f} bp")
            else:
                tprint(f"No read has TRC >= {vertex_x:.3f}, please double check the data or submit log to GitHub.")
        else:
            tprint("Not enough data points to recommend TRC cutoff.")

        
    return tprint("All telomere found, have a nice day.")

version_number = "1.0.0"
Topsicle_output_prefix = "Topsicle"

# if __name__ == "__main__":
def main():
    import sys
    start_time = time.time()
    parser = argparse.ArgumentParser(description='Topsicle - Telomere length estimation from long reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--inputDir','-i', type=str, help='Required, Path to the input file or directory', required=True)
    parser.add_argument('--outputDir','-o', type=str, help='Required, Path to the output directory', required=True)
    parser.add_argument('--pattern', type=str, help='Required, Telomere pattern, for example, human has TTAGGG', required=True)
    parser.add_argument('--minSeqLength', type=int, help='Minimum length required for long read, default = 9kbp', default=9000)
    parser.add_argument('--rawcountpattern', action='store_true', help='Output raw count of pattern abundance of each window')
    parser.add_argument('--telophrase', nargs='+', type=int, help=' k-mer of telomere pattern, can be 4, 5,... ')
    parser.add_argument('--cutoff',nargs='+', type=float, help='Threshold of TRC to be telomere, can be 0.4, 0.5,... ', default=0.7)
    parser.add_argument('--windowSize', type=int, help='Sliding window size', default=100)
    parser.add_argument('--slide', type=int, help=' Window sliding step, default is initial telomere length')
    parser.add_argument('--trimfirst', type=int, help='Trimming off first number of base pair to prevent adapter', default=100)
    parser.add_argument('--maxlengthtelo', type=int, help='Longest value can be for telomere or sequence', default=20000)
    parser.add_argument('--plot', action='store_true', help='Plot of changes in mean window and change point detected, boolean, presence=True')
    parser.add_argument('--rangecp', type=int, help='optional, set range of changepoint plot for visualization purpose, default is maxlengthtelo')
    parser.add_argument('--read_check', type=str, help='optional, to get telomere of a specific read')
    parser.add_argument('--override','-ov', action='store_true', help='Override output file')
    parser.add_argument('--threads','-t', type=int, help='Number of CPU cores to use (by default, use all cores there)', default=None)

    args = parser.parse_args()
    tprint.logfile = get_log_path(args)
    
    analysis_run(args)
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time(s): {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()