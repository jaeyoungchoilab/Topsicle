![TopsicleLogo](Topsicle/TopsicleLogo.png)

# Topsicle

This package is used to identify telomere boundaries and length using long-read sequencing data from ONT or PacBio platforms. 

Topsicle can analyze fasta or fastq data and outputs the length of each read within that input file in a .csv file and supplemental plots. 

## Table of contents

* [1. Getting started](#1-getting-started)
* [2. Running Topsicle](#2-running-topsicle)
  * [2.1.1 Quick example of running Topsicle](#211-quick-example-of-running-topsicle)
  * [2.1.2 Detailed explanation of running Topsicle](#212-detailed-explanation-of-running-topsicle)
  * [2.1.3 Explanation of output](#213-explanation-of-output)
  * [2.2. Plotting and visualization (Optional)](#22-plotting-and-visualization-optional)
  * [2.3. Topsicle workflow](#23-topsicle-workflow)
* [3. Troubleshooting](#3-troubleshooting)

## 1. Getting started

Topsicle is written in Python 3.6, but tested in Python 3.10 and Python 3.12 versions.

Let's get started with cloning this package from GitHub: 

### 1.1. From source (GitHub)

Make a new environment for this package 
```bash
python3 -m venv Topsicle   # minimum python version required is 3.6.8
source ./Topsicle/bin/activate
# update pip if necessary 
pip install --upgrade pip

```

Cloning the package [Topsicle](https://github.com/jaeyoungchoilab/Topsicle.git):

```bash
git clone https://github.com/jaeyoungchoilab/Topsicle.git # clone repo
```

### 1.2. Install requirements 
``` bash
cd Topsicle
pip install -e .
```

With upcoming [pip update](https://github.com/pypa/pip/issues/11457), cython might requires to be installed mannually. Also, to manually install dependencies:

``` 
biopython>=1.75
cython >=0.29.21 
matplotlib>=3.3.4
matplotlib-inline>=0.1.6
numpy>=1.22.4
pandas>=2.2.0
ruptures==1.1.9
seaborn>=0.11.2
```

Verify the installation:
```bash
topsicle --help

## can also call the main.py file to run Topsicle
# python3 Topsicle/main.py --help
```

## 2. Running Topsicle 

### 2.1.1: Quick example of running Topsicle

General example:
```bash
topsicle \
  --inputDir $input_dir \
  --outputDir $output_dir \
  --pattern $telo_pattern
```

Demo file example:
```bash
topsicle \
  --inputDir Topsicle_demo/data_col0_teloreg_chr \
  --outputDir Topsicle_demo/result_temp \
  --pattern AAACCCT
```
[Topsicle_demo](Topsicle_demo) contains A. thaliana Col-0 reads from chromosome 1R reference genome (TAIR10, GCF_000001735.4). 

Topsicle will output: 
- a .csv file (telolength_all.csv) with the telomere lengths of input reads
- a quadratic fit plot to predict optimal TRC threshold (Telomere Repeat Count statisitcs, reflect how confident Topsicle is in identifying telomere of a read) and corresponding telomere length.
- descriptive plots of the input data including heatmaps and the mean telomere k-mer counts per windows.

Note that when running Topsicle on your data, it might return more files than just 5 files as in the Demo folder.

### 2.1.2: Detailed explanation of running Topsicle

Detailed example run:
```bash
topsicle \
  --inputDir $input_dir \
  --outputDir $output_dir \
  --pattern $telo_pattern \
  --minSeqLength 9000 \
  --telophrase 4 \
  --cutoff 0.4 \
  --windowSize 100 \
  --slide 6 \
  --trimfirst 200 \
  --maxlengthtelo 20000 \
  --plot \
  --rawcountpattern \
  --threads 20 \
  --override
```

Explanation of each parameter (run topsicle --help):

```
Topsicle - Telomere length estimation from long reads

options:
  -h, --help            show this help message and exit
  --inputDir INPUTDIR, -i INPUTDIR
                        Required, Path to the input file or directory (default: None)
  --outputDir OUTPUTDIR, -o OUTPUTDIR
                        Required, Path to the output directory (default: None)
  --pattern PATTERN     Required, Telomere pattern, for example, human has TTAGGG (default: None)
  --minSeqLength MINSEQLENGTH
                        Minimum length required for long read, default = 9kbp (default: 9000)
  --rawcountpattern     Output raw count of pattern abundance of each window (default: False)
  --telophrase TELOPHRASE [TELOPHRASE ...]
                        k-mer of telomere pattern, can be 4, 5,... (default: None)
  --cutoff CUTOFF [CUTOFF ...]
                        Threshold of TRC to be telomere, can be 0.4, 0.5,... (default: 0.7)
  --windowSize WINDOWSIZE
                        Sliding window size (default: 100)
  --slide SLIDE         Window sliding step, default is initial telomere length (default: 6)
  --trimfirst TRIMFIRST
                        Trimming off first number of base pair to prevent adapter (default: 100)
  --maxlengthtelo MAXLENGTHTELO
                        Longest value can be for telomere or sequence (default: 20000)
  --plot                Plot of changes in mean window and change point detected, boolean, presence=True (default:
                        False)
  --rangecp RANGECP     optional, set range of changepoint plot for visualization purpose, default is maxlengthtelo
                        (default: None)
  --read_check READ_CHECK
                        optional, to get telomere of a specific read (default: None)
  --override, -ov       Override output file (default: False)
  --threads THREADS, -t THREADS
                        Number of CPU cores to use (by default, use all cores there) (default: None)
```

### 2.1.3 Explanation of output
Topsicle will output a .csv file containing the read ID and telomere length of all reads in the --inputDir that passed filtering.

#### Quick summary
Main outputs of interest.
- [telolengths_all.csv](Topsicle_demo/telolengths_all.csv): Output file with file number, IDs of reads in that file, and telomere length (default, always output this)
- [$output.fastq](Topsicle_demo/result_justone/Col-0-6909_GWHBDNP00000001.1_nano_right.fastq_trc_over_0.4.fastq): Reads that passed TRC threshold

Additional possible outputs based on flags: 
- [$read.png](Topsicle_demo/result_justone/plot_4_1.png): Plot showing mean telomere repeat count by window and the telomere-subtelomere boundary point for each read (flag **--plot**).
Example **mean window change plot** of a sequencing read:
![Mean window](Topsicle_demo/result_justone/plot_4_1.png)
The red line indicate the estimated telomere-subtelomere boundary point.

- [$read.csv](Topsicle_demo/result_justone/rawcount_4_1.csv): Raw count output used for calculating the sliding window and mean telomere repeat count (flag **--rawcountpattern**)

- log file with parameters Topsicle had and results

#### Detailed summary
Example output: [telolengths_all.csv](Topsicle_demo/telolengths_all.csv) 

Main output of Topsicle and updates in real time while Topsicle is running. 
- file_number: Name of the input file(s) in the directory
- phrase: The phase of the k-mer used for searching. By default, if the telomere pattern is 6-bp long, Topsicle will find 4-mer patterns (phrase = 4)
- trc: Telomere repeat count value of that read. This statistics is used for determining reads sequenced from the telomere (see the publication). 
- readID: ID of read
- telo_length: Estimated telomere length of read

Additional log file: [$output.log](Topsicle_demo/log_topsicle_demo.log)
- Information about resources used (number of cores, time, location of output)
- Real-time update
- Hard-choice TRC cutoff and median of telomere length if using this cutoff (line 11)
- Asymptotic TRC cutoff (line 12) and corresponding median telomere length (line 13). The asymptotic TRC is recommended if a hard-choice TRC cutoff can not be initially determined. 
=======

- [$read.csv](Topsicle_demo/result_justone/rawcount_4_1.csv): Raw count output used for calculating the sliding window and mean telomere repeat count (flag **--rawcountpattern**)

#### Detailed summary
Example output: [telolengths_all.csv](Topsicle_demo/telolengths_all.csv) 

Main output of Topsicle and updates in real time while Topsicle is running. 
- file_number: Name of the input file(s) in the directory
- phrase: The phase of the k-mer used for searching. By default, if the telomere pattern is 6-bp long, Topsicle will find 4-mer patterns (phrase = 4)
- trc: Telomere repeat count value of that read. This statistics is used for determining reads sequenced from the telomere (see the publication). 
- readID: ID of read
- telo_length: Estimated telomere length of read

Additional log file: [$output.log](Topsicle_demo/log_topsicle_demo.log)
- Information about resources used (number of cores, time, location of output)
- Real-time update
- Hard-choice TRC cutoff and median of telomere length if using this cutoff (line 11)
- Asymptotic TRC cutoff (line 12) and corresponding median telomere length (line 13). The asymptotic TRC is recommended if a hard-choice TRC cutoff can not be initially determined. 

If there is no line with "**All telomere found, have a nice day**" then Topsicle did not examine all possible reads in the raw data. The user can rerun the process or pick up the previous run by analyzing the smaller dataset containing reads that potentially have telomeres, called *Temporary fasta file*, as in line 8 of the demo log file. It is recommended to provide more resources and have a strict TRC cutoff value as well (any TRC > 0.6 will be strict). Also see section [3. Troubleshooting](#3-troubleshooting). 

### 2.2: Plotting and visualization of raw data (Optional)
Plot telomere k-mer matches in the sequencing read and a heatmap counting the different phases of the telomere k-mer.

As a note, this option is not developed to be called directly yet, we still need to call it using python3 /PATH/overview_plot.py as below:

```bash
python3 overview_plot.py \
  --inputDir "$input_dir" \
  --outputDir "$output_dir" \
  --pattern $telo_pattern \
  --minSeqLength 9000 \
  --telophrase 4 \
  --recfindingpattern \
  --rawcount
```

**Descriptive plot** displaying telomere repeat counts in the first 40 reads: 

![Descriptive plot](Topsicle_demo/result_justone/descriptive_plot_1.png)

**Heatmap** displaying telomere repeat count by telomere k-mers of different phases:

![Heatmap](Topsicle_demo/result_justone/heatmap_1.png)


### 2.3: Topsicle workflow 

1. We have a telomere pattern that we want to identify length (for example, the telomere pattern of Arabidopsis thaliana Col-0 strand is "CCCTAAA"). Since the initial telomere pattern has 7 base pairs (7-bp) and long read sequence methods (Oxford Nanopore Technologies, PacBio HiFi,...) can have  sequencing errors, identifying a k-mer (a subset) of that 7-bp pattern will be less specific than finding the whole 7-bp. Topsicle generates k-mers (or phrases, or subsets) of that pattern, for example, 4-bp or 5-bp from a 7-bp (--telophrase). Let's call them "k-mer patterns" and perform analysis on them. 

Also a note, with this strategy, if user does not 100% sure of the telomere pattern of a species, they can provide the most abundant pattern from tools such as Tandem Repeats Finder, tidk, or by checking reads that align to the chromosome ends of the reference genome, which are likely to have telomere. In short, Topsicle needs a telomere pattern, but does not require having too accurate ones. 

User can provide multiple telomere patterns by using the pipe ("|"), such as "AAACCCG|AAACCG", but this feature is not fully tested in this version yet. 

2. [Optional visualization step](#22-plotting-and-visualization-optional): This step provides an initial visual analysis of the input sequences.

This optional step is used to see if the input reads contain the telomere repeat. It will return plots with the location of the telomere pattern found in that read. A heatmap can also be generated for observation of tandem repeats that we can expect to have in the dataset. By default, it will return a description plot of positions of k-mer patterns for a read, and in addition, using -recfindingpattern option will generate a heatmap (k-mer profile) and let us know which k-mer from the telomere pattern is abundant and what errors are commonly found within this read. Along with the heatmap, the parameter -rawcount will extract the count of each k-mer and its match in case we want to know exactly how many matches are at each position. 

It is advised to run this supplemental function prior to running the main function to have an overview of observations of species and their reads before finding the telomere length.  

**Executing Topsicle**

3. Step 1. TRC filtering: If read is sequenced from the telomere, the first 1000bp of that read should contain telomere repeat k-mers. The count of this initial telomere repeat k-mers (i.e., TRC statistics) will be used for filtering candidate telomere sequencing reads. Reads with a threshold TRC value (--cutoff) will be analyzed for downstream. In case a threshold can't be determined, Topsicle can calculate an automatic threshold using the asymptotic method (see manuscript for details) and calculate the telomere length.

4. Step 2. Telomere length calculation: After identifying potential read that has telomere, Topsicle finds how long is that telomere by sliding (--slide) through window (--windowSize) and measuring the mean of number of patterns found within that window, and returning the boundary point between telomere and non-telomere regions based on changepoint algorithm.

5. Optional step 3: If we want to know what kmer-bp pattern is most found within a window (kmer has to smaller than initial length of telomere pattern), we use the flag --rawcountpattern to return a .csv file with position of window start, pattern, and number of pattern found in that window. 

7. See [2.1.3 Explanation of output](#213-explanation-of-output) for output explanations

## 3. Troubleshooting

### 3.1. The code runs but no output
1. Check pattern:

The pattern input is recommended to be the pattern at the 5' end of the read. For example, in Col-0 A. thaliana, this pattern should be 'CCCTAAA', which is different in human and Mimulus. Getting the heatmap profile is also recommended.

2. Check flags and input:

Sometimes, input can be missing or in the wrong format, and the code will not have any output then. A missing flag can be a reason for not being able to run as well. Providing a higher value of TRC (such as TRC > 0.6) can help, too. 

### 3.2. Run out of memory 
1. Check if you printed out so many plots or not (with the flag --plot)
2. Double-check the memory allowance 

### 3.3: Not enough resources
This issue usually appears when running whole genome analysis but using fewer than 6 cores in 24 hours for testing files that are more than 20GB and/or more than 1 million reads (observations based on testing trials on KU HPC)
1. It is recommended to have more resources allocated - maybe more core, more time, or both. If possible, breaking down the file into several 1GB and/or 0.2 million reads files, then submitting several jobs and putting them together can also help. 
2. If the analysis keeps cancelling after several attempts, please keep in mind that even though Topsicle returns some results in the telolength_all.csv file, it might not look through and contain results from every read with telomere and their length. However, this file should provide some information. We also recommend submitting an issue request on GitHub if this keeps happening.

**Note:** It will echo **"All telomere found, have a nice day."** when Topsicle checked all reads in the dataset. It is highly recommended to have a log file and look for this line when running Topsicle on a big dataset.
