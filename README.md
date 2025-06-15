# Topsicle

This package is used to identify telomere boundary in long read sequencing (ONT and PacBio). 

Topsicle gets inputs in format as either fasta (or fasta.gz) or fastq (or fastq.gz), along with other parameters about telomere patterns (see section 2), then output with the length of each read within that input file in a .csv file and supplemental plots. 

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

Make new environment for this package 
```bash
python3 -m venv Topsicle   # minimum python version required is 3.6.8
source ./Topsicle/bin/activate
# update pip if necessary 
pip install --upgrade pip

```

Cloning the package [Topsicle](https://github.com/jaeyoungchoilab/Topsicle.git):

```bash
git clone https://github.com/jaeyoungchoilab/Topsicle.git # clone repo
cd Topsicle

# verify the cloning process: 
python3 main.py -h
```

### 1.2. Install requirements 
``` bash
pip install -r requirements.txt
```

Make sure that we installed dependents in the **requirements.txt** files, but user might need to install cython manually.

To mannually install those packages instead of using requirement file:

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
Topsicle was developed in Python 3.6, but also passed tests in Python versions 3.10 and 3.12. 

## 2. Running Topsicle 

Using Topsicle, you can have an overview of where are telomere patterns within the sequence with [overview_plot.py](#22-descriptive-plots---overview_plotpy), or jump right into the **main analysis** to get telomere lengths in reads using [main.py](#21-telomere-length-finding---mainpy).  

The full directory with code and results are in [Topsicle_demo](Topsicle_demo). After running it, it will return a .csv file (telolength_all.csv) for all telomere length of reads in the input (Topsicle will write the result to this .csv file in real time), descriptive plots, heatmaps and mean window change visualizations.

### 2.1.1: Quick example of running Topsicle

General example:
```bash
python3 main.py \
  --inputDir $input_dir \
  --outputDir $output_dir \
  --pattern $telo_pattern
```

Demo file example:
```bash
python3 main.py \
  --inputDir Topsicle_demo/data_col0_teloreg_chr \
  --outputDir Topsicle_demo/result_justone \
  --pattern AAACCCT
```
[Topsicle_demo](Topsicle_demo) contains A. thaliana Col-0 reads from chromosome 1R reference genome (TAIR10, GCF_000001735.4). Note when running Topsicle on your own data, it might return more files than just 5 files as in the Demo folder. 

### 2.1.2: Detailed explanation of running Topsicle

Detailed example run:
```bash
python3 main.py \
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
  --rawcountpattern
```

Explanations of options from **main.py**: 
```
options:
  -h, --help            show this help message and exit
  --inputDir INPUTDIR   Required, Path to the input folder directory
  --outputDir OUTPUTDIR
                        Required, Path to the output folder directory
  --pattern PATTERN     Required, Telomere pattern, in human, TTAGGG
  --minSeqLength MINSEQLENGTH
                        Minimum of long read sequence, default = 9kbp
  --rawcountpattern     Print raw count of number of times see that pattern in each window
  --telophrase TELOPHRASE [TELOPHRASE ...]
                        Step 1 - Length of telomere cut, can be 4 or 5 or so on
  --cutoff CUTOFF [CUTOFF ...]
                        Step 1 - Cutoff of TRC value to have telomere, can be 0.4, 0.5 or so on
  --windowSize WINDOWSIZE
                        Step 2 - Window size for sliding
  --slide SLIDE         Step 2 - Window sliding step, default is initial telomere length
  --trimfirst TRIMFIRST
                        Step 2 - Trimming off first number of base pair in case of adapter
  --maxlengthtelo MAXLENGTHTELO
                        Step 2 - Longest value can be for telomere or sequence
  --plot                Step 2 - Plot of changes in mean window and change point detected, boolean, presence=True
  --rangecp RANGECP     optional, set range of changepoint plot for visualization purpose, default is maxlengthtelo
  --read_check READ_CHECK
                        optional, to get telomere of a specific read
```

Example **mean window change plot** of a sequencing read:

![Mean window](Topsicle_demo/result_justone/plot_4_1.png)


### 2.1.3 Explanation of output
Topsicle will output a .csv file containing read ID and telomere length of all reads in the inputDir that passed filtering.

#### Quick summary
Main outputs of interest.
- [$output.csv file](Topsicle_demo/telolengths_all.csv): Output file with file number, IDs of reads in that file, and telomere length (default, always output this)
- [$output.fastq file](Topsicle_demo/result_justone/Col-0-6909_GWHBDNP00000001.1_nano_right.fastq_trc_over_0.4.fastq): Reads that passed TRC threshold

Additional possible outputs based on flags: 
- [plots](Topsicle_demo/result_justone/plot_4_1.png): Plot showing telomere repeat mean count by window and the telomere-subtelomere boundary point for each read (flag **--plot**).
- [$rawcount.csv file](Topsicle_demo/result_justone/rawcount_4_1.csv): Raw count output used for calculating the sliding window and mean telomere repeat count (flag **--rawcountpattern**)

#### Detailed summary
Example output: [telolengths_all.csv](Topsicle_demo/telolengths_all.csv) 

Main output of Topsicle and updates in real time while Topsicle is running. 
- file_number: name of the input file(s) in the directory
- phrase: k-mer value. By default, if the telomere pattern is 6-bp long, Topsicle will find 4-mer patterns (phrase = 4)
- trc: Telomere repeat count value of that read, which seperate between reads with and without telomere (see the publication). Higher trc value, more confident the algorithm can be in the identified telomere length, but this value is just for reference and prediction of optimized TRC value. 
- readID: ID of read has telomere
- telo_length: Length of telomere

Also, there will be the [log file](Topsicle_demo/log_topsicle_demo.log), which contains:
- Information of resources used: number of cores, time, location of output
- Real time update
- Hard-choice TRC cutoff and median of telomere length if using this cutoff (line 11)
- In case a predicted TRC cutoff is desired, the asymptotic TRC cutoff (line 12) and corresponding median telomere length (line 13) can also be found. 

If there is no line "**All telomere found, have a nice day**", that means Topsicle did not examine all possible reads in the raw data. User can rerun the process or pick up previous run by supplying *Temporary fasta file*, as in line 8 of the demo log file with its path, and please supply more resources. Also see section [3. Troubleshooting](#3-troubleshooting). 


### 2.2: Plotting and visualization (Optional)
This code will output overview plot of locations of telomere pattern (or snippet of telomere pattern) in the sequence and heatmap. Run the code as below:

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

**Heatmap** displaying telomere repeat count by k-mer:

![Heatmap](Topsicle_demo/result_justone/heatmap_1.png)


### 2.3: Topsicle workflow 
A more detailed version is in the publication. 
1. We have an telomere pattern that we want to identify length (for example, the telomere pattern of Arabidopsis thaliana Col-0 strand is "CCCTAAA"). Since initial telomere pattern has 7 base pairs (7-bp) and long read sequence methods (Oxford Nanopore Technologies, PacBio HiFi,...) can have  sequencing errors, identifying a k-mer (a subset) of that 7-bp pattern will be less specificity than finding whole 7-bp. Topsicle generates k-mer (or phrase, or subset) of that pattern, for example 4-bp or 5-bp from 7-bp (--telophrase). Let's call them "k-mer patterns" and perform analysis on them. 

Also a note, with this strategy, if user does not 100% sure of telomere pattern of a species, they can provide the most abundance pattern from tools such as Tandem Repeats Finder, tidk or by checking reads align to chromosome ends of reference genome, which are likely to have telomere. In short, Topsicle needs telomere pattern, but does not require to have too accurate ones. 

User can provide multiple telomere patterns by using the pipe ("|"), such as "AAACCCG|AAACCG", but this feature is not fully tested in this version yet. 

2. [Optional visualization step](#22-plotting-and-visualization-optional): This step provides an initial visual analysis of the input sequences.

This optional step is used to see if the input reads contains the telomere repeat. It will return plots with location of telomere pattern found in that read. Heatmap can also be generated for observation of tandem repeats we can expect to have in the dataset. By default, it will return a description plot of positions of k-mer patterns for a read, and in addition, using -recfindingpattern option will generate a heatmap (k-mer profile) and let us know what k-mer from the telomere pattern is abundance and what errors are commonly found within this read. Along with the heatmap, the parameter -rawcount will extract the count of each k-mer and its match in case we want to know exact how many matches at each position. 

It is advised to run this supplemental function prior to running the main function to have the overview observations of species and its reads before finding telomere length.  

**Executing Topsicle**

3. Step 1. TRC filering: If read is sequenced from the telomere the first 1000bp of that read should contain telomere repeat k-mers. The count of this initial telomere repeat k-mers (i.e. TRC statistics) will be used for filtering candidate telomere sequencing reads. Reads with larger than a threshold TRC value (--cutoff) will be analyzed for downstream. In case a threshold can't be determined Topsicle can calculate an automatic threshold using the asymptotic method (see manuscript for detail) and calculate the telomere length.

4. Step 2. Telomere length calculation: After identifying potential read that has telomere, Topsicle finds how long is that telomere by sliding (--slide) through window (--windowSize) and measuring the mean of number of patterns found within that window, and returning the boundary point between telomere and non-telomere regions based on changepoint algorithm.

5. Optional step 3: If we want to know what kmer-bp pattern is most found within a window (kmer has to smaller than initial length of telomere pattern), we use the flag --rawcountpattern to return a .csv file with position of window start, pattern, and number of pattern found in that window. 

7. See [2.1.3 Explanation of output](#213-explanation-of-output) for output explanations

## 3. Troubleshooting

### 3.1. The code runs but no output
1. Check pattern:

The pattern input is recommended to be the pattern at 5 prime of the read. For example, in Col-0 A.thaliana, this pattern should be 'CCCTAAA', which is different in human and Mimulus. Getting the heatmap profile is also recommended as well.

2. Check flags and input:

Sometimes, input can be missing or in wrong format, and the code will not have any output then. Missing flag can be a reason for not being able to run as well.

### 3.2. Run out of memory 
1. Check if you printed out so many plots or not (with the flag --plot)
2. Double check the memory allowance 

### 3.3: Not enough resources
This issue usually appears when running whole genome analysis but use less than 6 cores in 24 hours for testing file that is more than 20GB and/or more than 1 million reads (observations based on testing trials on KU HPC)
1. It is recommended to have more resources allocate - maybe more core, more time or both. If possible, breaking down the file into several 1GB and/or 0.2 millions reads files, then submit several jobs and put them together can also help. 
2. If the analysis keeps cancelled after several attempts, please keep in mind that even though Topsicle returns some results in the telolength_all.csv file, it might not look through and contain results from every reads with telomere and their length. However, this file should provide some information. We also recommend to submit an issue request on GitHub if this keeps happening.

**Note:** It will echo **"All telomere found, have a nice day."** when Topsicle checked all reads in the dataset. It is highly recommended to have a log file and look for this line when running Topsicle on big dataset.


