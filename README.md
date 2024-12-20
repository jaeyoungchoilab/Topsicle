# Pocky

This package is used to identify telomere boundary in long read sequencing (ONT and PacBio). 

It will get inputs with format as either fasta (or fasta.gz) or fastq (or fastq.gz), along with other parameters about telomere patterns (see section 2), then output with the length of each read within that input file in a .csv file, and supplemental plots. 

## 1. Getting started

Let's get started with cloning this package
### 1.1. From source (GitHub)

Make new environment for this package 
```bash
python3 -m venv pocky   # python version = 3.6.8
source ./pocky/bin/activate
# update pip if necessary 
pip install --upgrade pip

```

Cloning the package

```bash
git clone ##put the github link here # clone repo
cd Pocky

# after you have cloned the github, you can run
python3 main.py -h
```

### 1.2. Install requirements 
``` bash
pip install -e .
```

Make sure that we installed dependents in the **requirements.txt** files, but maybe need to install cython manually.

Also can mannually install those packages instead of using requirement file:

``` 
biopython==1.75
cython
matplotlib==3.3.4
matplotlib-inline==0.1.6
numpy==1.19.5
pandas==2.2.0
ruptures==1.1.9
seaborn==0.11.2
```
Pocky was developed in Python 3.6, but still works well in Python versions 3.10 and 3.12. 

## 2. Usage 
### 2.1. Pocky workflow

In Pocky, you can have an overview of where are telomere patterns within the sequence with overview_plot.py, or jump right into the main analysis with main.py as below. 

1. We have an initial pattern that we want to look for (for example, the telomere pattern of Arabidopsis thaliana Col-0 strand is "CCCTAAA"). Since this desired pattern has 7 base pairs (7-bp), and long read sequence methods (Oxford Nanopore Technologies, PacBio HiFi,...) can have random sequencing errors, identifying a k-mer (a subset) of that 7-bp pattern will be less specificity than finding whole 7-bp. Pocky generates phrases of that strand to the length we want, for example 4-bp or 5-bp from 7-bp (--telophrase). Let's call them "k-mer patterns". 
2. Optional step 0: Overview descriptive plot to see if input read has telomere or not by marking location of telomere pattern found in that read (section 2.2)
3. Step 1: If that read has telomere at 5' end, first 1000bp should also have telomere-phrase patterns there. Checking first 1kbp will have the proportion of telomere-phrase patterns (TRC value), and the read will be marked as having telomere if its TRC value is larger than the TRC cutoff value (--cutoff). 
4. Step 2: After identifying which read has telomere, Pocky finds how long is that telomere by sliding (--slide) through window (--windowSize), counting the mean of number of patterns found within that window, and returning the boundary point between telomere and non-telomere region as the window with the mean value drops significantly. 
5. Optional step 3: If we want to know what kmer-bp pattern (kmer < initial length of telomere pattern), we use the flag --rawcountpattern to return a .csv file with location of window start, pattern, number of pattern in that window.   

**Flags in Pocky (both overview_plot and main):**
```
  -h, --help            show this help message and exit
  --inputDir INPUTDIR   Path to the input folder directory
  --outputDir OUTPUTDIR
                        Path to the output folder directory
  --pattern PATTERN     Telomere pattern, in Mver, AAACCG
  --minSeqLength MINSEQLENGTH
                        Minimum of long read sequence, default = 9kbp
  --rawcountpattern     Print raw count of number of times see that pattern in
                        each window
  --telophrase TELOPHRASE [TELOPHRASE ...]
                        Step 1 - Length of telomere cut, can be 4 or 5 or so
                        on
  --cutoff CUTOFF       Step 1 - Cutoff of TRC value to have telomere
  --windowSize WINDOWSIZE
                        Step 2 - Window size for sliding
  --slide SLIDE         Step 2 - Window sliding step
  --trimfirst TRIMFIRST
                        Step 2 - Trimming off first number of base pair in
                        case of adapter
  --maxlengthtelo MAXLENGTHTELO
                        Step 2 - Longest value can be for telomere or sequence
  --plot                Step 2 - Plot of changes in mean window and change
                        point detected, boolean, presence=True
```

### 2.2. Descriptive plot - overview_plot.py
This code will output overview plot of locations of telomere pattern (or snippet of telomere pattern) in the sequence, run the code below:

The script can be run from the command line using the following command:

```bash
python3 overview_plot.py [-h] [--inputDir INPUTDIR] [--outputDir OUTPUTDIR] [--singlefilePath] [--pattern PATTERN]
                        [--minSeqLength MINSEQLENGTH]
```

An example code can be like this:

```bash
python3 overview_plot.py \
  --inputDir "$input_dir" \
  --outputDir "$output_dir" \
  --pattern $telo_pattern \
  --minSeqLength 9000 \
  --telophrase 4 \
  --recfindingpattern
```

### 2.3. Telomere length finding - main.py
The script will output a .csv file containing read ID and telomere length of each read for all the reads in that inputDir which passed basic filtering (>9kbp).

The script can be run from the command line as below, as an example of using all flags, which will output:
- a .csv file with file number, IDs of reads in that file, and telomere length (default, always output this)
- plots of mean window changes and boundary points for each read tail, either start or end tail or both (flag --plot).
- a .csv file of rawcount 


```bash
python3 main.py \
  --inputDir "$input_dir" \
  --outputDir "$output_dir" \
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

There will be a .csv output file with file number, IDs of reads in that file, and telomere length, and plots of mean window changes and boundary points if flag --plot was called.

Also, if you want to know what specific patterns contribute to the mean window changes, the flag --rawcountpattern will provide information of read ID

The full directory with code and results are in **Pocky_usage** folder. After running it, it will return a .csv file (telolength_all.csv) for all telomere length of reads in the input, descriptive plots, heatmaps and the mean window change visualizations.

## 3. Troubleshooting

### 3.1. The code runs but no output
#### Check pattern
The pattern input should be the pattern at 5 prime of the read. For example, in Col-0 A.thaliana, this pattern should be 'CCCTAAA' (telomere pattern found at 5 prime), not 'TTTAGGG' (telomere pattern found at 3 prime). 

#### Check flags and input 
Sometimes, input can be missing or in wrong format, and the code will not have any output then. Missing flag can be a reason for 

### 3.2. Run out of memory 
- Check if you printed out so many plots or not (--plot)
- Also double check the memory allowance 


