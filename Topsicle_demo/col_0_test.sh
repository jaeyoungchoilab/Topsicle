#!/bin/bash
#SBATCH --job-name=test_col0_%j  # Job name
#SBATCH --partition=sixhour     # Partition Name 
#SBATCH --ntasks=4              # Run on a single CPU
#SBATCH --mem=10gb               # Job memory request
#SBATCH --time=6:00:00       # Time limit hrs:min:sec
#SBATCH --output=test_col0_%j.log    # Standard output and error log
#SBATCH --error test_col0_%j.err

# activate environment 
environment_name=/kuhpc  # replace with your environment name
source $environment_name/topsicle/bin/activate  # replace your actual topsicle path

# initial directories and patterns
input_dir=Topsicle_demo/data_col0_teloreg_chr   # can input the file or the directory, replace with your actual data path
output_dir=Topsicle_demo/result_all             # need to specific the output directory
telo_pattern=CCCTAAA                                      # telomere pattern                    

# 1. running main.py

# Output: telolengths_all.csv for telomere in all reads 
# also mean window plots, .csv files of raw count of each kmer pattern

python3 main.py \
  --inputDir $input_dir \
  --outputDir $output_dir \
  --pattern $telo_pattern \
  --plot \
  --rawcountpattern

# overview_plot.py run 
# output: heatmap of patterns and the matchs after it 
python3 overview_plot.py \
  --inputDir "$input_dir" \
  --outputDir "$output_dir" \
  --pattern $telo_pattern \
  --recfindingpattern \
  --rawcount

# then, submit sbatch col_0_test.sh, or run it directly in terminal (see log of running in terminal in log_topsicle_demo.log)
# output will be in the output_dir directory


# In case you curious, below are all possible parameters:

# 1. running main.py

# Output: telolengths_all.csv for telomere in all reads 
# also mean window plots, .csv files of raw count of each kmer pattern

# python3 main.py \
#   --inputDir $input_dir \
#   --outputDir $output_dir \
#   --pattern $telo_pattern \
#   --minSeqLength 9000 \
#   --telophrase 4 \
#   --cutoff 0.4 \
#   --windowSize 100 \
#   --slide 6 \
#   --trimfirst 200 \
#   --maxlengthtelo 20000 \
#   --plot \
#   --rawcountpattern

# # overview_plot.py run 
# # output: heatmap of patterns and the matchs after it 
# python3 overview_plot.py \
#   --inputDir "$input_dir" \
#   --outputDir "$output_dir" \
#   --pattern $telo_pattern \
#   --minSeqLength 9000 \
#   --telophrase 5 \
#   --recfindingpattern \
#   --rawcount
