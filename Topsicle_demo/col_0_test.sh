#!/bin/bash
#SBATCH --job-name=test_col0  # Job name
#SBATCH --partition=jychoi     # Partition Name (Required)
#SBATCH --ntasks=4              # Run on a single CPU
#SBATCH --mem=10gb               # Job memory request
#SBATCH --time=36:00:00       # Time limit hrs:min:sec
#SBATCH --output=test_col0_%j.log    # Standard output and error log
#SBATCH --error test_col0_%j.err

# activate environment 
source /topsicle/bin/activate

input_dir=/topsicle/topsicle_demo/data_col0_teloreg_chr
output_dir=/topsicle/topsicle_demo/result_all
telo_pattern=CCCTAAA

# main.py run

# Output: telolengths_all.csv for telomere in all reads 
# and additional: mean window plots, .csv files of raw count of each kmer pattern

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

# overview_plot.py run 
# output: heatmap of patterns and the matchs after it 
python3 overview_plot.py \
  --inputDir "$input_dir" \
  --outputDir "$output_dir" \
  --pattern $telo_pattern \
  --minSeqLength 9000 \
  --telophrase 4 \
  --recfindingpattern

# then, sbatch col_0_test.sh in the terminal
