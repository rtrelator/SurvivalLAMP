
Survival LAMP is an extended version of LAMP (Terada et al 2013), a multiple-testing correction method,
for performing log-rank test in survival analysis.

LAMP is a code for multiple testing correction to discover combinatorial effects. The lcm53.zip is the redistributed file that is provided in http://research.nii.ac.jp/~uno/codes.htm. 
The usage of LAMP is described in http://a-terada.github.io/lamp/.

## Input Files:
1. ITEM FILE
- similar format to the original LAMP item file, this file contains association information of samples/individuals and markers
- first row as header (commented first entry, contains marker identifiers)
- first column contains the list of sample/individual IDs (one sample/individual per row)
- succeeding columns represent markers (e.g. gene, SNPs, etc.) (one marker per column)
- entries are 1 or 0 (integer type)
- file format: csv
example: 
sample/sample_logrank_item.csv


2. STATUS FILE
- similar format to the original LAMP value file, except that it contains the status of samples/individuals
- first row as header (commented first entry similar to item file)
- first column must be the same as the item file 
- second column contains status of corresponding samples/individuals (1 = event, 0 = censored)
- file format: csv
example:
sample/sample_logrank_status.csv


3. TIME FILE
- contains survival time information of samples/individuals
- first row as header (commented first entry similar to item and status files)
- first column must be the same as the item and value files 
- second column contains survival time of corresponding samples/individuals (e.g. in months, years, or days)
- file format: csv
example:
sample/sample_logrank_time.csv


## Output File:
- similar format to the original LAMP output file
example:
sample/sample_logrank_output.txt


## How to run LAMP for survival analysis:
compile: `make`
run: `python lampSA_ver1.py -p logrank [item_file] [status_file] [significance_level] -t [time_file] > [output_file]`

EXAMPLE: `python lampSA_ver1.py -p logrank sample/sample_logrank_item.csv sample/sample_logrank_status.csv 0.05 -t sample/sample_logrank_time.csv > sample/sample_logrank_output.txt`
