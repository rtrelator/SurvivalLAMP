
Survival LAMP is an extended version of LAMP (Terada et al 2013) for performing log-rank test in survival analysis.

LAMP is a code for multiple testing correction to discover combinatorial effects. The lcm53.zip is the redistributed file that is provided in [http://research.nii.ac.jp/~uno/codes.htm.](http://research.nii.ac.jp/~uno/codes.htm)
The usage of LAMP is described in [http://a-terada.github.io/lamp/.](http://a-terada.github.io/lamp/)


## How to run LAMP for survival analysis:
1. Clone or download from [GitHub](https://github.com/rtrelator/SurvivalLAMP).
2. Compile: run `make` in corresponding directory
3. Run: 
```
python lampSA.py -p logrank [item_file] [status_file] [significance_level] -t [time_file] > [output_file]
```
EXAMPLE: 
```
python lampSA.py -p logrank sample/sample_logrank_item.csv sample/sample_logrank_status.csv 0.05 -t sample/sample_logrank_time.csv > sample/sample_logrank_output.txt
```


## Input Files:
1. ITEM FILE
- similar format to the original LAMP item file, this file contains association information of samples/individuals and markers
- first row as header (commented first entry, contains marker identifiers)
- first column contains the list of sample/individual IDs (one sample/individual per row)
- succeeding columns represent markers (e.g. gene, SNPs, etc.) (one marker per column)
- entries are 1 or 0 (integer type)
- file format: csv
- example: sample/sample_logrank_item.csv


2. STATUS FILE
- similar format to the original LAMP value file, except that it contains the status of samples/individuals
- first row as header (commented first entry similar to item file)
- first column must be the same as the item file 
- second column contains status of corresponding samples/individuals (1 = event, 0 = censored)
- file format: csv
- example: sample/sample_logrank_status.csv


3. TIME FILE
- contains survival time information of samples/individuals
- first row as header (commented first entry similar to item and status files)
- first column must be the same as the item and value files 
- second column contains survival time of corresponding samples/individuals (e.g. in months, years, or days)
- file format: csv
- example: sample/sample_logrank_time.csv


## Output File:
- similar format to the original LAMP output file
- example: sample/sample_logrank_output.txt
```
# Survival LAMP ver. 1.0
# item-file: ../LAMP_SurvivalAnalysis/lampSA_ver1/sample/sample_logrank_item.csv
# value-file: ../LAMP_SurvivalAnalysis/lampSA_ver1/sample/sample_logrank_status.csv
# time-file: ../LAMP_SurvivalAnalysis/lampSA_ver1/sample/sample_logrank_time.csv
# significance-level: 0.05
# P-value computing procedure: logrank
# # of tested elements: 9, # of samples: 292, # of positive samples: 101
# Adjusted significance level: 0.0005618, Correction factor: 89 (# of target rows >= 1)
# # of significant combinations: 5
Rank	Raw p-value	Adjusted p-value	Combination	Arity	# of target rows	# of failed targets
1	1.738e-05	0.0015468	r60_n9,G3PDH_570	2	14	10
2	2.4931e-05	0.0022188	r60_n9,G3PDH_570,r60_1,r60_3	4	9	7
3	8.6793e-05	0.0077246	r60_n9,G3PDH_570,r60_3	3	11	8
4	0.00039433	0.035096	r60_n9,r60_a22,G3PDH_570	3	10	7
5	0.00042816	0.038106	r60_n9,Pro25G,G3PDH_570	3	10	7
Time (sec.): Computing correction factor 0.743, P-value 1.304, Total 2.046
```

- *#* of tested elements = number of marker columns in item file
- *#* of samples = number of samples/individuals included (not counting censored samples/individuals after first failure time)
- *#* of positive samples = number of failed samples/individuals (i.e. status = 1)
- *#* of target rows = number of samples/individuals affected by the combination
- *#* of failed targets = number of samples/individuals affected by the combination whose status is 1 (i.e. failures)
