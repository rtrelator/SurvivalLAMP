# summarize results of LAMPSA

from __future__ import division
import os, glob

#nfiles = 71; d = "ov"
nfiles = 59; d = "brca"
flist = ["/work2/raissa/logrank_results/%s/%s_logrank_set%d_250.txt"%(d,d,i+1) for i in range(nfiles)]
thres = 0.05

raw_p = []; adj_p = []; comb = []; cf_list = []; cf = []
ntargets = []; nrisk = []
sig_res = []
for file in flist:
	f = open(file).read().splitlines()
	nlines = len(f)
	for line_idx,line in enumerate(f):
		if "Correction factor" in line:
			cf_line = line.split(" ")
			cf_val = int(cf_line[7])
			cf.append(cf_val)
		if "Rank" in line: res_line_idx = line_idx + 1
	line_range = range(res_line_idx, nlines - 1)
	for line_idx in line_range:
		res_line = f[line_idx]
		sig_res.append(res_line)
		res = res_line.split("\t")
		raw_p.append(float(res[1]))
		adj_p.append(float(res[2]))
		comb.append(res[3].split(","))
		ntargets.append(int(res[5]))
		nrisk.append(int(res[6]))
		cf_list.append(cf_val)

# idx_sort = sorted(range(len(raw_p)), key=lambda k: raw_p[k])		
# raw_p_sorted = [raw_p[i] for i in idx_sort]
# adj_p_sorted = [adj_p[i] for i in idx_sort]
# comb_sorted = [comb[i] for i in idx_sort]
# ntargets_sorted = [ntargets[i] for i in idx_sort]
# nrisk_sorted = [nrisk[i] for i in idx_sort]
# sig_res_sorted = [sig_res[i] for i in idx_sort]
# 
# header_line = f[res_line_idx-1]
# res_file = [header_line] + sig_res_sorted
# fo = open("/work2/raissa/logrank_results/%s_all_v3.txt"%d, "w")
# for line in res_file:
#    fo.write("{}\n".format(line))
# 
# cf_total = sum(cf)
# new_thres = thres/cf_total
# idx_filter = [i for i,e in enumerate(raw_p_sorted) if e<new_thres]
# sig_res_filtered = [sig_res_sorted[i] for i in idx_filter]			
# 
# res_file_filtered = [header_line] + sig_res_filtered
# fo = open("/work2/raissa/logrank_results/%s_filtered_%.2g_%d_v3.txt"%(d,thres,cf_total), "w")
# for line in res_file_filtered:
#    fo.write("{}\n".format(line))

fo = open("/work2/raissa/logrank_results/%s_cf.txt"%(d), "w")
for e in cf: fo.write("{}\n".format(e))




