# summarize results of LAMPSA

from __future__ import division
import os, glob

d = "aml"; chr = 12
qfile = "/work2/raissa/data/qlist_%s/qlist_outfile_dom_chr%d.txt"%(d,chr)
qlist = open(qfile).read().splitlines()
nfiles = len(qlist) # number of sets for corresponding chromosome

raw_p = []; adj_p = []; comb = []; cf_list = []; cf = []; file_idx = [];
ntargets = []; nrisk = []
sig_res = []; ctime = []

for i in range(nfiles):
	file = qlist[i]
	f = open(file).read().splitlines()
	nlines = len(f)
	if nlines > 0:
		for line_idx,line in enumerate(f):
			if "Correction factor" in line:
				cf_line = line.split(" ")
				cf_val = int(cf_line[7])
				cf.append(cf_val)
			if "Rank" in line: 
				res_line_idx = line_idx + 1
				hline_idx = line_idx
				header_line = f[hline_idx]+"\tCorrection Factor\tTotal Time\tSet#"
				break
			else: 
				res_line_idx = 0
		t = float(f[-1].split(" ")[-1])
		if res_line_idx != 0:
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
				file_idx.append(i+1)
				ctime.append(t)
		else:
			sig_res.append("NA\tNA\tNA\tNA\tNA\tNA\tNA\n")
			raw_p.append("NA")
			adj_p.append("NA")
			comb.append("NA")
			ntargets.append("NA")
			nrisk.append("NA")
			cf_list.append(cf_val)
			file_idx.append(i+1)
			ctime.append(t) 
		
		
# for i in range(nfiles):
# 	for j in range(3):
# 		file = "/work2/raissa/logrank_results/%s/%s_logrank_set%d_sub%d.txt"%(d,d,i+1,j+1)
# 		f = open(file).read().splitlines()
# 		nlines = len(f)
# 		for line_idx,line in enumerate(f):
# 			if "Correction factor" in line:
# 				cf_line = line.split(" ")
# 				cf_val = int(cf_line[7])
# 				cf.append(cf_val)
# 			if "Rank" in line: res_line_idx = line_idx + 1
# 		line_range = range(res_line_idx, nlines - 1)
# 		t = float(f[-1].split(" ")[-1])
# 		for line_idx in line_range:
# 			res_line = f[line_idx]
# 			sig_res.append(res_line)
# 			res = res_line.split("\t")
# 			raw_p.append(float(res[1]))
# 			adj_p.append(float(res[2]))
# 			comb.append(res[3].split(","))
# 			ntargets.append(int(res[5]))
# 			nrisk.append(int(res[6]))
# 			cf_list.append(cf_val)
# 			file_idx.append(i+1)
# 			ctime.append(t)

idx_sort = sorted(range(len(raw_p)), key=lambda k: raw_p[k])		
raw_p_sorted = [raw_p[i] for i in idx_sort]
adj_p_sorted = [adj_p[i] for i in idx_sort]
comb_sorted = [comb[i] for i in idx_sort]
ntargets_sorted = [ntargets[i] for i in idx_sort]
nrisk_sorted = [nrisk[i] for i in idx_sort]
cf_list_sorted = [cf_list[i] for i in idx_sort]
file_idx_sorted = [file_idx[i] for i in idx_sort]
sig_res_tmp = [sig_res[i] for i in idx_sort]
ctime_sorted = [ctime[i] for i in idx_sort]
sig_res_sorted = ["%s\t%s\t%.03f\t%d"%(e[0],e[1],e[2],e[3]) for e in zip(sig_res_tmp,cf_list_sorted,ctime_sorted,file_idx_sorted)]

# header_line = f[hline_idx]+"\tCorrection Factor\tTotal Time\tSet#"
res_file = [header_line] + sig_res_sorted
fo = open("/work2/raissa/logrank_results/%s_summary/%s_dom_chr%d.txt"%(d,d,chr), "w")
for line in res_file:
   fo.write("{}\n".format(line))

# cf_total = sum(cf)
# new_thres = thres/cf_total
# idx_filter = [i for i,e in enumerate(raw_p_sorted) if e<new_thres]
# sig_res_filtered = [sig_res_sorted[i] for i in idx_filter]			

# res_file_filtered = [header_line] + sig_res_filtered
# fo = open("/work2/raissa/logrank_results/%s_filtered_%.2g_%d_v3.txt"%(d,thres,cf_total), "w")
# for line in res_file_filtered:
#    fo.write("{}\n".format(line))

# fo = open("/work2/raissa/logrank_results/%s_wFilter_cf.txt"%(d), "w")
# for e in cf: fo.write("{}\n".format(e))




