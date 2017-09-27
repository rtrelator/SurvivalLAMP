#!/usr/bin/env python

"""
Copyright (c) 2013, LAMP development team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP development team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Define fuctions that is used in multiple_test.py
# This source includes calculate P-value and MASL of the log-rank test.
# @author Terada, 16, Apr, 2013

from __future__ import division
from scipy import stats
import sys, os, operator, math
import functionsSuper as fs
import pvalTable

pardir = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
sys.path.append(pardir)
import readFileSA

##
# Define class
# This class calculate function f that define in paper. (C(n1, x)/C(n0+n1, x))
# transaction_list: list of transactions
##
class FunctionOfX(fs.FunctionsSuper):
	def __init__(self, transaction_list, row_size):
		fs.FunctionsSuper.__init__(self)
		self.__t_size = len(transaction_list) # all transaction size
		self.__f_size = self.sumValue(transaction_list) # transaction size which have flag = 1 (n1)
		self.__pvalTable = pvalTable.PvalTable( row_size ) # P-value table
		self.__chiTable = pvalTable.PvalTable( row_size ) # P-value table
		if self.__f_size == 0:
			sys.stdout.write("Error: There is no up-regulate gene.\n")
			sys.exit()
		# Check the transaction value.
		# If the value is not 1 or 0, raise error.
		# Because fisher's exact test does not handle numerical value.
		for t in transaction_list:
			if not (t.value == 1.0 or t.value == 0.0):
				sys.stderr.write("Error: \"" + t.name + "\" value is " + str(t.value)+".\n")
				sys.stderr.write("       But value is 1 or 0 if you test by fisher's exact test.\n")
				sys.exit()

		# check the support size.
		# If support size larger than half of all data size, raise error.
		# Because this version does not treat x > (n1+n0)/2.
		"""
		if self.__f_size > (self.__t_size/2):
			e_out = "The support size larger than half of all transaction size.\n"
			e_out = e_out + "                 This version does not treat this case."
			sys.exit()
		"""

	def getN1(self):
		return self.__f_size

	##
	# calclate MASL
	##
	def funcF(self, x):
		p1 = p2 = 1.0
		chi1 = chi2 = 0.0
		total_row1 = self.__f_size
		total = self.__t_size
		# when x < n_u
		if x < total_row1:
			ovalues = [[x, 0], [total_row1 - x, total - total_row1]]
			p1, chi1 = self.__probabilityTable( ovalues )
			ovalues = [[0, x], [total_row1, total - total_row1 - x]]
			p2, chi2 = self.__probabilityTable( ovalues )
		# when x >= n_u
		else:
			ovalues = [[total_row1, x-total_row1], [0, total - x]]
			p1, chi1 = self.__probabilityTable( ovalues )
			ovalues = [[0, x], [total_row1, total - total_row1 - x]]
			p2, chi2 = self.__probabilityTable( ovalues )
		if p1 < p2:
			return p1
		else:
			return p2
	
	# RRedit
	def funcF4logrank(self, x, transaction_list):
		p = 1.0
		transaction_list.sort( key=operator.attrgetter('failuretime') )
		failure_times, time_idx = getFailureTimes( transaction_list )
		for time_j in failure_times:
			subtransaction_list = getSubTrans( transaction_list, time_j )
			Yj = len(subtransaction_list) # all transaction size: t_size 
			dj = 0 # transaction size which have flag = 1 (n1): f_size and time = time_j
			for sub_t in subtransaction_list: 
				if sub_t.value == 1 and sub_t.failuretime == time_j:
					dj = dj + 1

			pj = 1.0
			if x < dj:
				pj = nCr(dj, x)/nCr(Yj, x)
			else:
				#pj = 1/nCr(Yj, x)
				pj = 1/nCr(Yj, dj)
			p = p*pj
			
# 			if x < dj:
# 				ovalues_j = [[x, 0], [dj - x, Yj - dj]]
# 				p1j, chi1j = probabilityTable4logrank( ovalues_j, Yj, dj )
# 				ovalues_j = [[0, x], [dj, Yj - dj - x]]
# 				p2j, chi2j = probabilityTable4logrank( ovalues_j, Yj, dj )
# 			else:				
# 				ovalues_j = [[dj, x - dj], [0, Yj - x]]
# 				p1j, chi1j = probabilityTable4logrank( ovalues_j, Yj, dj )
# 				ovalues_j = [[0, x], [dj, Yj - dj - x]]
# 				p2j, chi2j = probabilityTable4logrank( ovalues_j, Yj, dj )
# 			if p1j < p2j:
# 				p = p*p1j
# 			else:
# 				p = p*p2j
# 			if p1j < 0: sys.exit() #rrdebug
# 			if p2j < 0: sys.exit() #rrdebug
# 			testfile.write("p1j = %.3g, p2j = %.3g\n"%(p1j,p2j))

		#testfile.write("p = %.3g\n"%p)
		return p
			
	##
	# Calculate p-value by using fisher's exact test.
	# transaction_list: List of transactions
	# flag_itemset_id: Transactions which have items
	##
	def calPValue(self, transaction_list, flag_transactions_id):
		ovalues = self.contingencyTable( transaction_list, flag_transactions_id, self.__t_size, self.__f_size )
#		print ovalues
		total_col1 = self.__f_size
		total_row1 = sum( ovalues[0] )
		p = self.__pvalTable.getValue( total_row1, ovalues[0][0] )
		chi = self.__chiTable.getValue( total_row1, ovalues[0][0] )
		if p < 0: # calculate P-value and save to the table
			p, chi = self.__probabilityTable(ovalues)
			self.__pvalTable.putValue( total_row1, ovalues[0][0], p )
			self.__chiTable.putValue( total_row1, ovalues[0][0], chi )
		return p, chi
	
	def calPValue4logrank(self, transaction_list, flag_transactions_id):
		transaction_list.sort( key=operator.attrgetter('failuretime') )
		failure_times, time_idx = getFailureTimes( transaction_list )
		#treatment = []
		#for i in range(len(transaction_list)):
			#if i in flag_transactions_id: treatment.append(1)
			#else: treatment.append(0)
		chi2N = 0; chi2D = 0; n_risk = 0
		for time_j in failure_times:
			subtransaction_list, flag_subtransactions_id = getSubTransAndFlagSubTrans( transaction_list, flag_transactions_id, time_j )
			Yj = len( subtransaction_list )
			Y1j = len( flag_subtransactions_id ) #Y1j = total_row1 flagged 1 (nu)
			Y0j = Yj - Y1j
			# compute contingency table values for subtransaction list
			dj = 0; d1j = 0
			for sub_t in subtransaction_list:
				if sub_t.value == 1 and sub_t.failuretime == time_j:
					dj = dj + 1
			for id in flag_subtransactions_id:
				flag_sub_t = subtransaction_list[id]
				if flag_sub_t.value == 1 and flag_sub_t.failuretime == time_j:
					d1j = d1j + 1
			if Yj > 1: 
				chi2N = chi2N + d1j - (dj*Y1j/Yj)
				chi2D = chi2D + Y0j*Y1j*dj*(Yj - dj)/((Yj-1)*Yj**2)
			n_risk = n_risk + d1j
			# chi2D = chi2D + (dj*Y1j/Yj)
		chi2 = chi2N**2/chi2D
		p = chi2pval4logrank(chi2)
		ovalues = contingencyTable4logrank( transaction_list, flag_transactions_id, failure_times[0])
		total_col1 = ovalues[0][0] + ovalues[1][0] # d1
		total_row1 = sum( ovalues[0] )
		self.__pvalTable.putValue( total_row1, ovalues[0][0], p )
		self.__chiTable.putValue( total_row1, ovalues[0][0], chi2 )
		return p, chi2, n_risk
		
	def __calMeans(self, ovalues):
		total = self.__t_size
		total_col1 = self.__f_size # the number of all flag 1 transaction (n1)
		total_col2 = total - total_col1 # the number of all flag 0 transactio (n0)
		total_row1 = ovalues[0][0] + ovalues[0][1]
		total_row2 = ovalues[1][0] + ovalues[1][1]
		means = []; means.append([0]*2); means.append([0]*2)
		means[0][0] = float(total_row1 * total_col1) / total
		means[0][1] = float(total_row1 * total_col2) / total
		means[1][0] = float(total_row2 * total_col1) / total
		means[1][1] = float(total_row2 * total_col2) / total
		return means

	##
	# Calculate probability of occurrence probability about table.
	# a: Top left of table
	# b: Top right of table
	# n1: Sum of top and bottom left (a + c)
	# n0: Sum of top and bottom right (b + d)
	##
	def __probabilityTable(self, ovalues):
		means = self.__calMeans(ovalues) # calculate the exception value
#		print ovalues
#		print means

		# Yate continuity correction
		yate_corr = 0
		for i in means:
			for j in i:
				if j < 5:
					yate_corr = 0.5
					break

		chi = 0
		for i in xrange(0, len(ovalues)):
			row = ovalues[i]
			for j in xrange(0, len(row)):
				chi = chi + (abs(row[j] - means[i][j]) - yate_corr)**2/means[i][j]

		return self.__chi2pval( chi ), chi

	def __chi2pval(self, chi):
		if (chi == 0.0):
			return 1.0
		else: # dimension = 1
			return (self.stdNorDistribution(chi**0.5)) * 2.0


def maxLambda(transaction_list):
	# Count each item size
	item_sizes = {}
	for t in transaction_list:
#		print t.itemset
		for item in t.itemset:
#			print item
			# If item does not exist in item_size, then make mapping to 0
			if not item_sizes.has_key(item):
				item_sizes[item] = 0
			item_sizes[item] = item_sizes[item] + 1
	# Get max value in item_sizes
	max_value = 0
	for i in item_sizes.itervalues():
		if i > max_value:
			max_value = i
	return max_value

# RRedit: Find failure times: the times when at least one subject failed and the first transaction index for that time
def getFailureTimes(transaction_list):
	failure_times = []; time_idx = []
	for i in xrange( len(transaction_list) ):
		t = transaction_list[i]
		if ( t.value == 1 ) and not ( t.failuretime in failure_times ):
			failure_times.append( t.failuretime )
			time_idx.append( i )
	return failure_times, time_idx

# RRedit: CONFIRM THIS --> output as input of funcF
# Find sub transactions from jth failure time
def getSubTrans(transaction_list, time_j):
	subtransaction_list = []
	for i in xrange( len(transaction_list) ):
		t = transaction_list[i]
		if ( t.failuretime >= time_j ):
				subtransaction_list.append(t)
	return subtransaction_list

# Find sub transactions and flagged sub transactions from jth failure time
def getSubTransAndFlagSubTrans(transaction_list, flag_transactions_id, time_j):
	subtransaction_list = []; flag_subtransactions_id = []
	new_idx = 0
	for i in xrange( len(transaction_list) ):
		t = transaction_list[i]
		if ( t.failuretime >= time_j ):
				subtransaction_list.append(t)
				if i in flag_transactions_id: 
					flag_subtransactions_id.append( new_idx )
				new_idx = new_idx + 1
	return subtransaction_list, flag_subtransactions_id

# RRedit: tailored from original code
def calMeans4logrank(ovalues, Yj, dj):
	total = Yj
	total_col1 = dj # the number of all flag 1 transaction (n1)
	total_col2 = total - total_col1 # the number of all flag 0 transactio (n0)
	total_row1 = ovalues[0][0] + ovalues[0][1]
	total_row2 = ovalues[1][0] + ovalues[1][1]
	means = []; means.append([0]*2); means.append([0]*2)
	means[0][0] = float(total_row1 * total_col1) / total
	means[0][1] = float(total_row1 * total_col2) / total
	means[1][0] = float(total_row2 * total_col1) / total
	means[1][1] = float(total_row2 * total_col2) / total
	return means
				
def probabilityTable4logrank(ovalues, Yj, dj):
	means = calMeans4logrank(ovalues, Yj, dj) # calculate the exception value
#		print ovalues
#		print means

	# Yate continuity correction
	yate_corr = 0
	for i in means:
		for j in i:
			if j < 5:
				yate_corr = 0.5
				break

	chi = 0
	for i in xrange(0, len(ovalues)):
		row = ovalues[i]
		for j in xrange(0, len(row)):
			chi = chi + (abs(row[j] - means[i][j]) - yate_corr)**2/means[i][j]

	return chi2pval4logrank( chi ), chi

def chi2pval4logrank(chi):
	if (chi == 0.0):
		return 1.0
	else: # dimension = 1
		return (stdNorDistribution4logrank(chi**0.5)) * 2.0

def stdNorDistribution4logrank(x):
	range_20_1 = range(1, 21)
	range_20_1.reverse()
	pi2 = 0.398942280401432677940
	is_value = -1
	y = abs(x)
	c = y*y
	p = 0.0
	z = math.exp(-c*0.5)*pi2
	if (y < 2.5):
		for i in range_20_1:
			p = i*c/(i*2+1+is_value*p)
			is_value = -is_value
		p = 0.5-z*y/(1.0-p)
	else:
		for i in range_20_1:
			p = i/(y+p)
		p = z/(y+p)
#		p = 2 * p # double p-value because returens about two-sided test.
#		print str(x) + " " + str(p)
	return p

def contingencyTable4logrank( transaction_list, flag_transactions_id, failure_time_j):
	transaction_list.sort( key=operator.attrgetter('failuretime') )
	failure_times, time_idx = getFailureTimes( transaction_list )
	subtransaction_list = getSubTrans(transaction_list, failure_time_j)
	Yj = len(subtransaction_list); total = Yj
	dj = 0
	for sub_t in subtransaction_list:
		if sub_t.value == 1 and sub_t.failuretime == failure_time_j:
			dj = dj + 1
	total_col1 = dj
	ovalues = [ [0, 0], [0, 0] ]
	total_col2 = total - total_col1 # the number of all flag 0 transactio (n0)
	# count trahsaction which contains itemset and flag is 1. (This is indicate a of paper.)
	total_row1 = len(flag_transactions_id) # count all size that flag = 1 (x of paper)
	for i in flag_transactions_id:
		t = transaction_list[i]
		# If t flag = 1, then sum_has_flag ++.
		if t.value == 1 and t.failuretime == failure_time_j:
			ovalues[0][0] = ovalues[0][0] + t.value
	ovalues[0][0] = int(ovalues[0][0])
	ovalues[0][1] = total_row1 - ovalues[0][0] # the number of transaction which contains itemset and flag is 0 (This is indicate b of paper)
	ovalues[1][0] = total_col1 - ovalues[0][0]
	ovalues[1][1] = total_col2 - ovalues[0][1]
	return ovalues

def nCr(n,r):
    f = math.factorial
    return f(n) / (f(r)*f(n-r))

def run(xls_file, value_file, itemset_str_lst, delimiter):
	transaction_list, columnid2name, lcm2transaction_id = readFileSA.readFiles(xls_file, value_file, delimiter, time_file)
	transaction_list.sort( key=operator.attrgetter('failuretime') )
	max_lambda = maxLambda(transaction_list)
	func = FunctionOfX(transaction_list, max_lambda)
	colname2id_dict = readFileSA.colname2id(columnid2name)

	itemset = set()
	for i in itemset_str_lst:
		item_id = colname2id_dict[i]
		itemset.add(item_id + 1)

	flag_transactions_id = []
	for i in xrange( len(transaction_list) ):
		t = transaction_list[i]
		if len( itemset & t.itemset ) == len(itemset):
			flag_transactions_id.append( i )
	p_value, stat_value = func.calPValue(transaction_list, flag_transactions_id)
	n = len(transaction_list)
	n1 = func.getN1()
	sys.stdout.write("p-value: %s (N: %s, n1: %s, x: %s, chi: %s)\n"
					 % (p_value, n, n1, len(flag_transactions_id), stat_value))
	return p_value, len(flag_transactions_id)

if __name__ == "__main__":
	if (len(sys.argv) < 4):
		sys.stderr.write("Error: functions4logrank.py [item-file] [value_file] [itemset] [time_file]\n")
		sys.exit()

	xls_file = sys.argv[1]
	value_file = sys.argv[2]
	itemset_str_lst = sys.argv[3].split(',')
	delimiter = ','
	p_value, down_size = run(xls_file, value_file, itemset_str_lst, delimiter)
