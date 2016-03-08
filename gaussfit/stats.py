#!/usr/bin/env python3
'''
Copyright (C) 2015 Ryan Chiechi <r.c.chiechi@rug.nl>
Description:
        This program parses raw current-voltage data obtained from
        molecular tunneling junctions. It is specifically designed
        with EGaIn in mind, but may be applicable to other systems.

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


from gaussfit.Output import Writer,WriteStats
#from gaussfit.StatArgs import Opts
from gaussfit.Parser import Parse
#from gaussfit.Output import Writer

import sys,os,logging,warnings,csv
from gaussfit.colors import *

try:
	#from scipy.optimize import curve_fit,OptimizeWarning
	#import scipy.interpolate 
	from scipy.stats import gmean,kstest,ttest_ind,ttest_rel,ttest_1samp
	#import scipy.misc 
	import numpy as np
	np.seterr(invalid='raise')
	# SciPy throws a useless warning for noisy J/V traces
	#warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	sys.exit()
#warnings.filterwarnings('error','.*Mean of empty slice.*', RuntimeWarning)
#warnings.filterwarnings('error','.*Degrees of freedom <= 0 for slice.*', RuntimeWarning)
#warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*',UserWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)

class Stats:

	def __init__(self,opts):

		self.opts = opts
		self.SetA = []
		self.SetB = []
		self.PopA = []
		self.PopB = []
		self.dataset = {}
		self.fnstats = {}
		logging.info("Gathering Statistics")

		
		self.SetA, self.PopA = self.__getsetpop(self.opts, self.opts.setA)

		#if not len(self.opts.setB):
		#	return

		self.SetB, self.PopB = self.__getsetpop(self.opts, self.opts.setB)
		if self.SetA.keys()  != self.SetB.keys():
			logging.error("Are you trying to compare two datasets with different voltage steps?")
			print("Set A:")
			print("|",end='')
			for k in self.SetA.keys(): print(k, end='|')
			print("\nSet B:")
			print("|",end='')
			for k in self.SetB.keys(): print(k, end='|')
			print("\n")
			sys.exit()
		self.Ttest('R')
		self.Ttest('J')
		self.TtestFN()
		if self.opts.write:
			writer = Writer(self)
			writer.WriteGNUplot("statplot")
		
	def __getsetpop(self,opts,pop_files):
		Set, Pop = {}, {}
		maxfev = opts.maxfev
		opts.maxfev=1000
		parser = Parse(opts)
		#parser.ReadFiles(pop_files)
		#for x in parser.X:
		#	if x == 0:
		#		continue
		#	Pop[x] = {'J':{'mean':0,'std':0},'R':{'mean':0,'std':0}, 'FN':{'pos':{'mean':0,'std':0},'neg':{'mean':0,'std':0}}}
		#	Pop[x]['J']['mean'] = parser.XY[x]['hist']['mean']
		#	Pop[x]['J']['std'] = parser.XY[x]['hist']['std']
		#	Pop[x]['R']['mean'] = parser.R[abs(x)]['hist']['mean']
		#	Pop[x]['R']['std'] = parser.R[abs(x)]['hist']['std']
		#	Pop[x]['FN']['neg']['mean'] = parser.FN['neg']['mean']
		#	Pop[x]['FN']['neg']['std'] = parser.FN['neg']['std']
		#	Pop[x]['FN']['pos']['mean'] = parser.FN['pos']['mean']
		#	Pop[x]['FN']['pos']['std'] = parser.FN['pos']['std']
		traces = 0
		trace_error = 90
		for f in pop_files:
			opts.maxfev=maxfev
			parser = Parse(opts)
			try:
				parser.ReadFiles([f])
			except FloatingPointError as msg:
				logging.warning("Skipping %s because of %s." % (f, str(msg)))
				continue
			for x in parser.X:
				if x == 0:
					continue
				if x not in Set:
					Set[x] = {'J':{'mean':[],'std':[], 'Y':[]},'R':{'mean':[],'std':[], 'Y':[]}, \
							'FN':{'pos':{'mean':[],'std':[]},'neg':{'mean':[],'std':[]}}}
				Set[x]['J']['mean'].append(parser.XY[x]['hist']['Gmean'])
				Set[x]['J']['std'].append(parser.XY[x]['hist']['Gstd'])
				Set[x]['J']['Y'].append(parser.XY[x]['LogY'])
				Set[x]['R']['mean'].append(parser.R[abs(x)]['hist']['Gmean'])
				Set[x]['R']['std'].append(parser.R[abs(x)]['hist']['Gstd'])
				Set[x]['R']['Y'].append(parser.R[abs(x)]['r'])
				Set[x]['FN']['pos']['mean'].append(parser.FN['pos']['Gmean'])
				Set[x]['FN']['pos']['std'].append(parser.FN['pos']['Gstd'])
				Set[x]['FN']['neg']['mean'].append(parser.FN['neg']['Gmean'])
				Set[x]['FN']['neg']['std'].append(parser.FN['neg']['Gstd'])
				#print(parser.R[abs(x)]['r'])
				#if traces == 0:
				#	traces = len(parser.R[abs(x)]['r'])
				#elif len(parser.R[abs(x)]['r']) != traces:
				#	trace_error = (x, traces, len(parser.R[abs(x)]['r']))
			#		print(trace_error)
			#if trace_error:
			#	logging.warning("%s has a different numbers of traces (%s vs %s)." % trace_error)
			#	traces = 0
		return Set, Pop	

	def Ttest(self, key):
		
		logging.info("Performing independent T-test on %s" % key)
		logging.info("Gathering mean %s-values" % key)
		
		dataset = ([],[],[],[],[])
		
		p_vals = []
		aA_vals = []
		aB_vals = []

		for x in self.SetA:
			t,p= ttest_ind(self.SetA[x][key]['mean'], self.SetB[x][key]['mean'], equal_var=False)
			p_vals.append(p)
			
			#if p < 0.001: c = GREEN
			#else: c = RED
			#print("P-value at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,c,str(p),RS) )
			#print(self.SetA[x][key]['Y'])
			
			ssetsA = self.Sortsets(self.SetA[x][key]['Y'])
			ssetsB = self.Sortsets(self.SetB[x][key]['Y'])
			AlphaA = []
			AlphaB = []
			for k in ssetsA:
				A = self.Cronbach(ssetsA[k])
				aA_vals.append(A)
				AlphaA.append(A)
			for k in ssetsB:
				A = self.Cronbach(ssetsB[k])
				aB_vals.append(A)
				AlphaB.append(A)

			AlphaA = np.array(AlphaA)
			AlphaB = np.array(AlphaB)
			#AlphaA = self.Cronbach( self.SetA[x][key]['Y'] )
			#AlphaB = self.Cronbach( self.SetB[x][key]['Y'] )
			
			#aA_vals.append(AlphaA)
			#aB_vals.append(AlphaB)

			#if AlphaA > 0.7: c = GREEN
			#else: c = RED
			#print("α SetA at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,c,str(AlphaA),RS) )
			#if AlphaB > 0.7: c = GREEN
			#else: c = RED
			#print("α SetB at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,c,str(AlphaB),RS) )


			dataset[0].append(x)
			dataset[1].append(p)
			dataset[2].append(t)
			dataset[3].append(AlphaA[AlphaA > 0].mean())
			dataset[4].append(AlphaB[AlphaB > 0].mean())
			
			#dataset[3].append(AlphaA)
			#dataset[4].append(AlphaB)

		if self.opts.write:
			logging.info("Writing stats to %s" % self.opts.outfile+"_Ttest"+key+".txt")
			WriteStats(self.opts.out_dir, self.opts.outfile, dataset, \
				"Ttest"+key, ["Voltage", "P-value", "T-stat", "AlphaA", "AlphaB"])
		
		p_vals = np.array(p_vals)
		aA_vals = np.array(aA_vals)
		aB_vals = np.array(aB_vals)
		
		if p_vals.mean() < 0.001: c = GREEN
		else: c = RED
		print("p-value Mean: %s%s%s" % (GREEN, p_vals.mean(), RS))
		if aA_vals[aA_vals > 0].mean() > 0.7: c = GREEN
		else: c = RED
		print("α  SetA Mean: %s%s%s" % (c,aA_vals[aA_vals > 0].mean(),RS))
		if aB_vals[aB_vals > 0].mean() > 0.7: c = GREEN
		else: c = RED
		print("α  SetB Mean: %s%s%s" % (c,aB_vals[aB_vals > 0].mean(),RS))
		
		self.dataset[key] = dataset

	def TtestFN(self):
		x = list(self.SetA.keys())[-1]
		t_pos, p_pos = ttest_ind(self.SetA[x]['FN']['pos']['mean'], \
				self.SetB[x]['FN']['pos']['mean'], equal_var=False)	
		t_neg, p_neg = ttest_ind(self.SetA[x]['FN']['neg']['mean'], \
				self.SetB[x]['FN']['neg']['mean'], equal_var=False)
			
		if p_pos < 0.001: c = GREEN
		else: c = RED
		print("P-Value Vtrans(+): %s%s%s" % (c,str(p_pos),RS))
		
		if p_neg < 0.001: c = GREEN
		else: c = RED
		print("P-Value Vtrans(–): %s%s%s" % (c,str(p_neg),RS))
		self.fnstats = {"p_pos":p_pos, "p_neg":p_neg, "t_pos":t_pos, "t_neg":t_neg}

	def Sortsets(self, traces):
		ssets = {}
		for t in traces:
			if len(t) in ssets:
				ssets[len(t)].append(t)
			else:
				ssets[len(t)] = [t]
		return ssets

	def Cronbach(self, muA):
		'''
		Compute Cronbach's Alpha
		'''
		muI = np.asarray(muA)
		K = len(muI)
		try:
			tscores = muI.sum(axis=0)
			muvars = muI.var(axis=1,ddof=1)
			Alpha =  K/(K-1.) * (1 - muvars.sum() / tscores.var(ddof=1))
		except ZeroDivisionError:
			logging.error("Division by zero error computing Alpha!")
			Alpha = np.inf
		except ValueError as msg:
			logging.debug("%s while computing Cronbach's Alpha." % str(msg))
			Alpha = np.inf
		return Alpha
