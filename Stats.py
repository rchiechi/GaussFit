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


#from gaussfit.Output import Writer,Plotter
from gaussfit.StatArgs import Opts
from gaussfit.Parser import Parse
from gaussfit.Output import Writer

import sys,os,logging,warnings,csv
from gaussfit.colors import *

try:
	from scipy.optimize import curve_fit,OptimizeWarning
	import scipy.interpolate 
	from scipy.stats import gmean,kstest,ttest_ind,ttest_rel,ttest_1samp
	import scipy.misc 
	import numpy as np
	# SciPy throws a useless warning for noisy J/V traces
	warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	sys.exit()

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*',UserWarning)
warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)

def Go(opts):
	'''
	Call this function to execute the parsing engine
	i.e., as main()
	'''

	if len(opts.setA) and len(opts.setB):
		SetA, SetB, PopA, PopB = GetStats(opts)
		Ttest(SetA,SetB,PopA,PopB,'R')
		Ttest(SetA,SetB,PopA,PopB,'J')
		TtestFN(SetA,SetB,PopA,PopB)
		WriteGNUplot(os.path.join(opts.out_dir,opts.outfile+"_Ttest"))


	elif len(opts.setA):
		SigTest(opts)
	else:
		logging.error("No input files to parse!")

def GetStats(opts):
	
	logging.info("Gathering Statistics")
		
	def Getsetpop(opts, pop_files):
		Set, Pop = {}, {}
		maxfev = opts.maxfev
		opts.maxfev=1000
		parser = Parse(opts)
		parser.ReadFiles(pop_files)
		for x in parser.X:
			if x == 0:
				continue
			Pop[x] = {'J':{'mean':0,'std':0},'R':{'mean':0,'std':0}, 'FN':{'pos':{'mean':0,'std':0},'neg':{'mean':0,'std':0}}}
			Pop[x]['J']['mean'] = parser.XY[x]['hist']['mean']
			Pop[x]['J']['std'] = parser.XY[x]['hist']['std']
			Pop[x]['R']['mean'] = parser.R[abs(x)]['hist']['mean']
			Pop[x]['R']['std'] = parser.R[abs(x)]['hist']['std']
			Pop[x]['FN']['neg']['mean'] = parser.FN['neg']['mean']
			Pop[x]['FN']['neg']['std'] = parser.FN['neg']['std']
			Pop[x]['FN']['pos']['mean'] = parser.FN['pos']['mean']
			Pop[x]['FN']['pos']['std'] = parser.FN['pos']['std']
		for f in pop_files:
			opts.maxfev=maxfev
			parser = Parse(opts)
			parser.ReadFiles([f])
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
		return Set, Pop

	SetA, PopA = Getsetpop(opts, opts.setA)

	if not len(opts.setB):
		return SetA, SetB, PopA, PopB

	SetB, PopB = Getsetpop(opts, opts.setB)
		
	if SetA.keys()  != SetB.keys():
		logging.error("Are you trying to compare two datasets with different voltage steps?")
		sys.exit()
	return SetA, SetB, PopA, PopB
	

def Ttest(SetA, SetB, PopA, PopB, key):
	
	logging.info("Performing independent T-test on %s" % key)
	logging.info("Gathering mean %s-values" % key)
	
	dataset = ([],[],[],[],[])
	
	p_vals = []
	aA_vals = []
	aB_vals = []

	for x in SetA:
		t,p= ttest_ind(SetA[x][key]['mean'], SetB[x][key]['mean'], equal_var=False)
		p_vals.append(p)
		if p < 0.001: c = GREEN
		else: c = RED
		print("P-value at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,c,str(p),RS) )
	
		AlphaA = Cronbach( SetA[x][key]['Y'] )
		AlphaB = Cronbach( SetB[x][key]['Y'] )
		
		aA_vals.append(AlphaA)
		aB_vals.append(AlphaB)

		if AlphaA > 0.7: c = GREEN
		else: c = RED
		print("α SetA at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,c,str(AlphaA),RS) )
		if AlphaB > 0.7: c = GREEN
		else: c = RED
		print("α SetB at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,c,str(AlphaB),RS) )

		dataset[0].append(x)
		dataset[1].append(p)
		dataset[2].append(t)
		dataset[3].append(AlphaA)
		dataset[4].append(AlphaB)

	WriteGeneric(dataset, "Ttest"+key, opts, ["Voltage", "P-value", "T-stat", "AlphaA", "AlphaB"])
	
	p_vals = np.array(p_vals)
	aA_vals = np.array(aA_vals)
	aB_vals = np.array(aB_vals)

	print("p-value Mean: %s" % p_vals.mean())
	print("α  SetA Mean: %s" % aA_vals[aA_vals > 0].mean())
	print("α  SetB Mean: %s" % aB_vals[aB_vals > 0].mean())

def TtestFN(SetA, SetB, PopA, PopB):
	x = list(SetA.keys())[-1]
	t_pos, p_pos = ttest_ind(SetA[x]['FN']['pos']['mean'], SetB[x]['FN']['pos']['mean'], equal_var=False)	
	t_neg, p_neg = ttest_ind(SetA[x]['FN']['neg']['mean'], SetB[x]['FN']['neg']['mean'], equal_var=False)
		
	if p_pos < 0.001: c = GREEN
	else: c = RED
	print("P-Value Vtrans(+): %s%s%s" % (c,str(p_pos),RS))
	
	if p_neg < 0.001: c = GREEN
	else: c = RED
	print("P-Value Vtrans(–): %s%s%s" % (c,str(p_neg),RS))


def WriteGeneric(dataset, bfn, opts, labels=[]):
	'''Output for a generic set of data expecting an n-dimensional array'''

	if len(labels) and len(labels) != len(dataset):
		logging.error("Length of column labels does not match number of data columns for WriteGeneric!")
		return

	lencola = len(dataset[0])
	for d in dataset:
		if len(d) != lencola:
			logging.error("Length of columns differs for WriteGeneric!")
			return

	fn = os.path.join(opts.out_dir,opts.outfile+"_"+bfn+".txt")
	with open(fn, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, dialect='JV')
		if len(labels):
			writer.writerow(labels)

		for n in range(0, lencola):
			row = []
			for i in range(0,len(dataset)):
				row.append(dataset[i][n])
			writer.writerow(row)

def Cronbach(muA):
	'''
	Compute Cronbach's Alpha
	'''
	muI = np.asarray(muA)
	K = len(muI)
	tscores = muI.sum(axis=0)
	muvars = muI.var(axis=1,ddof=1)
	Alpha =  K/(K-1.) * (1 - muvars.sum() / tscores.var(ddof=1))
	return Alpha

def CronbachBroken(sigmaX, sigmaI, quiet=True):
	K = len(sigmaI)
	sI = np.array(sigmaI)**2
	sX = float( ((np.array(sigmaI)-sigmaX)**2).sum()/K-1 )
	if K < 2:
		Alpha = np.nan
	else:
		Alpha = float((K/(K-1))*(1-(sI.sum()/sX)))
	if not quiet:
		print("--- Cronbach\'s Alpha ---")
		print("       α=%s K=%s       " % (Alpha,K) )
		print("-------------------------")
	return Alpha

def WriteGNUplot(bfn):
	txt = '''
FONT="Arial,22"
DBLUE='#3300FF'
PURPLE='#FF00FF'
BLACK='#000000'

#set terminal pdfcairo color enhanced font FONT
set terminal wxt enhanced font FONT

#set output "Ttest.pdf"
set key outside horizontal font "Arial, 12"
set pointsize 1
set ylabel "P-value" offset 1,0
set xlabel "Potential (V)" offset 0,0.5
set format y "10^{%%T}"
set logscale y

f(x) = 0.01

plot "%sJ.txt" using 1:($1 == 0 ? '-':$2) title "J" with p lc rgb DBLUE lw 2 pt 5, "%sR.txt" using 1:($1 == 0 ? '-':$2) title "R" with p lc rgb PURPLE lw 2 pt 5, f(x) lc rgb BLACK dt "-" notitle

pause -1
exit''' % (bfn, bfn)
	fn = open('ttest_plot.gpin', 'wt')
	fn.write(txt)
	fn.close()




opts = Opts
if opts.GUI:
		from gui import statsfilebrowser
		gui = statsfilebrowser.ChooseFiles(opts,Go)
else:
		Go(opts)
