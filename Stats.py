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
		SetA, SetB = GetStats(opts)
		Ttest(SetA,SetB,'R')
		Ttest(SetA,SetB,'J')
		TtestFN(SetA,SetB)
		WriteGNUplot(os.path.join(opts.out_dir,opts.outfile+"_Ttest"))


	elif len(opts.setA):
		SigTest(opts)
	else:
		logging.error("No input files to parse!")

def SigTest(opts):
	parser_full = Parse(opts)
	parser_full.ReadFiles(opts.setA)

	parsers = {}
	for n in range(1,len(opts.setA)):
		parsers[n] = Parse(opts)
		parsers[n].ReadFiles(opts.setA[0:n])
	pvalues = {}
	for x in parser_full.X:
		pvalues[x] = []
		for n in parsers:
			t,p= ttest_ind(parser_full.XY[x]['Y'], parsers[n].XY[x]['Y'])
			pvalues[x].append(p)
	n = 0
	for p in pvalues[list(pvalues.keys())[-1]]:
		n += 1
		print("%s : %s" % (str(n), str(p)))


def GetStats(opts):
	
	logging.info("Gathering Statistics")
	
	SetA = {}
	SetB = {}

	for f in opts.setA:
		parser = Parse(opts)
		parser.ReadFiles([f])
		for x in parser.X:
			if x not in SetA:
				SetA[x] = {'J':[],'R':[], 'FN':{'pos':[],'neg':[]}}
			SetA[x]['J'].append(parser.XY[x]['hist']['Gmean'])
			if x > 0:
				SetA[x]['R'].append(parser.R[x]['hist']['Gmean'])
			SetA[x]['FN']['neg'].append(parser.FN['neg']['Gmean'])
			SetA[x]['FN']['pos'].append(parser.FN['pos']['Gmean'])
		

	for f in opts.setB:
		parser = Parse(opts)
		parser.ReadFiles([f])
		for x in parser.X:
			if x not in SetB:
				SetB[x] = {'J':[],'R':[], 'FN':{'pos':[],'neg':[]}}
			SetB[x]['J'].append(parser.XY[x]['hist']['Gmean'])
			if x > 0:
				SetB[x]['R'].append(parser.R[x]['hist']['Gmean'])
			SetB[x]['FN']['neg'].append(parser.FN['neg']['Gmean'])
			SetB[x]['FN']['pos'].append(parser.FN['pos']['Gmean'])
		
	if SetA.keys()  != SetB.keys():
		logging.error("Are you trying to compare two datasets with different voltage steps?")
		sys.exit()

	return SetA, SetB
	
def Ttest(SetA, SetB, key):
	
	
	logging.info("Performing independent T-test on %s" % key)
	logging.info("Gathering mean %s-values" % key)
	
	mua = {}
	mub = {}

	dataset = ([],[],[])
	for x in SetA:
		t,p= ttest_ind(SetA[x][key], SetB[x][key], equal_var=False)
		if p < 0.01:
			c = GREEN
		else:
			c = RED
		print("P-value at %s%s V%s: %s%s%s" % (TEAL,str(x),RS,str(c),p,RS) )
		dataset[0].append(x)
		dataset[1].append(p)
		dataset[2].append(t)

	WriteGeneric(dataset, "Ttest"+key, opts, ["Voltage", "P-value", "T-stat"])

def TtestFN(SetA, SetB):
	x = list(SetA.keys())[-1]
	t_pos, p_pos = ttest_ind(SetA[x]['FN']['pos'], SetB[x]['FN']['pos'], equal_var=False)	
	t_neg, p_neg = ttest_ind(SetA[x]['FN']['neg'], SetB[x]['FN']['neg'], equal_var=False)
	
	if p_pos < 0.01:
		c = GREEN
	else:
		c = RED
	print("P-Value Vtrans(+): %s%s%s" % (c,str(p_pos),RS))
	
	if p_neg < 0.01:
		c = GREEN
	else:
		c = RED
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
Go(opts)
