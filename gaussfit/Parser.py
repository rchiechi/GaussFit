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

import sys,os,logging,warnings,csv
from gaussfit.Parseopts import Opts,ShowUsage 
from gaussfit.colors import *

try:
	from scipy.optimize import curve_fit,OptimizeWarning
	from scipy.interpolate import UnivariateSpline
	import numpy as np
	# SciPy throws a useless warning for noisy J/V traces
	warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	ShowUsage()

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*',UserWarning)
#warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)

class Parse():
	'''
	This is the main parsing class that takes input data
	in the form of text files, parses them into dictionary
	and list attributes and provides methods for performing
	operations: Gaussian fits, F-N calculations, Rectification,
	Vtrans, and dJ/DV.
	'''
	def __init__(self,opts):
		self.opts = opts
		self.parsed = {}
		self.XY ={}
		self.X = np.array([])
		self.FN = {}
		self.compliance_traces = []
		self.ohmic = []
		self.DJDV = {}
		self.GHists = {}
		self.filtered = []
		self.R = {}
	def isfloat(self,f):
		try:
			float(f)
			return True
		except ValueError:
			return False
	def tofloat(self, f):
		try:
			f = float(f)
		except ValueError:
			f = np.NAN
		return f
		
	def ReadFiles(self, fns):
		''' Walk through input files and parse
		them into attributes '''
		Ycol = self.opts.Ycol
		Xcol = self.opts.Xcol
		line_idx = 0
		uniqueX = {}
		for fn in fns:
			logging.info("Parsing %s%s%s", TEAL,fn,YELLOW)
			try:
				rows = []
				with open(fn, 'rt', newline='') as csvfile:
					for row in csv.reader(csvfile, dialect='JV'):
						rows.append(row)
			except FileNotFoundError:
				logging.error("%s not found.", fn)
				continue
			except csv.Error:
				logging.error("Error parsing %s", fn)
				continue
			labels = []
			numcol = 0
			for label in rows[0]:
				if label:
					numcol +=1 # Only count non-empty columns
					if not self.isfloat(label):
						labels.append(label)
			if numcol < Ycol:
				logging.critical("Ycol > total number of columns! (%d > %d)", Ycol, numcol)
				logging.warn("%sNot importing:%s %s!", RED, YELLOW,fn)
				continue
			logging.debug("Found "+str(numcol)+" columns in "+fn)
			if 0 < len(labels) >= Ycol:
				logging.debug("Assuming first row is column labels")
				logging.info('Y column is "%s"' % labels[Ycol])
				del(rows[0])
			else:
				logging.debug("No labels found in first row.")
			for row in rows:
				if len(row) < numcol:
					logging.warning("Mismatch in %s, expecting %d colums, got %d", labels[Ycol], numcol, len(row))
					continue
				x, y = self.tofloat(row[Xcol]), self.tofloat(row[Ycol])
				if y == np.NAN:
					logging.warn("Treating invalid Y %f as NAN", y)
				if x == np.NAN:
					logging.warn("Treating invalid X %f as NAN", x)
				try:
					z = np.log( abs(y)/x**2 )
				except:
					z = np.NAN
					if x != 0.0: logging.error("Couldn't comput FN for %f, %f", x, y, exc_info=True)
				self.parsed[line_idx]=(x,y,z)
				if x in uniqueX.keys(): 
					uniqueX[x][line_idx]=(y,z)
				else: 
					uniqueX[x]={line_idx:(y,z)}
				line_idx += 1

		if not len(uniqueX):
			logging.error("No files parsed.")
			sys.exit() # Otherwise we loop over an empty uniqueX forever

		for x in uniqueX:
			# Gather each unique X (usually Voltage) value
			# and pull the Y (usually current-density) values
			# associated with it
			logging.debug('Pulling Y where X=%0.2f', x)
			y, fn = [], []
			for i in sorted(uniqueX[x].keys()): 
				y.append(uniqueX[x][i][0])
				fn.append(uniqueX[x][i][1])
			y = np.array(y)
			logy = np.log10(abs(y))
			self.XY[x] = { "Y":y, 
				   "LogY":logy, 
				   "hist":self.dohistogram(logy,"J"), 
				   "FN":np.array(fn) }

		self.X = np.array(sorted(self.XY.keys()))
		logging.info("Done parsing input data")
		print("* * * * * * Computing dY/dX  * * * * * * * *")
		self.DJDV, self.GHists, self.filtered = self.dodjdv() # This must come first for self.ohmic to be populated!
		print("* * * * * * Computing Vtrans * * * * * * * *")
		self.FN["neg"], self.FN["pos"] = self.findmin()
		print("* * * * * * Computing |R|  * * * * * * * * *")
		self.R = self.dorect()
		print("* * * * * * * * * * * * * * * * * * *")

	def dorect(self):
		''' 
		Divide each value of Y at +V by Y at -V
		and build a histogram of rectification, R
		WARNING: There is no way to ensure that
		each value of Y is being divided by the
		corresponding value from the same J/V trace
		because we cannot guarantee that all input
		files are formatted identically and have 
		identical, complete sweeps with the same
		spacing of V.
		'''
		R = {}
		for x in self.X:
			R[x] = {'r':np.array([]),'hist':{"bin":[],"freq":0,"mean":0,"std":0,"var":0,"bins":[],"fit":[]}}
			if x <= 0: continue
			elif -1*x not in self.X:
				logging.warn("(Rectification) Didn't find a negative voltage for %d.", x)
				continue
			ypos, yneg = self.XY[x]["Y"], self.XY[-1*x]["Y"]
			try:
				R[x]={'r':ypos/yneg,'hist':self.dohistogram(ypos/yneg,"R")}
			except ValueError:
				# TODO: This should never be allowed to happen by the input parser!
				logging.warn("(Rectification) Length of Y values differs for +/- %f (%s != %s).", x, len(ypos), len(yneg) )
			# NOTE: We do not filter for maxr here, rather in the Gaussian calcualtion
		R['X'] = np.array(sorted(R.keys()))
		return R
	
	def dodjdv(self):
		'''
		Fit a spline function to X/Y data and 
		compute dY/dX and normalize 
		'''
		vfilterneg,vfilterpos = np.array( [self.X.min()]) , np.array( [self.X.max()] )
		if self.opts.vcutoff > 0:
			vfilterpos = self.X[self.X >= self.opts.vcutoff] 
			vfilterneg = self.X[self.X <= -1*self.opts.vcutoff]
		spls = {}
		splhists = {}
		filtered = [('Potential', 'dY/dV', 'Y')]
		for x in np.linspace(self.X.min(), self.X.max(), 200): 
			spls[x] = []
			splhists[x] = {'spl':[],'hist':{}}
		i = -1
		while True:
			i += 1
			try:
				y,yf = [],[]
				for x in self.X: y.append(self.XY[x]['Y'][i])
				# k (smoothing)  must be 4 for 
				# the derivative to be cubic (k=3)
				spl = UnivariateSpline( self.X, y, k=4).derivative()
				maxY = 0
				for x in spls:
					if abs(spl(x)) > maxY:
						maxY = abs(spl(x))
				for x in spls: 
					spls[x].append(spl(x)/maxY)
					#splhists[x]['spl'].append(np.log10(abs(spl(x))))
				
				#d = np.array(spl(vfilter))
				X = np.linspace(self.X.min(), self.X.max(), 100)
				
				dd =  UnivariateSpline(X, spl(X), k=4).derivative()
				d = dd(vfilterpos) #Compute d2J/dV2
				d += -1*dd(vfilterneg) #Compute d2J/dV2
				if len(d[d < 0]): # Hackish because any() wasn't working
					# record in the index where dY/dX is < 0 at vcutoff
					self.ohmic.append(i)  
				else:
					for x in spls: 
						splhists[x]['spl'].append(np.log10(abs(spl(x))))
					for x in self.X:
						# filtered is a list containing only "clean" traces
						filtered.append( (x, spl(x), self.XY[x]['Y'][i]) )
			except IndexError:
				break
		logging.info("Non-tunneling traces: %s (out of %s)" % 
					( len(self.ohmic), len( self.XY[ list(self.XY.keys())[0]]['Y']) ) )
		for x in splhists:
			splhists[x]['hist'] = self.dohistogram(np.array(splhists[x]['spl']))
		return spls, splhists, filtered

	def getminroot(self, spl):
		'''
		Return the root of the first derivative of a spline function
		that results in the lowest value of ln(Y/X^2) vs. 1/X
		'''
		splvals = {}
		for r in spl.derivative().roots():
			splvals[float(spl(r))] = r
		if not splvals:
			return np.NAN
			#return False
		# Because UnivariateSpline crashes when we use 
		# 1/X, the plots are flipped so we take the max
		return splvals[np.nanmax(list(splvals.keys()))]

	def findmin(self):
		'''
		Find the toughs of ln(Y/X^2) vs. 1/X plots
		i.e., Vtrans, by either interpolating the data with
		a spline function and finding X where dY/dX = 0
		that gives the most negative value of Y (opts.smooth)
		or simply the most negative value of Y (! opts.smooth)
		'''
		neg_min_x, pos_min_x = [],[]
		i = -1
		tossed = 0
		if self.opts.skipohmic:
			# Vtrans has no physical meaning for curves with negative derivatives
			logging.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation." % 
					( len(self.ohmic), len( self.XY[ list(self.XY.keys())[0]]['Y']) ) )
		while True:
			i += 1
			if i in self.ohmic and self.opts.skipohmic:
				continue
			try:
				pos,neg = {},{}
				x_neg, y_neg, x_pos, y_pos = [],[],[],[]
				for x in self.X:
					# Without smoothing, we have to toss shorts or we get nonsense values
					if abs(self.XY[x]['Y'][i]).max() >= self.opts.compliance and not self.opts.smooth:
						tossed += 1
						continue
					y = self.XY[x]['FN'][i]
					if x < 0: 
						neg[y] = x
						x_neg.append(x)
						y_neg.append(y)
					if x > 0: 
						pos[y] = x
						x_pos.append(x)
						y_pos.append(y)

				if not len(neg.keys()) or not len(pos.keys()):
					logging.warn("Skipping empty column in FN calculation.")
					continue
				if self.opts.smooth:
					logging.debug("Using interpolation on FN")
					rootneg = self.getminroot(UnivariateSpline( x_neg, y_neg, k=4 ))
					if ~np.isnan(rootneg):
						neg_min_x.append(rootneg)
					rootpos = self.getminroot(UnivariateSpline( x_pos, y_pos, k=4 ))
					if ~np.isnan(rootpos):
						pos_min_x.append(rootpos)
					if np.NAN in (rootneg,rootpos):
						logging.warn("No minimum found in FN derivative (-):%s, (+):%s" % (rootneg, rootpos) )
				else:
					neg_min_x.append(neg[np.nanmin(list(neg.keys()))])
					pos_min_x.append(pos[np.nanmin(list(pos.keys()))])
			except IndexError:
				break
		
		if tossed: logging.warn("Tossed %d compliance traces during FN calculation.", tossed)
		neg_min_x = np.array(neg_min_x)
		pos_min_x = np.array(pos_min_x)
		return self.dohistogram(neg_min_x, "Vtrans(-)"), self.dohistogram(pos_min_x, "Vtrans(+)")

	def gauss(self, x, *p):
		'''
		Return a gaussian function
		'''
		A, mu, sigma = p
		return A*np.exp(-(x-mu)**2/(2.*sigma**2))
	
	def dohistogram(self, Y, label=""):
		'''
		Return a histogram of Y-values and a gaussian
		fit of the histogram, excluding values that
		exceed either the compliance limit (for current
		or current-density) or the ceiling for R. We 
		would like to include all data in the histogram,
		but outliers sometimes confuse the fitting
		routine, which defeats the purpose of machine-fitting
		'''
		j_compliance = np.nonzero(Y > self.opts.compliance) #BROKEN
		r_compliance = np.nonzero(abs(Y) > self.opts.maxr)
		if len(j_compliance[0]) and label == "J":
			logging.warn("Tossing %d data points > %0.1f for %s histogram!", len(j_compliance[0]), self.opts.compliance, label)
			Y = Y[np.nonzero(Y <= self.opts.compliance)]
		if len(r_compliance[0]) and label == "R":
			logging.warn("Tossing %d data points > %0.1f for %s histogram!", len(r_compliance[0]), self.opts.maxr, label)
			Y = Y[np.nonzero(abs(Y) <= self.opts.maxr)]
		logging.debug("%d points to consider.", len(Y))
		freq, bins = np.histogram(Y, bins=self.opts.bins, density=False)      
		p0 = [1., Y.mean(), Y.std()]
		bin_centers = (bins[:-1] + bins[1:])/2
		try:
			# NOTE: maxfev is hard-coded at 10000 because larger values are slow
			# and don't offer any advantage, but we still want to iterate to find the best fit
			coeff, var_matrix = curve_fit(self.gauss, bin_centers, freq, p0=p0, maxfev=10000)
			hist_fit = self.gauss(bin_centers, *coeff)
		except RuntimeError as msg:
			logging.warning("|%s| Fit did not converge (%s)", label, str(msg), exc_info=False)
			coeff = p0
			hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])
		return {"bin":bin_centers, "freq":freq, "mean":coeff[1], "std":coeff[2], "var":coeff[2], "bins":bins, "fit":hist_fit}

	def PrintFN(self):
		'''
		Print Vtrans values to the command line for convinience
		'''
		for key in ('pos', 'neg'):
			print("|Vtrans %s| mean: %0.4f Standard Deviation: %f" % (key, self.FN[key]['mean'], self.FN[key]['var']) )
		print("* * * * * * * * * * * * * * * * * * *")

