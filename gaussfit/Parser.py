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

import sys,os,logging,warnings,csv,threading
from gaussfit.colors import *

try:
	import pandas as pd
	from pandas import Series, DataFrame
	from scipy.optimize import curve_fit,OptimizeWarning
	import scipy.interpolate 
	from scipy.stats import gmean,norm
	import scipy.misc 
	import numpy as np
	# SciPy throws a useless warning for noisy J/V traces
	warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/pandas/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	sys.exit()

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*',UserWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*invalid value encountered.*',RuntimeWarning)
warnings.filterwarnings('ignore','.*Mean of empty slice.*', RuntimeWarning)
#warnings.filterwarnings('error','.*Degrees of freedom <= 0 for slice.*', RuntimeWarning)
#warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)



class Parse():
	'''
	This is the main parsing class that takes input data
	in the form of text files, parses them into dictionary
	and list attributes and provides methods for performing
	operations: Gaussian fits, F-N calculations, Rectification,
	Vtrans, and dJ/DV.
	'''
	def __init__(self,opts,handler=None,lock=None):
		self.opts = opts
		self.df = DataFrame()
		self.XY = {}
		self.X = np.array([])
		self.FN = {}
		self.compliance_traces = []
		self.ohmic = []
		self.DJDV = {}
		self.GHists = {}
		self.filtered = []
		self.R = {}
		self.traces = {}
		if lock:
			self.lock = lock
		else:
			self.lock = threading.Lock()
		if not handler:
			self.loghandler = logging.StreamHandler()
			self.loghandler.setFormatter(logging.Formatter(\
				fmt=GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
		else:
			self.loghandler = handler
		self.logger = logging.getLogger('parser')
		self.logger.addHandler(self.loghandler)
		self.logger.setLevel(getattr(logging,self.opts.loglevel.upper()))

		if self.opts.nomin:
			self.logger.debug("Using interpolation on FN")

	@classmethod
	def doOutput(cls,writer):
		writer.logger.info("Writing files...")
		writer.WriteParseInfo()
		writer.WriteVtrans()
		writer.WriteGNUplot('Vtransplot')
		writer.WriteFN()
		writer.WriteGauss()
		writer.WriteGNUplot('JVplot')
		writer.WriteData()
		writer.WriteDJDV()
		writer.WriteFiltered()
		writer.WriteData(True)
		writer.WriteRData()
		writer.WriteGNUplot('Rplot')
		try:
			writer.WriteHistograms()
			writer.WriteGNUplot('JVhistplot')
		except IndexError:
			print("Error outputting histrograms")
		try:
			writer.WriteGHistogram()
		except IndexError:
			print("Error outputting Ghistrograms")
		try:
			writer.WriteGMatrix()
			writer.WriteGNUplot('Gplot', ['parula.pal']) 
		except IndexError:
			print("Error outputting GMatrix")
		writer.logger.info("Done!")

#	def isfloat(self,f):
#		try:
#			float(f)
#			return True
#		except ValueError as msg:
#			self.logger.debug("%s: %s" %(f, str(msg)) )
#			return False
#	def tofloat(self, f):
#		try:
#			f = float(f)
#		except ValueError:
#			f = np.NAN
#		return f
#	def cleanrow(self, inrow):
#		# Hack to get rid of emtpy columns
#		return inrow
#		outrow = []
#		for i in inrow:
#			outrow.append(list(filter(None,i)))
#		return outrow
	def ReadFiles(self, fns):
		''' Walk through input files and parse
		them into attributes '''
		if self.opts.Ycol > 0:
			self.logger.info("Parsing two columns of data.")
			self.df = pd.concat((pd.read_csv(f,sep=self.opts.delim,usecols=(self.opts.Xcol,self.opts.Ycol),names=('V','J'),header=0) for f in fns),ignore_index=True)
		else:
			self.logger.info("Parsing all columns of data.")
			self.df = pd.concat((pd.read_csv(f,sep=self.opts.delim,header=None,skiprows=1) for f in fns),ignore_index=True)
			X,Y = [],[]
			for row in self.df.iterrows():
				for y in row[1][1:]:
					X.append(row[1][0])
					Y.append(y)
			self.df = DataFrame({'V':X,'J':Y})
		self.df['FN'] = np.log(abs(self.df.J)/self.df.V**2)
		uniqueX = {}
		for v in self.df.V.unique(): uniqueX[v]={}
		for row in self.df.iterrows(): uniqueX[row[1].V][row[0]]=(row[1].J,row[1].FN)
		if not len(uniqueX):
			self.logger.error("No files parsed!?")
			sys.exit() # Otherwise we loop over an empty uniqueX forever
		self.findTraces()
		x = self.df.V.unique()
		lowy = uniqueX[np.nanmin(x[x>0])][ list(uniqueX[np.nanmin(x[x>0])].keys())[0]  ][0]/10 # typical for electrometer
		for x in self.df.V.unique():
			# Gather each unique X (usually Voltage) value
			# and pull the Y (usually current-density) values
			# associated with it
			self.logger.debug('Pulling Y where X=%0.2f', x)
			y, fn = [], []
			for i in sorted(uniqueX[x].keys()): 
				y.append(uniqueX[x][i][0])
				fn.append(uniqueX[x][i][1])
			y = np.array(y)
			ynz = [] # Make a special Y-array with no zeros for dealing with fits
			for _y in y:
				if np.isnan(_y) or _y == 0: 
					ynz.append(lowy)
					#print("%s --> %s" %(_y,lowy))
				else: 
					ynz.append(_y)
			logy =  np.nan_to_num(np.log10(abs(np.array(ynz))))
			self.XY[x] = { "Y":y, 
				   "LogY":logy, 
				   "hist":self.dohistogram(logy,"J"), 
				   "FN":np.array(fn) }
		self.X = np.array(sorted(self.XY.keys()))
		self.logger.info("Done parsing input data")
		self.logger.info("* * * * * * Computing dY/dX  * * * * * * * *")
		self.DJDV, self.GHists, self.filtered = self.dodjdv() # This must come first for self.ohmic to be populated!
		self.logger.info("* * * * * * Computing Vtrans * * * * * * * *")
		self.FN["neg"], self.FN["pos"] = self.findmin()
		self.logger.info("* * * * * * Computing |R|  * * * * * * * * *")
		self.R = self.dorect()
		self.logger.info("* * * * * * * * * * * * * * * * * * *")
		self.PrintFN()

	def findTraces(self):
	
		def __checktraces(traces):
			self.logger.debug("Checking V starting from slice %s:%s" % (traces[0][0],traces[0][1]) )
			#V = self.df.V[traces[0][0]]
			lt = len(self.df.V[ traces[0][0]:traces[0][1] ])
			for trace in traces:
				if lt != len(self.df.V[trace[0]:trace[1]]):
					self.logger.warn("Unequal voltage steps in dataset!")
					return False
				self.logger.debug("Trace: %s -> %s" % (self.df.V[trace[0]],self.df.V[trace[-1]]) )
			self.logger.info("Traces look good.")
			return True
		
		traces = []
		ntraces = 0
		if self.df.V.value_counts().index[0] == 0.0:
			try:
				ntraces = int(self.df.V.value_counts()[0]/3) # Three zeros in every trace!
				for t in zip(*(iter(self.df[self.df.V == 0.00].V.index),) * 3):
					traces.append( (t[0],t[2]) )
			except ValueError:
				self.logger.warn("Did not find three-zero (EGaIn) traces!")
		if not ntraces:
			self.logger.warn("This does not look like an EGaIn dataset.")
			try:
				#ntraces = int(self.df.V.value_counts()[0]/2) # Two end-pointss in every trace!
				ntraces = int(self.df.V.value_counts()[self.df.V[0]]/2) # Two end-pointss in every trace!
				for t in zip(*(iter(self.df[self.df.V == self.df.V[0]].V.index),) * 2):
					traces.append( (t[0],t[1]) )
				#for t in zip(*(iter(self.df[self.df.V == self.df.V.value_counts().index[0]].V.index),) * 2):
				#	traces.append( (t[0],t[1]) )
			except ValueError:
				self.logger.warn("Did not find three-zero (EGaIn) traces!")

		if not ntraces or not __checktraces(traces):
			self.logger.warn("Recomputing traces based on repeat values")
			traces = []
			Vinit = self.df.V[0]
			trace = [0]
			for row in self.df[trace[0]:].iterrows():
				if row[1].V == Vinit:
					if len(trace) == 1:
						trace.append(row[0])
					elif len(trace) == 2:
						traces.append(trace)
						#print(self.df[trace[0]:trace[1]+1])
						trace = [row[0]]
			ntraces = len(traces)
		if not __checktraces(traces):
			self.logger.error("Problem with traces: FN and derivative probably will not work correctly!")
		self.logger.info("Found %s traces (%s)." % (ntraces,len(traces)) )
		self.traces = traces

	def dorect(self):
		''' 
		Divide each value of Y at +V by Y at -V
		and build a histogram of rectification, R
		'''
		R = {}
		r = {}
		for trace in self.traces:
			rows = {}
			for row in self.df[trace[0]:trace[1]+1].iterrows():
				if row[1].V in rows:
					#self.logger.warn("Traces are out of sync, do not trust R values!")
					#TODO Don't just average out hysterysis	
					rows[row[1].V] = (rows[row[1].V]+row[1].J)/2
				else:
					rows[row[1].V] = row[1].J
				if row[1].V not in r:
					r[row[1].V] = []
			for x in rows:
				#TODO AFM data break this!
				if x not in rows or -1*x not in rows:
					self.logger.warn("Rectification data missing voltages.")
					continue
				if x == 0.0:
					r[x].append(1.)
					continue
				if self.opts.logr:
					r[x].append(np.log10(abs(rows[x]/rows[-1*x])))
				else:
					r[x].append(abs(rows[x]/rows[-1*x]))
		
#		for trace in self.traces['trace']:
#			xy = {}
#			#print("X-range: %d, %d" % (self.traces['XY'][0][trace[0]],self.traces['XY'][0][trace[1]]) )
#			for i in range(trace[0],trace[1]+1):
#				x = self.traces['XY'][0][i]
#				xy[x] = self.traces['XY'][1][i]
#				if x not in r:
#					r[x] = []
#			for x in xy:
#				if x <= 0: continue
#				if -1*x not in xy:
#					self.logger.warn("%f not found in trace while computing R", -1*x)
#					continue
#				if self.opts.logr:
#					r[x].append( np.log10(abs(xy[x]/xy[-1*x])) )
#				else:
#					r[x].append(abs(xy[x]/xy[-1*x]))
		for x in r:
			if x < 0: continue
			R[x]={'r':np.array(r[x]),'hist':self.dohistogram(np.array(r[x]),"R")}
			R[-1*x] = R[x]
		R['X'] = np.array(sorted(R.keys()))
		return R

	def dodjdv(self):
		'''
		Fit a spline function to X/Y data and 
		compute dY/dX and normalize 
		'''
		if self.opts.vcutoff > 0:
			vfilterneg,vfilterpos = np.linspace(-1*self.opts.vcutoff,0,200), np.linspace(0,self.opts.vcutoff.max(),200)
		else:
			vfilterneg,vfilterpos = np.linspace(self.df.V.min(),0,200), np.linspace(0,self.df.V.max(),200)
		spls = {}
		splhists = {}
		filtered = [('Potential', 'Fit', 'Y')]
		linx = np.linspace(self.df.V.min(), self.df.V.max(), 200)
		for x in linx: 
			spls[x] = []
			splhists[x] = {'spl':[],'hist':{}}
		if self.opts.vcutoff > 0:
			vfilterneg,vfilterpos = linx[-1*self.opts.vcutoff < linx < 0],linux[0 < linx < self.opts.vcutoff]
		else:
			vfilterneg,vfilterpos = linx[linx < 0], linx[linx > 0]
		for col in range(0,len(self.traces)):
			#TODO Yuck!
			V,J = [],[]
			avg = {}
			fbtrace = self.df[self.traces[col][0]:self.traces[col][1]+1]
			for row in fbtrace.iterrows():
				if row[1].V in avg: avg[row[1].V].append(row[1].J)
				else: avg[row[1].V] = [row[1].J]
			for x in sorted(avg.keys()): 
				V.append(x)
				J.append(np.mean(avg[x]))
			try:
				spl = scipy.interpolate.UnivariateSpline(V,J, k=5, s=self.opts.smooth )
				dd =  scipy.interpolate.UnivariateSpline(V, J, k=5, s=None).derivative(2)
			except Exception as msg:
				continue
				self.logger.error('Error in derivative calulation: %s' % str(msg))

			spldd = dd(vfilterpos) #Compute d2J/dV2
			spldd += -1*dd(vfilterneg) #Compute d2J/dV2
			if len(spldd[spldd<0]):
					# record in the index where dY/dX is < 0 within vcutoff range
					self.ohmic.append(col)  
					if self.opts.skipohmic:
						continue
			else:
				for row in fbtrace.iterrows():
					# filtered is a list containing only "clean" traces			
					filtered.append( (row[1].V, spl(row[1].V), row[1].J) )
			err = None
			for x in sorted(spls.keys()):
				try:
					d = spl.derivatives(x)
				except ValueError as msg:
					err = str(msg)
					#self.logger.error('Error computing derivative: %s' % str(msg))
					continue
				if np.isnan(d[self.opts.heatmapd]):
					self.logger.warn("Got NaN computing dJ/dV")
					continue
				spls[x].append(d[self.opts.heatmapd])
				splhists[x]['spl'].append(np.log10(abs(d[self.opts.heatmapd])))
			if err:
				self.logger.error("Error while computing derivative: %s" % str(err))

		self.logger.info("Non-tunneling traces: %s (out of %0d)" % 
					( len(self.ohmic), len(self.traces) ) )
		for x in splhists:
			splhists[x]['hist'] = self.dohistogram(np.array(splhists[x]['spl']), label='DJDV')
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
		Find the troughs of ln(Y/X^2) vs. 1/X plots
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
			self.logger.info("Skipping %s (out of %s) non-tunneling traces for Vtrans calculation." % 
					( len(self.ohmic), (len(self.traces))))
			#len( self.XY[ list(self.XY.keys())[0]]['Y']) ) )
		for col in range(0,len(self.traces)):
			if self.opts.skipohmic and col in self.ohmic:
				continue
			fbtrace = self.df[self.traces[col][0]:self.traces[col][1]+1]
			avg = {}
			for row in fbtrace.iterrows():
				# Without smoothing, we have to toss shorts or we get nonsense values
				if abs(row[1].J).max() >= self.opts.compliance and not self.opts.nomin:
					tossed += 1
					continue
				if row[1].V in avg:
					avg[row[1].V].append(row[1].FN)
				else:
					avg[row[1].V] = [row[1].FN]
			Vpos,FNpos = [],[]
			Vneg,FNneg = [],[]
			for x in sorted(avg.keys()):
				if x > 0:
					Vpos.append(x)
					FNpos.append(np.mean(avg[x]))
				elif x < 0:
					Vneg.append(x)
					FNneg.append(np.mean(avg[x]))
			if self.opts.nomin and len(Vpos) and len(FNpos):
				try:
					rootpos = self.getminroot(scipy.interpolate.UnivariateSpline(Vpos,FNpos, k=4, s=None))
					rootneg = self.getminroot(scipy.interpolate.UnivariateSpline(Vneg,FNneg, k=4, s=None))
				except Exception:
					continue
				if np.isfinite(rootneg):
					neg_min_x.append(rootneg)
				if np.isfinite(rootpos):
					pos_min_x.append(rootpos)
				if np.NAN in (rootneg,rootpos):
					self.logger.warn("No minimum found in FN derivative (-):%s, (+):%s" % (rootneg, rootpos) )
			elif len(Vpos) and len(FNpos):
				neg_min_x.append( fbtrace.V[ fbtrace.FN[fbtrace.V < 0].idxmin() ] )
				pos_min_x.append( fbtrace.V[ fbtrace.FN[fbtrace.V > 0].idxmin() ] )

		if tossed: self.logger.warn("Tossed %d compliance traces during FN calculation.", tossed)
		neg_min_x = np.array(neg_min_x)
		pos_min_x = np.array(pos_min_x)
		return self.dohistogram(neg_min_x, "Vtrans(-)"), self.dohistogram(pos_min_x, "Vtrans(+)")

	def gauss(self, x, *p):
		'''
		Return a gaussian function
		'''
		A, mu, sigma = p
		return A*np.exp(-(x-mu)**2/(2.*sigma**2))
	
	def lorenz(self, x, *p):
		'''
		Return a lorenzian function
		'''
		A, mu, B = p
		return A/( (x-mu)**2 + B**2)

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
			self.logger.warn("Tossing %d data points > %0.1f for %s histogram!", len(j_compliance[0]), self.opts.compliance, label)
			Y = Y[np.nonzero(Y <= self.opts.compliance)]
		if len(r_compliance[0]) and label == "R":
			self.logger.warn("Tossing %d data points > %0.1f for %s histogram!", len(r_compliance[0]), self.opts.maxr, label)
			Y = Y[np.nonzero(abs(Y) <= self.opts.maxr)]
		if len(Y) < 10:
			self.logger.debug("< 10 points! %d points to consider.", len(Y))
		if label == "J":
			yrange = (Y.min()-1,Y.max()+1)
		else:
			yrange = None
		if label=='DJDV':
			nbins = self.opts.heatmapbins
		else:
			nbins = self.opts.bins
		
		try:
			freq, bins = np.histogram(Y, range=yrange, bins=nbins, density=False)
		except ValueError as msg:
			self.logger.warning("Encountered this error while constructing histogram: %s", str(msg), exc_info=False)
			bins=np.array([0.,0.,0.,0.])
			freq=np.array([0.,0.,0.,0.])
		Ym,Ys = 0.0,0.0
		if len(Y):	
			# Compute the geometric mean and give it the
			# correct sign
			if len(Y[Y<0]):
				Ym = -1*gmean(abs(Y))
			else:
				Ym = gmean(abs(Y))
			
			# Somehow this produces negative sigmas
			Ys = abs(Y.std())
		# This will trigger if the mean happens to be zero
		#if not Ys or not Ym:
		#	self.logger.error("Failed to find G-mean and G-std!")
		#	print(Ys,Ym)
		# Initital conditions for Gauss-fit
		p0 = [1., Ym, Ys]
		bin_centers = (bins[:-1] + bins[1:])/2
		try:
			with self.lock:
				if self.opts.lorenzian:
					coeff, covar = curve_fit(self.lorenz, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev)
					#hist_fit = self.lorenz(bin_centers, *coeff)
				else:
					coeff, covar = curve_fit(self.gauss, bin_centers, freq, p0=p0, maxfev=self.opts.maxfev)
					#hist_fit = self.gauss(bin_centers, *coeff)
				hist_fit = self.gauss(bin_centers, *coeff)
		except RuntimeError as msg:
			if self.opts.maxfev > 10:
				self.logger.warning("|%s| Fit did not converge", label, exc_info=False)
			coeff = p0
			hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])
		except ValueError as msg:
			self.logger.warning("|%s| Skipping data with ridiculous numbers in it (%s)", label, str(msg), exc_info=False )
			coeff=p0
			hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])

		return {"bin":bin_centers, "freq":freq, "mean":coeff[1], "std":coeff[2], \
				"var":coeff[2], "bins":bins, "fit":hist_fit, "Gmean":Ym, "Gstd":Ys}

	def PrintFN(self):
		'''
		Print Vtrans values to the command line for convinience
		'''
		for key in ('pos', 'neg'):
			self.logger.info("|Vtrans %s| Gauss-mean: %0.4f Standard Deviation: %f" % (key, self.FN[key]['mean'], self.FN[key]['std']) )
			self.logger.info("|Vtrans %s| Geometric-mean: %0.4f Standard Deviation: %f" % (key, self.FN[key]['Gmean'], self.FN[key]['Gstd']) )
		#print("* * * * * * * * * * * * * * * * * * *")

