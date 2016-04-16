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
from gaussfit.Parser import Parse

import sys,os,logging,warnings,csv,datetime,threading,time,math,tempfile
from gaussfit.colors import *

try:
	from scipy.stats import gmean,kstest,ttest_ind,ttest_rel,ttest_1samp
	import numpy as np
	np.seterr(invalid='raise')

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

class ParserThread(threading.Thread):

	def __init__(self,limiter,alive,lock,opts,files,Set,Pop,handler):
		threading.Thread.__init__(self)
		self.alive = alive
		self.files = files
		self.Set, self.Pop = Set,Pop
		self.lock = lock
		self.limiter = limiter
		self.parser = Parse(opts,handler,lock)
		self.logger = logging.getLogger('thread')
		self.logger.setLevel(getattr(logging,opts.loglevel.upper()))
		self.logger.addHandler(handler)
	def run(self):
		self.limiter.acquire()
		self.logger.info("Starting thread...")
		try:
			self.parser.ReadFiles([self.files])
			for x in self.parser.X:
				if not self.alive.isSet():
					break
				if x == 0:
					continue
				with self.lock:
					Stats.DoSet(self.parser, self.Set, x)
		except FloatingPointError as msg:
			self.logger.warning("Skipping %s because of %s." % (f, str(msg)))
		finally:
			self.limiter.release()
		
class Stats:

	def __init__(self,opts,handler=None):
		self.opts = opts
		self.SetA = {}
		self.SetB = {}
		self.PopA = {}
		self.PopB = {}
		self.dataset = {}
		self.fnstats = {}
		self.extrainfo=''
		if not handler:
			self.loghandler = logging.StreamHandler()
			self.loghandler.setFormatter(logging.Formatter(\
				fmt=GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
		else:
			self.loghandler = handler

		self.logger = logging.getLogger('statparser')
		self.logger.addHandler(self.loghandler)
		self.logger.setLevel(getattr(logging,self.opts.loglevel.upper()))
		self.logger.info("Gathering Statistics")
			
	def parse(self):
		if not self.opts.autonobs:
			self.SetA, self.PopA = self.__getsetpop(self.opts, self.opts.setA)
			self.SetB, self.PopB = self.__getsetpop(self.opts, self.opts.setB)
		else:
			self.SetA, self.PopA = self.__autonobs(self.opts, self.opts.setA)
			self.SetB, self.PopB = self.__autonobs(self.opts, self.opts.setB)

		if self.SetA.keys()  != self.SetB.keys():
			self.logger.error("Are you trying to compare two datasets with different voltage steps?")
			print("Set A:")
			print("|",end='')
			for k in self.SetA.keys(): print(k, end='|')
			print("\nSet B:")
			print("|",end='')
			for k in self.SetB.keys(): print(k, end='|')
			print("\n")
			sys.exit()
		
		if not self.SetA.keys() or not self.SetB.keys():
			print("Didn't parse any input files?")
			return
		
		self.Ttest('R')
		self.Ttest('J')
		self.TtestFN()

		if self.opts.write:
			writer = Writer(self)
			writer.WriteGNUplot("statplot")
			writer.WriteParseInfo(self.extrainfo)


	def __getsetpop(self,opts,pop_files):
		Set, Pop = {}, {}
		self.logger.info("Using %s threads." % opts.threads)
		alive = threading.Event()
		lock = threading.Lock()
		limiter = threading.BoundedSemaphore(opts.threads)
		alive.set()
		children = []
		for f in pop_files:
			children.append(ParserThread(limiter,alive,lock,opts,f,Set,Pop,self.loghandler))
			children[-1].start()
		for c in children:
			c.join()
		return Set,Pop

	@classmethod
	def DoSet(cls,parser,Set,x):
		if x not in Set:
			Set[x] = {'J':{'mean':[],'std':[], 'Y':[]},'R':{'mean':[],'std':[], 'Y':[]}, \
					'FN':{'pos':{'mean':[],'std':[]},'neg':{'mean':[],'std':[]}}}
		if abs(x) not in parser.R:
			parser.R[abs(x)] = {'hist':{'Gmean':1.0, 'Gstd':0.0},'r':[]}
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

	@classmethod
	def DoPop(cls,parser,Set,x):
		Pop[x] = {'J':{'mean':0,'std':0},'R':{'mean':0,'std':0}, 'FN':{'pos':{'mean':0,'std':0},'neg':{'mean':0,'std':0}}}
		Pop[x]['J']['mean'] = parser.XY[x]['hist']['mean']
		Pop[x]['J']['std'] = parser.XY[x]['hist']['std']
		Pop[x]['R']['mean'] = parser.R[abs(x)]['hist']['mean']
		Pop[x]['R']['std'] = parser.R[abs(x)]['hist']['std']
		Pop[x]['FN']['neg']['mean'] = parser.FN['neg']['mean']
		Pop[x]['FN']['neg']['std'] = parser.FN['neg']['std']
		Pop[x]['FN']['pos']['mean'] = parser.FN['pos']['mean']
		Pop[x]['FN']['pos']['std'] = parser.FN['pos']['std']

	def __autonobs(self,opts,pop_files):
		Set,Pop = {},{}
		cTimes = {}
		for f in pop_files:
			try:
				with open(f.replace('_data.txt','')) as fh:
					lines = fh.readlines()
					dt = datetime.datetime.strptime('%s; %s' % (lines[0].strip(), lines[1].strip()),\
							"%A, %B %d, %Y; %I:%M %p")
					cTimes[f] = dt
			except FileNotFoundError:
				self.logger.debug("%s does not exist, not doing autonobs." % f.replace('_data.txt',''))
				return Set, Pop	
		diffdays = {}
		diffmins = {}
		min_factor = 60
		for f in cTimes:
			for s in cTimes:
				td = cTimes[s]-cTimes[f]
				if td.total_seconds() == 0:
					continue
				dd = abs(td.days)
				if dd not in diffdays:
					diffdays[dd] = [(f,s)]
				elif (s,f) not in diffdays[dd] and (f,s) not in diffdays[dd]:
					diffdays[dd].append((f,s))
				dm = math.floor(abs(td.total_seconds())/(min_factor*60))
				if dm not in diffmins:
					diffmins[dm] = [(f,s)]
				elif (s,f) not in diffmins[dm] and (f,s) not in diffmins[dm]:
					diffmins[dm].append((f,s))
			
		def __maketmpfiles(popfiles):
				outfile = tempfile.NamedTemporaryFile(mode='w+t',delete=False)
				header = False
				for f in popfiles:
					with open(f, 'rt', newline='') as csvfile:
						writer = csv.writer(outfile,dialect='JV')
						for row in csv.reader(csvfile,dialect='JV'):
							try:
								float(row[0])
							except ValueError:
								if not header:
									header = True
									writer.writerow(row)
								continue
							writer.writerow(row)
				#tmpfiles[outfile.name]=outfile
				return outfile.name
		
		tmpfiles = []
		if len(diffdays.keys()) < 3:
			days, minutes = False, True
		elif len(diffmins.keys()) < 100:
			days, minutes = True, False
		else:
			self.logger.error("Too few days apart and too many minutes!")
			sys.exit()
		if days:
			for d in diffdays:
				popfiles = []
				for f in diffdays[d]:
					if f[0] not in popfiles:
						popfiles.append(f[0])
					if f[1] not in popfiles:
						popfiles.append(f[1])
				n = __maketmpfiles(popfiles)
				tmpfiles.append(n)
				self.logger.info("Measurements done %s days apart: %s " % (d,n))
				self.extrainfo += "Measurements done %s days apart: %s\n" % (d,n)	
				for p in popfiles:
					self.extrainfo += p+"\n"
		if minutes:
			for m in diffmins:
				popfiles = []
				for f in diffmins[m]:
					if f[0] not in popfiles:
						popfiles.append(f[0])
					if f[1] not in popfiles:
						popfiles.append(f[1])
				n = __maketmpfiles(popfiles)
				tmpfiles.append(n)
				self.logger.info("Measurements done %s minutes apart: %s " % (m*min_factor,n))
				self.extrainfo += "Measurements done %s minutes apart: %s\n" % (m*min_factor,n)
				for p in popfiles:
					self.extrainfo += p+"\n"
		return self.__getsetpop(opts,tmpfiles)

	def Ttest(self, key):
		
		self.logger.info("Performing independent T-test on %s" % key)
		self.logger.info("Gathering mean %s-values" % key)
		
		dataset = ([],[],[],[],[],[],[])
		
		p_vals = []
		Yp_vals = []
		aA_vals = []
		aB_vals = []

		for x in self.SetA:
			# Welch's T-test
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
			ya,yb = np.array([]),np.array([])
			for y in self.SetA[x][key]['Y']:
				ya = np.append(ya,y)
			for y in self.SetB[x][key]['Y']:
				yb = np.append(yb,y)
			if not self.opts.nobs:
				Yab = [ya,yb]
			else:
				Yab = [[],[]]
				for n in np.array_split(ya,self.opts.nobs):
					#Yab[0].append(n.mean())
					Yab[0].append(gmean(abs(n)))
				for n in np.array_split(yb,self.opts.nobs):
					#Yab[1].append(n.mean())
					Yab[1].append(gmean(abs(n)))
			Yt,Yp = ttest_ind(np.asarray(Yab[0]),np.asarray(Yab[1]), equal_var=False)
			Yp_vals.append(Yp)
			if not Yp > 0:
				#print("Negaive P-value (%s) setting to 1.0." % Yp)
				Yp = 1.0
				#sys.exit()
			
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

			try:
				meanAlphaA = AlphaA[AlphaA > 0].mean()
			except FloatingPointError:
				meanAlphaA = np.inf
			try:
				meanAlphaB = AlphaB[AlphaB > 0].mean()
			except FloatingPointError:
				meanAlphaB = np.inf

			dataset[0].append(x)
			dataset[1].append(p)
			dataset[2].append(t)
			dataset[3].append(Yp)
			dataset[4].append(Yt)
			dataset[5].append(meanAlphaA)
			dataset[6].append(meanAlphaB)

		if self.opts.write:
			self.logger.info("Writing stats to %s" % self.opts.outfile+"_Ttest"+key+".txt")
			WriteStats(self.opts.out_dir, self.opts.outfile, dataset, \
				"Ttest"+key, ["Voltage", "P-value (Gmean)", "T-stat (Gmean)", "P-value (J)","T-stat (J)", "AlphaA", "AlphaB"])
		self.dataset[key] = dataset

		p_vals = np.array(p_vals)
		aA_vals = np.array(aA_vals)
		aB_vals = np.array(aB_vals)
		if self.opts.GUI:
			green,red,rs = '','',''
		else:
			green,red,rs = GREEN,RED,RS
		try:	
			if p_vals.mean() < 0.001: c = green
			else: c = red
			self.logger.info("p-value Mean: %s%s%s" % (green, p_vals.mean(), rs))
			if aA_vals[aA_vals > 0].mean() > 0.7: c = green
			else: c = red
			self.logger.info("α  SetA Mean: %s%s%s" % (c,aA_vals[aA_vals > 0].mean(),rs))
			if aB_vals[aB_vals > 0].mean() > 0.7: c = green
			else: c = red
			self.logger.info("α  SetB Mean: %s%s%s" % (c,aB_vals[aB_vals > 0].mean(),rs))
		except FloatingPointError:
			self.logger.error("Error outputting mean stats.")

	def TtestFN(self):
		x = list(self.SetA.keys())[-1]
		t_pos, p_pos = ttest_ind(self.SetA[x]['FN']['pos']['mean'], \
				self.SetB[x]['FN']['pos']['mean'], equal_var=False)	
		t_neg, p_neg = ttest_ind(self.SetA[x]['FN']['neg']['mean'], \
				self.SetB[x]['FN']['neg']['mean'], equal_var=False)
		#TODO Can't we figure this out at import?
		if self.opts.GUI:
			green,red,rs = '','',''
		else:
			green,red,rs = GREEN,RED,RS
		if p_pos < 0.001: c = green
		else: c = red
		self.logger.info("P-Value Vtrans(+): %s%s%s" % (c,str(p_pos),rs))
		
		if p_neg < 0.001: c = green
		else: c = red
		self.logger.info("P-Value Vtrans(–): %s%s%s" % (c,str(p_neg),rs))
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
			self.logger.debug("Division by zero error computing Alpha!")
			Alpha = np.inf
		except ValueError as msg:
			self.logger.debug("%s while computing Cronbach's Alpha." % str(msg))
			Alpha = np.inf
		except FloatingPointError as msg:
			self.logger.debug("%s while computing Cronbach's Alpha." % str(msg))
			Alpha = np.inf
		return Alpha
