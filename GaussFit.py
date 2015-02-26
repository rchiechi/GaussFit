#!/usr/bin/python2.7
'''
Version: !!! DEPRECATED DO NOT USE !!!
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
import sys,os,logging,warnings
from getopt import gnu_getopt, GetoptError

print("DEPRECATED DO NOT USE!!! Use GaussFit3 instead!!!")
sys.exit()
try:
	from scipy.optimize import curve_fit
	import numpy as np
except ImportError as msg:
	print(msg)
	print("Error importing numpy/scipy, you probably need to install them.")

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)


YELLOW="\033[1;33m"
WHITE="\033[0m"
RED="\033[1;31m"
TEAL="\033[1;36m"
GREEN="\033[1;32m"
BLUE="\033[1;34m"
RS="\033[0m"
CL="\033[2K"


def ShowUsage():
	print '''
		%(path)s opts -o <outfile> <infiles>

		%(y)sThis program expects all X values in one column and all Y values in another.
		Ideally, feed it *_data.txt files and it will take V and J. It can extract X- and 
		Y-values from any two columns. Setting a compliance limit excludes Y > compliance
		from gaussian fits and rescales plots for display, but it rarely influences fits.%(rs)s

		%(b)scmd	command		help (default)%(rs)s
		%(g)s-h	--help		This help
		-b	--bins		Number of bins (default 50)
		-l	--loglevel 	Logging level (info)
		-d	--delimeter	Delimeter in input files (default: tab)
		-X	--Xcol		The column with X-values (default:1)
		-Y	--Ycol		The column wiht Y-Values (default:3)
		-o	--output	Outputfile (taken from first input)
		-p	--plot		Plot data save to a png file
		-n	--nowrite	Don't write output files (implies -p)
		-c	--compliance	set compliance limit for gaussian fits (default: inf)%(rs)s
	''' % {'path':os.path.basename(sys.argv[0]) ,'rs':RS,'y':YELLOW,'b':BLUE,'r':RED,'t':TEAL,'g':GREEN,'w':WHITE}
	sys.exit()

class Opts:

        def __init__(self):
                try:
			opts, self.in_files = gnu_getopt(sys.argv[1:], "hb:l:d:o:X:,Y:pc:n", ["help" , "bins", \
									"loglevel=",\
									"delimeter=","output=", "Xcol", "Ycol",
									"plot", "compliance", "nowrite"])
                except GetoptError:
                        error("Invalid option(s)")
                        ShowUsage()
		if not len(self.in_files):
			print RED+"\n\t\t> > > No input files! < < < "+RS
			ShowUsage()

		#Set Defaults
		
		LOGLEVEL = logging.INFO
		LOG = False
		self.Xcol = 0
		self.Ycol = 2
		self.bins = 50
		self.delim='\t'
		self.plot=False
		self.write=True
		self.compliance=np.inf

		if len(self.in_files):
			template = os.path.basename(self.in_files[0])
			if template[-4] == '.':
				template = template[:-4]
			self.outfile = template
		else:
			self.outfile = ""
		self.maxnumcol = 0

		for opt, arg in opts:
			if opt in ('-h', '--help'):
				ShowUsage()
				sys.exit()
			if opt in ('-b', '--bins'):
				self.bins =int(arg)
			if opt in ('-l', '--loglevel'):
				if arg.lower() == 'debug':
					LOGLEVEL = logging.DEBUG
				elif arg.lower() == 'warn':
					LOGLEVEL = logging.WARN
				elif arg.lower() == 'error':
					LOGLEVEL = logging.ERROR
				elif arg.logwer() == 'info':
					LOGLEVEL = logging.INFO
				else:
					print RED+"%s is not a logging level."+RS % arg
					ShowUsage()
					sys.exit()
			if opt in ('-d', '--delimeter'):
				self.delim = str(arg)
			if opt in ('-o', '--output'):
				self.outfile = arg
			if opt in ('-X', '--Xcol'):
				self.Xcol = int(arg)-1
			if opt in ('-Y', '--Ycol'):
				self.Ycol = int(arg)-1
			if opt in ('-p', '--plot'):
				self.plot=True
			if opt in ('-c', '--compliance'):
				self.compliance=float(arg)
			if opt in ('-n', '--nowrite'):
				self.write=False
				self.plot=True

		if LOG:
			logging.basicConfig(level=LOGLEVEL,format = '%(asctime)s %(process)d %(levelname)s %(message)s', filename=LOG,filemode='a+')
		else:
			logging.basicConfig(level=LOGLEVEL,format = GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE)

	def __numerror(self,s):
		print RED+s+" needs to be a number."+RS
		sys.exit()


class Parse():

	def __init__(self,opts):
		self.opts = opts
		self.parsed = {}
		self.XY ={}
		self.X = np.array([])
		self.FN = {}
		self.compliance_traces = []	
	def isfloat(self,f):
		try:
			float(f)
			return True
		except ValueError:
			return False
	def splitline(self,l):
		columns = []
		for n in l.strip('\n').strip('\r').split(self.opts.delim):
			if n:
				columns.append(n)
		return columns
	def tofloat(self, f):
		try:
			f = float(f)
		except ValueError:
			f = np.NAN
		return f
		
	def ReadFiles(self, fns):
		Ycol = self.opts.Ycol
		Xcol = self.opts.Xcol
		line_idx = 0
		uniqueX = {}
		for fn in fns:
			logging.info("Parsing %s%s%s", TEAL,fn,YELLOW)
			fh = open(fn, 'rt')
			firstline = self.splitline(fh.readline())
			numcol = len(firstline)
			if numcol < Ycol:
				logging.critical("Ycol > total number of columns! (%d > %d)", Ycol, numcol)
				logging.warn("%sNot importing:%s %s!", RED, YELLOW,fn)
				continue
			logging.debug("Found "+str(numcol)+" columns in "+fn)
			labels = []
			if self.isfloat(firstline[0]):
				logging.debug("No labels found in first row.")
				fh.seek(0)
			else:
				logging.debug("Assuming first row is column labels")
				for l in firstline:
					labels.append(l)
				logging.info("Y column is %s." , labels[Ycol])
			for l in fh.readlines():
			    row = self.splitline(l)
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

		for x in uniqueX:
			logging.debug('Pulling Y where X=%0.2f', x)
			y, fn = [], []
			idx = uniqueX[x].keys()
			idx.sort()
			for i in idx: y.append(uniqueX[x][i][0]), fn.append(uniqueX[x][i][1])
			y = np.array(y)
			logy = np.log10(abs(y))
			self.XY[x] = { "Y":y, \
				       "LogY":logy, \
				       "hist":self.dohistogram(logy,"J"), \
				       "FN":np.array(fn) }
		self.X = np.array(self.XY.keys())
		self.X.sort()
		self.FN["neg"], self.FN["pos"] = self.findmin()
		self.R = self.dorect()

	def dorect(self):
		R = {}
		for x in self.X:
			if x == 0:
				R[x]['r']=[1]
				R[x]['hist']={}
			if x < 0:
				continue
			if -1*x not in X:
				logging.warn("(Rectification) Didn't find a negative voltage for %d.", x)
				continue
			else:
				R[x]['r'] = []
				if len(self.XY[x]) != len(self.XY[-1*x]):
					logging.warn("(Rectification) Length of Y values differs for +/-%d.",x)
					continue
				for i in range(0, len(self.XY[x])):
					R[x]['r'].append( self.XY[x][i]/self.XY[-1*x][i] )
		for x in R:
			if x == 0:
				continue
			R[x]['hist'] = self.dohistogram(R[x]['r'],"R")
			print(R[x]['r'])
		return R

	def findmin(self):
		neg_min_x, pos_min_x = [],[]
		i = -1
		tossed = 0
		while True:
			i += 1
			try:
				pos,neg = {},{}
				for x in self.X:
					if abs(self.XY[x]['Y'][i]).max() >= self.opts.compliance:
						tossed += 1
						continue
					y = self.XY[x]['FN'][i]
					if x < 0: neg[y] = x
					if x > 0: pos[y] = x
				if not len(neg.keys()) or not len(pos.keys()):
					logging.warn("Skipping empty column in FN calculation.")
					continue
				neg_min_x.append(neg[np.nanmin(neg.keys())])
				pos_min_x.append(pos[np.nanmin(pos.keys())])
			except IndexError:
				break
		
		if tossed: logging.warn("Tossed %d compliance traces during FN calculation.", tossed)
		neg_min_x = np.array(neg_min_x)
		pos_min_x = np.array(pos_min_x)
		return self.dohistogram(neg_min_x, "Vtrans(-)"), self.dohistogram(pos_min_x, "Vtrans(+)")

	def gauss(self, x, *p):
                A, mu, sigma = p
                return A*np.exp(-(x-mu)**2/(2.*sigma**2))
        
        def dohistogram(self, Y, label=""):
		y_compliance = np.nonzero(abs(Y) > self.opts.compliance)
		if len(y_compliance[0]) and label == "J":
			logging.warn("Tossing %d data points for %s histogram!", len(y_compliance[0]), label)
			Y = Y[y_compliance]
		logging.debug("%d points to consider.", len(Y))
		freq, bins = np.histogram(Y, bins=self.opts.bins, density=False)      
                p0 = [1., Y.mean(), Y.std()]
                bin_centers = (bins[:-1] + bins[1:])/2
		try:
                    coeff, var_matrix = curve_fit(self.gauss, bin_centers, freq, p0=p0, maxfev=10000)
                    hist_fit = self.gauss(bin_centers, *coeff)
                except RuntimeError as msg:
                    logging.warning("|%s| Fit did not converge (%s)", label, str(msg), exc_info=False)
                    coeff = p0
                    hist_fit = np.array([x for x in range(0, len(bin_centers))])
		return {"bin":bin_centers, "freq":freq, "mean":coeff[1], "std":coeff[2], "var":coeff[2]**2, "bins":bins, "fit":hist_fit}

	def PrintFN(self):
		for key in ('pos', 'neg'):
			print("|Vtrans %s| mean: %0.4f variance: %f" % (key, self.FN[key]['mean'], self.FN[key]['var']) )

	def WriteHistograms(self):
		ofh = file(self.opts.outfile+"_Histograms.txt", 'wb')
		for x in self.X: ofh.write("\t".join( ("Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x) )+"\t" )
		ofh.write('\n')
		for i in range(0, len( self.XY[self.XY.keys().index(0)]['hist']['bin'] ) ):
			for x in self.X: ofh.write("\t".join( ("%0.4f"%self.XY[x]['hist']['bin'][i], \
						"%0.2d"%self.XY[x]['hist']['freq'][i], \
						"%0.2d"%self.XY[x]['hist']['fit'][i]) )+"\t" )
			ofh.write('\n')
		ofh.close()

	def WriteVtrans(self):
		for key in ('pos', 'neg'):
			ofh = file(self.opts.outfile+"_Vtrans_"+key+".txt", 'wb')
			ofh.write("\t".join( ("Vtrans (eV)", \
					"Frequency", \
					"Gauss Fit (mean: %0.4f, variance: %f)"%(self.FN[key]['mean'], self.FN[key]['var'])) )+"\n")
			data = {}
			for i in range(0, len(self.FN[key]['bin'])):
				data[self.FN[key]['bin'][i]] = (self.FN[key]['freq'][i],self.FN[key]['fit'][i])
			v = data.keys()
			v.sort()
			for x in v:
				ofh.write("\t".join(('%0.4f'%x,'%d'%data[x][0],'%0.2f'%data[x][1]))+"\n")
	def WriteFN(self):
		ofh = file(self.opts.outfile+"_FN.txt", 'wb')
		cols = ["1/V"] + ['Y_%d'%x for x in range(1,len(self.XY[self.XY.keys().index(0)]['FN'] )+1)]
		ofh.write("\t".join(cols)+"\n")
		for x in self.X:
			if x == 0.0:
				continue
			y = map(str, self.XY[x]['FN'])
			ofh.write("%0.4f\t" % (1/x))
			ofh.write("\t".join(y)+"\n")
		ofh.close()

	def WriteGauss(self):
		ofh = file(self.opts.outfile+"_Gauss.txt", 'wb')
		ofh.write("\t".join(("Potential (V)","Log|J|","Variance"))+"\n")
		Y = []
		Yerr = []
		for x in self.X:
			ofh.write("\t".join(('%f'%x,'%f'%self.XY[x]['hist']['mean'],'%f'%self.XY[x]['hist']['var']))+"\n")
		ofh.close()

	def WriteData(self, log=False):
		if log:	key,label ='LogY','LogJ'
		else:	key, label ='Y','data'
		ofh = file(self.opts.outfile+"_"+label+".txt", 'wb')
		cols = ["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.XY[self.XY.keys().index(0)][key] )+1)]
		ofh.write("\t".join(cols)+"\n")
		for x in self.X:
			ofh.write("\t".join(["%0.4f"%x]+map(str,self.XY[x][key]))+"\n")
		ofh.close()

	def PlotData(self, key, ax, sym, **kw):
		xax = self.X
		if key == "FN":
			xax = 1/xax
			ax.set_title("Fowler Nordheim Plot of Initial Data")
			ax.set_xlabel("1/V")
			ax.set_ylabel("ln(I/V^2)")
		if key == 'Y':
			if self.opts.compliance != np.inf: ax.set_ylim( (-1*self.opts.compliance, self.opts.compliance) )
			ax.set_title("Initial Data")
			ax.set_xlabel("Potenial (V)")
			ax.set_ylabel("Current Density (A/cm^2)")
		if key == 'LogY':
			#if self.opts.compliance != np.inf:ax.set_ylim( ymax=np.log10(self.opts.compliance) )
			ax.set_title("Semilog Plot of Initial Data")
			ax.set_xlabel("Potenial (V)")
			ax.set_ylabel("Current Density log10(A/cm^2)")
		i = -1
		while True:
			i += 1
			try:
				ax.plot(xax,[self.XY[x][key][i] for x in self.X], sym, **kw)
			except IndexError:
				break

	def PlotHist(self,ax):
		ax.set_title("Gaussian Fit and Intial Data")
		ax.set_xlabel('Potential (V)')
		ax.set_ylabel('Current Density (Log10 |A/cm^2|)')
		Y, Yerr = [],[]
		for x in self.X:
			Y.append(self.XY[x]["hist"]["mean"])
			Yerr.append(self.XY[x]["hist"]["var"])
		ax.errorbar(self.X, Y, yerr=Yerr, lw=3.0, color='k')

	def PlotVtrans(self,ax):
		ax.set_title("Histogram and fit of Vtrans")
		ax.set_xlabel('Vtrans')
		ax.set_ylabel('Frequency')
		for key in ('pos','neg'):
			ax.bar(self.FN[key]['bin'], self.FN[key]['freq'], width=0.01, color='g')
			ax.plot(self.FN[key]['bin'], self.FN[key]['fit'], lw=2.0, color='b', label='Fit')

if __name__ == "__main__":
	opts = Opts()
	parser = Parse(opts)
	parser.ReadFiles(opts.in_files)
	if opts.write:	
		logging.info("Writing files...")
		parser.WriteVtrans()
		parser.WriteFN()
		parser.WriteGauss()
		parser.WriteData()
		parser.WriteData(True)
		parser.WriteHistograms()
	if opts.plot:
		logging.info("Generating plots...")
		try:
			import matplotlib.pyplot as plt
		except:
			logging.error("Cannot import matplotlib!", exc_info=True)
		fig = plt.figure(figsize=(15,10))
		ax1 = fig.add_axes([0.05, 0.55, 0.4, 0.4])
		ax2 = fig.add_axes([0.55, 0.55, 0.4, 0.4])
		ax3 = fig.add_axes([0.05, 0.05, 0.4, 0.4])
		ax4 = fig.add_axes([0.55, 0.05, 0.4, 0.4])
		parser.PlotData('Y', ax1, '-')
		parser.PlotData('LogY',ax2,':',lw=0.25, color='c')
		parser.PlotData('FN', ax3, 'x', ms=2)
		parser.PlotHist(ax2)
		parser.PlotVtrans(ax4)
		fig.savefig(opts.outfile+"_fig.png", format="png")
		plt.show()
	parser.PrintFN()
	sys.exit()
