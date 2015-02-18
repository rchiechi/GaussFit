#!/usr/bin/env python3

YELLOW="\033[1;33m"
WHITE="\033[0m"
RED="\033[1;31m"
TEAL="\033[1;36m"
GREEN="\033[1;32m"
BLUE="\033[1;34m"
RS="\033[0m"
CL="\033[2K"

import sys,os,logging,warnings
from getopt import gnu_getopt, GetoptError

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)

def ShowUsage():
	print('''
		%(path)s opts -o <outfile> <infiles>

		%(y)sThis program expects all X values in one column and all Y values in another.
		Ideally, feed it *_data.txt files and it will take V and J. It can extract X- and 
		Y-values from any two columns. Setting a compliance limit excludes Y > compliance
		from gaussian fits and rescales plots for display, but it rarely influences fits.
		The maxr parameter is needed to filter out huge values that occur when a junction shorts.%(rs)s

		%(b)scmd	command		help (default)%(rs)s
		%(g)s-h	--help		This help
		-b	--bins		Number of bins (default 50)
		-l	--loglevel 	Logging level (info)
		-d	--delimeter	Delimeter in input files (default: tab)
		-X	--Xcol		The column with X-values (default:1)
		-Y	--Ycol		The column wiht Y-Values (default:3)
		-m	--maxr		Maximum allowable value of R (default:10)
		-o	--output	Outputfile (taken from first input)
		-p	--plot		Plot data save to a png file
		-n	--nowrite	Don't write output files (implies -p)
		-c	--compliance	set compliance limit for gaussian fits (default: inf)
		-D      --dir           Output directory (merged with --output)
		-G      --GUI           Launch the GUI
		-M	--min		Compute Vtrans from the min(Y) instead of the derivitve of the cubic spline%(rs)s

	''' % {'path':os.path.basename(sys.argv[0]) ,'rs':RS,'y':YELLOW,'b':BLUE,'r':RED,'t':TEAL,'g':GREEN,'w':WHITE})
	sys.exit()
try:
	from scipy.optimize import curve_fit
	from scipy.interpolate import UnivariateSpline
	import numpy as np

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	ShowUsage()


class Opts:

	def __init__(self):
		try:
			opts, self.in_files = gnu_getopt(sys.argv[1:], "hb:l:d:o:X:,Y:m:pc:nD:GM", ["help" , "bins", \
									"loglevel=",\
									"delimeter=","output=", "Xcol", "Ycol", "maxr"
									"plot", "compliance", "nowrite","dir:","GUI","min"])
		except GetoptError:
			error("Invalid option(s)")
			ShowUsage()
		

		#Set Defaults
		
		LOGLEVEL = logging.INFO
		LOG = False
		self.Xcol = 0
		self.Ycol = 2
		self.maxr = 10.0
		self.bins = 50
		self.delim='\t'
		self.plot=False
		self.write=True
		self.compliance=np.inf
		self.GUI=False
		self.out_dir = os.environ['PWD']
		self.smooth = True

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
					print(RED+"%s is not a logging level."+RS % arg)
					ShowUsage()
			if opt in ('-d', '--delimeter'):
				self.delim = str(arg)
			if opt in ('-o', '--output'):
				self.outfile = arg
			if opt in ('-X', '--Xcol'):
				self.Xcol = int(arg)-1
			if opt in ('-Y', '--Ycol'):
				self.Ycol = int(arg)-1
			if opt in ('-m', '--maxr'):
				self.maxr = float(arg)
			if opt in ('-p', '--plot'):
				self.plot=True
			if opt in ('-c', '--compliance'):
				self.compliance=float(arg)
			if opt in ('-n', '--nowrite'):
				self.write=False
				self.plot=True
			if opt in ('-d', '--dir'):
                                if os.path.exists(arg):
                                    self.out_dir = arg
                                else:
                                    print(RED+"\n\t\t> > > Output directory "+arg+" does not exist! < < <"+RS)
                                    ShowUsage()
			if opt in ('-G', '--GUI'):
                                self.GUI = True
			if opt in ('-M', '--min'):
                                self.smooth = True

##                self.logger = logging.getLogger('CLI')
		if LOG:
			logging.basicConfig(level=LOGLEVEL,format = '%(asctime)s %(process)d %(levelname)s %(message)s', filename=LOG,filemode='a+')
##		elif self.GUI:
##			print("Loading GUI settings")
##			logging.basicConfig(level=LOGLEVEL,format = os.path.basename(sys.argv[0])+' %(levelname)s %(message)s')
		else:
			logging.basicConfig(level=LOGLEVEL,format = GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE)

		if not len(self.in_files) and not self.GUI:
			print(RED+"\n\t\t> > > No input files! < < < "+RS)
			ShowUsage()

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
				if x in uniqueX.keys(): uniqueX[x][line_idx]=(y,z)
				else: uniqueX[x]={line_idx:(y,z)}
				line_idx += 1

		for x in uniqueX:
			logging.debug('Pulling Y where X=%0.2f', x)
			y, fn = [], []
			for i in sorted(uniqueX[x].keys()): y.append(uniqueX[x][i][0]), fn.append(uniqueX[x][i][1])
			y = np.array(y)
			logy = np.log10(abs(y))
			self.XY[x] = { "Y":y, \
				   "LogY":logy, \
				   "hist":self.dohistogram(logy,"J"), \
				   "FN":np.array(fn) }
		self.X = np.array(sorted(self.XY.keys()))
		self.X.sort()
		self.FN["neg"], self.FN["pos"] = self.findmin()
		self.R = self.dorect()
		self.DJDV = self.dodjdv()

	def dorect(self):
		R = {}
		for x in self.X:
			if x <= 0: continue
			elif -1*x not in self.X:
				logging.warn("(Rectification) Didn't find a negative voltage for %d.", x)
				continue
			ypos, yneg = self.XY[x]["Y"], self.XY[-1*x]["Y"]
			if len(ypos) != len(yneg):
				logging.warn("(Rectification) Length of Y values differs for +/-%d.",x)
				continue
			r = []
			for i in range(0, len(ypos)):
				_r = abs(ypos[i]/yneg[i])
				#if _r > self.opts.maxr:
				#	logging.warn("Rejecting R=%0.4f at V=%0.2f because it exceeds maxR (%0.1f)", _r, x, self.opts.maxr)
				#	continue
				r.append( abs(ypos[i]/yneg[i]) )	
			R[x]={'r': np.array(r)}
		for x in R:
			R[x]['hist'] = self.dohistogram(R[x]['r'],"R")
			#print(R[x]['hist'])
		R['X'] = np.array(sorted(R.keys()))
		return R
	
	def dodjdv(self):
		spls = {}
		xs = np.linspace(self.X.min(), self.X.max(), 100)
		for x in xs:
			spls[x] = []
		i = -1
		while True:
			try:
				i += 1
				y = []
				for x in self.X:
					y.append(self.XY[x]['Y'][i])
				spl = UnivariateSpline( self.X, y, k=4).derivative()
				for x in spls:
					spls[x].append(spl(x))
			except IndexError:
				break
		return spls

	def findmin(self):
		neg_min_x, pos_min_x = [],[]
		i = -1
		tossed = 0
		while True:
			i += 1
			try:
				pos,neg = {},{}
				x_neg, y_neg, x_pos, y_pos = [],[],[],[]
				for x in self.X:
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

				splneg = UnivariateSpline( x_neg, y_neg, k=4  )
				splpos = UnivariateSpline( x_pos, y_pos, k=4 )
				#print("Spline:",splneg(x_neg))
				#print("Data:",y_neg)
				#print("DI/DV roots (-):",splneg.derivative().roots())
				#print("DI/DV roots (+):",splpos.derivative().roots())
				if self.opts.smooth:
					rootneg = np.nanmin(splneg.derivative().roots())
					neg_min_x.append(rootneg[0])
					rootpos = np.nanmin(splpos.derivative().roots())
					pos_min_x.append(rootpos[0])
					if not rootneg[0] or rootpos[0]:
						logging.warn("No minimum found in FN derivative.")
					logging.debug("FN neg roots: "+str(len(rootneg)))
					logging.debug("FN pos roots: "+str(len(rootpos)))
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
		A, mu, sigma = p
		return A*np.exp(-(x-mu)**2/(2.*sigma**2))
	
	def dohistogram(self, Y, label=""):
		j_compliance = np.nonzero(abs(Y) > self.opts.compliance)
		r_compliance = np.nonzero(abs(Y) > self.opts.maxr)
		if len(j_compliance[0]) and label == "J":
			logging.warn("Tossing %d data points > %0.1f for %s histogram!", len(j_compliance[0]), self.opts.compliance, label)
			Y = Y[np.nonzero(abs(Y) <= self.opts.compliance)]
		if len(r_compliance[0]) and label == "R":
			logging.warn("Tossing %d data points > %0.1f for %s histogram!", len(r_compliance[0]), self.opts.maxr, label)
			Y = Y[np.nonzero(abs(Y) <= self.opts.maxr)]
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
			hist_fit = np.array([x*0 for x in range(0, len(bin_centers))])
		return {"bin":bin_centers, "freq":freq, "mean":coeff[1], "std":coeff[2], "var":coeff[2]**2, "bins":bins, "fit":hist_fit}

	def PrintFN(self):
		for key in ('pos', 'neg'):
			print("|Vtrans %s| mean: %0.4f variance: %f" % (key, self.FN[key]['mean'], self.FN[key]['var']) )

	def WriteHistograms(self):
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_Histograms.txt"), 'wt')
		for x in self.X: ofh.write("\t".join( ("Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x) )+"\t" )
		ofh.write('\n')
		for i in range(0, len( self.XY[list(self.XY.keys())[0]]['hist']['bin'] ) ):
			for x in self.X: ofh.write("\t".join( ("%0.4f"%self.XY[x]['hist']['bin'][i], \
						"%0.2d"%self.XY[x]['hist']['freq'][i], \
						"%0.2d"%self.XY[x]['hist']['fit'][i]) )+"\t" )
			ofh.write('\n')
		ofh.close()
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_R_Histograms.txt"), 'wt')
		for x in self.R['X']: ofh.write("\t".join( ("|R| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x) )+"\t" )
		ofh.write('\n')
		for i in range(0, len( self.R[list(self.R.keys())[0]]['hist']['bin'] ) ):
			for x in self.R['X']: ofh.write("\t".join( ("%0.4f"%self.R[x]['hist']['bin'][i], \
						"%0.2d"%self.R[x]['hist']['freq'][i], \
						"%0.2d"%self.R[x]['hist']['fit'][i]) )+"\t" )
			ofh.write('\n')
		ofh.close()

	def WriteVtrans(self):
		for key in ('pos', 'neg'):
			ofh = open(os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_"+key+".txt"), 'wt')
			ofh.write("\t".join( ("Vtrans (eV)", \
					"Frequency", \
					"Gauss Fit (mean: %0.4f, variance: %f)"%(self.FN[key]['mean'], self.FN[key]['var'])) )+"\n")
			data = {}
			for i in range(0, len(self.FN[key]['bin'])):
				data[self.FN[key]['bin'][i]] = (self.FN[key]['freq'][i],self.FN[key]['fit'][i])
			for x in sorted(data.keys()):
				ofh.write("\t".join(('%0.4f'%x,'%d'%data[x][0],'%0.2f'%data[x][1]))+"\n")
	def WriteFN(self):
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_FN.txt"), 'wt')
		cols = ["1/V"] + ['Y_%d'%x for x in range(1,len( self.XY[list(self.XY.keys())[0]]['FN'] )+1)]
		ofh.write("\t".join(cols)+"\n")
		for x in self.X:
			if x == 0.0:
				continue
			y = map(str, self.XY[x]['FN'])
			ofh.write("%0.4f\t" % (1/x))
			ofh.write("\t".join(y)+"\n")
		ofh.close()

	def WriteGauss(self):
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_Gauss.txt"), 'wt')
		ofh.write("\t".join(("Potential (V)","Log|J|","Variance"))+"\n")
		Y = []
		Yerr = []
		for x in self.X:
			ofh.write("\t".join(('%f'%x,'%f'%self.XY[x]['hist']['mean'],'%f'%self.XY[x]['hist']['var']))+"\n")
		ofh.close()
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_R_Gauss.txt"), 'wt')
		ofh.write("\t".join(("Potential (V)","|R|","Variance"))+"\n")
		Y = []
		Yerr = []
		for x in self.R['X']:
			ofh.write("\t".join(('%f'%x,'%f'%self.R[x]['hist']['mean'],'%f'%self.R[x]['hist']['var']))+"\n")
		ofh.close()

	def WriteData(self, log=False):
		if log:	key,label ='LogY','LogJ'
		else:	key, label ='Y','data'
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt"), 'wt')
		cols = ["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.XY[list(self.XY.keys())[0]][key] )+1)]
		ofh.write("\t".join(cols)+"\n")
		for x in self.X:
			ofh.write("\t".join(["%0.4f"%x]+list(map(str,self.XY[x][key])))+"\n")
		ofh.close()

	def WriteDJDV(self):
		#if log:	key,label ='LogY','LogJ'
		#else:	key, label ='Y','data'
		label = 'DJDV'
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt"), 'wt')
		cols = ["Potential (V)"] + ['DJDV_%d'%x for x in range(1,len(self.DJDV[list(self.DJDV.keys())[0]])+1)]
		ofh.write("\t".join(cols)+"\n")
		X = list(self.DJDV.keys())
		X.sort()
		for x in X:
			ofh.write("\t".join(["%0.4f"%x]+list(map(str,self.DJDV[x])))+"\n")
		ofh.close()

	def WriteRData(self):
		key,label = 'r', 'R_data'
		ofh = open(os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt"), 'wt')
		cols = ["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.R[list(self.R.keys())[0]][key] )+1)]
		ofh.write("\t".join(cols)+"\n")
		for x in self.R['X']:
			ofh.write("\t".join(["%0.4f"%x]+list(map(str,self.R[x][key])))+"\n")
		ofh.close()

	def PlotData(self, key, ax, sym, **kw):
		xax = self.X
		if key == "FN":
			xax = 1/xax
			ax.set_title("Fowler Nordheim Plot of Initial Data")
			ax.set_xlabel(r'$V^{-1}$')
			ax.set_ylabel(r'$ln(\frac{I}{V^2})$')
		if key == 'Y':
			if self.opts.compliance != np.inf: ax.set_ylim( (-1*self.opts.compliance, self.opts.compliance) )
			ax.set_title("Initial Data")
			ax.set_xlabel("Potenial (V)")
			ax.set_ylabel(r'Current Density ($A cm^{-2}$)')
		if key == 'LogY':
			ax.set_title("Semilog Plot of Initial Data")
			ax.set_xlabel("Potenial (V)")
			ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm^{-2}})|$')
		i = -1
		while True:
			i += 1
			try:
				ax.plot(xax,[self.XY[x][key][i] for x in self.X], sym, **kw)
			except IndexError:
				break
	def PlotDJDV(self,ax):
		xax = list(self.DJDV.keys())
		xax.sort()
		ax.set_title("Derivitive of Initial Data")
		ax.set_xlabel("Potential (V)")
		ax.set_ylabel(r'$\frac{dJ}{dV}$')
		ax.set_ylim(-0.05,0.25)
		i = -1
		while True:
			i += 1
			try:
				ax.plot(xax, [self.DJDV[x][i] for x in xax], "-")
			except IndexError:
				break

	def PlotHist(self,ax):
		ax.set_title("Gaussian Fit and Intial Data")
		ax.set_xlabel('Potential (V)')
		ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm}^{-2})|$')
		Y, Yerr = [],[]
		for x in self.X:
			Y.append(self.XY[x]["hist"]["mean"])
			Yerr.append(self.XY[x]["hist"]["var"])
		ax.errorbar(self.X, Y, yerr=Yerr, lw=3.0, color='k')

	def PlotVtrans(self,ax):
		ax.set_title(r'Histogram and fit of $V_{trans}$')
		ax.set_xlabel(r'$V_{trans}$')
		ax.set_ylabel('Counts')
		for key in ('pos','neg'):
			ax.bar(self.FN[key]['bin'], self.FN[key]['freq'], width=0.01, color='g')
			ax.plot(self.FN[key]['bin'], self.FN[key]['fit'], lw=2.0, color='b', label='Fit')

	def DoPlots(self, plt):
		fig = plt.figure(figsize=(16,10))
		ax1 = fig.add_axes([0.06, 0.55, 0.4, 0.4])
		ax2 = fig.add_axes([0.56, 0.55, 0.4, 0.4])
		ax3 = fig.add_axes([0.06, 0.05, 0.4, 0.4])
		ax4 = fig.add_axes([0.56, 0.05, 0.4, 0.4])
		#self.PlotData('Y', ax1, '-')
		self.PlotDJDV(ax1)
		self.PlotData('LogY',ax2,':',lw=0.25, color='c')
		self.PlotData('FN', ax3, 'x', ms=2)
		self.PlotHist(ax2)
		self.PlotVtrans(ax4)
		fig.savefig(opts.outfile+"_fig.png", format="png")


def Go(opts):
		parser = Parse(opts)
		parser.ReadFiles(opts.in_files)
		if opts.write:	
				logging.info("Writing files...")
				parser.WriteVtrans()
				parser.WriteFN()
				parser.WriteGauss()
				parser.WriteData()
				parser.WriteDJDV()
				parser.WriteData(True)
				parser.WriteRData()
				parser.WriteHistograms()
		if opts.plot:
				logging.info("Generating plots...")
				try:
						import matplotlib.pyplot as plt
						parser.DoPlots(plt)
						plt.show()
				except ImportError as msg:
						logging.error("Cannot import matplotlib! %s", str(msg), exc_info=False)
		parser.PrintFN()
	

if __name__ == "__main__":
		opts = Opts()
		if opts.GUI:
				from gui import filebrowser
				gui = filebrowser.ChooseFiles(opts,Go)

		else:
				Go(opts)
