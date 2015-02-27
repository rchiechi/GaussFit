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

import os,csv
import numpy as np

class Writer():
	def __init__(self,parser):
		self.parser = parser
		self.opts = self.parser.opts

	def __getattr__(self, name):
		try:
    			return getattr(self.parser, name)
		except AttributeError as e:
    			raise AttributeError("Child' object has no attribute '%s'" % name)

	def WriteHistograms(self):
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Histograms.txt")
		with open(fn, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			headers = []
			for x in self.X: headers += ["Log |J| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
			writer.writerow(headers)
			for i in range(0, len( self.XY[list(self.XY.keys())[0]]['hist']['bin'] ) ):
				row = []
				for x in self.X: row += ["%0.4f"%self.XY[x]['hist']['bin'][i], 
						         "%0.2d"%self.XY[x]['hist']['freq'][i],
							 "%0.2d"%self.XY[x]['hist']['fit'][i]]
				writer.writerow(row)
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RHistograms.txt")
		with open(fn, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			headers = []
			for x in self.R['X']: headers += ["|R| (%0.4f)"%x, "Frequency (%0.4f)"%x, "Fit (%0.4f)"%x]
			writer.writerow(headers)
			for i in range(0, len( self.R[list(self.R.keys())[0]]['hist']['bin'] ) ):
				row = []
				for x in self.R['X']: row += ["%0.4f"%self.R[x]['hist']['bin'][i],
						         "%0.2d"%self.R[x]['hist']['freq'][i],
						         "%0.2d"%self.R[x]['hist']['fit'][i]]
				writer.writerow(row)


	def WriteVtrans(self):
		for key in ('pos', 'neg'):
			fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_Vtrans_"+key+".txt")
			with open(fn, 'w', newline='') as csvfile:
				writer = csv.writer(csvfile, dialect='JV')
				writer.writerow(["Vtrans (eV)","Frequency",
					"Gauss Fit (mean: %0.4f, variance: %f)"%(self.FN[key]['mean'], self.FN[key]['var'])])
				data = {}
				for i in range(0, len(self.FN[key]['bin'])):
					data[self.FN[key]['bin'][i]] = (self.FN[key]['freq'][i],self.FN[key]['fit'][i])
				for x in sorted(data.keys()):
					writer.writerow(['%0.4f'%x,'%d'%data[x][0],'%0.2f'%data[x][1]])
					
	def WriteFN(self):
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_FN.txt")
		with open(fn, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			writer.writerow(["1/V"] + ['Y_%d'%x for x in range(1,len( self.XY[list(self.XY.keys())[0]]['FN'] )+1)])
			for x in self.X:
				if x == 0.0:
					continue
				y = map(str, self.XY[x]['FN'])
				writer.writerow(["%0.4f\t" % (1/x)] + list(y))
	
	def WriteGauss(self):
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_Gauss.txt")
		with open(fn, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			writer.writerow(["Potential (V)","Log|J|","Variance"])
			Y = []
			Yerr = []
			for x in self.X:
				writer.writerow(['%f'%x,'%f'%self.XY[x]['hist']['mean'],'%f'%self.XY[x]['hist']['var']])
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_RGauss.txt")
		with open(fn, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			writer.writerow(["Potential (V)","|R|","Variance"])
			Y = []
			Yerr = []
			for x in self.R['X']:
				writer.writerow(['%f'%x,'%f'%self.R[x]['hist']['mean'],'%f'%self.R[x]['hist']['var']])


	def WriteData(self, log=False):
		if log:	key,label ='LogY','LogJ'
		else:	key, label ='Y','J'
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
		with open(fn,'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			writer.writerow(["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.XY[list(self.XY.keys())[0]][key] )+1)])
			for x in self.X:
				writer.writerow(["%0.4f"%x]+list(map(str,self.XY[x][key])))

	def WriteDJDV(self):
		label = 'DJDV'
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
		with open(fn,'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			writer.writerow(["Potential (V)"] + ['DJDV_%d'%x for x in range(1,len(self.DJDV[list(self.DJDV.keys())[0]])+1)])
			X = list(self.DJDV.keys())
			X.sort()
			for x in X:
				writer.writerow(["%0.4f"%x]+list(map(str,self.DJDV[x])))

	def WriteFiltered(self):
		label = 'filtered'
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
		with open(fn,'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			for l in self.filtered:
				writer.writerow(list(map(str,l)))

	def WriteRData(self):
		key,label = 'r', 'Rdata'
		fn = os.path.join(self.opts.out_dir,self.opts.outfile+"_"+label+".txt")
		with open(fn,'w', newline='') as csvfile:
			writer = csv.writer(csvfile, dialect='JV')
			writer.writerow(["Potential (V)"] + ['Y_%d'%x for x in range(1,len(self.R[list(self.R.keys())[0]][key] )+1)])
			for x in self.R['X']:
				writer.writerow(["%0.4f"%x]+list(map(str,self.R[x][key])))


class Plotter():
	def __init__(self,parser):
		self.parser = parser
		self.opts = self.parser.opts

	def __getattr__(self, name):
		try:
    			return getattr(self.parser, name)
		except AttributeError as e:
    			raise AttributeError("Child' object has no attribute '%s'" % name)

	def PlotData(self, key, ax, sym, **kw):
		xax = self.X
		if key == "FN":
			xax = 1/xax
			ax.set_title("Fowler Nordheim Plot of Initial Data")
			ax.set_xlabel(r'$V^{-1}$')
			ax.set_ylabel(r'$\mathregular{ln(\frac{J}{V^2})}$')
		if key == 'Y':
			if self.opts.compliance != np.inf: ax.set_ylim( (-1*self.opts.compliance, self.opts.compliance) )
			ax.set_title("Initial Data")
			ax.set_xlabel("Potenial (V)")
			ax.set_ylabel(r'Current Density ($A cm^{-2}$)')
		if key == 'LogY':
			ax.set_title("Semilog Plot of Initial Data")
			ax.set_xlabel("Potenial (V)")
			ax.set_ylabel(r'Current Density $\mathregular{log_{10}|J(\mathrm{A cm^{-2}})|}$')
		i = -1
		allY = np.array([])
		while True:
			i += 1
			try:
				allY = np.append(allY,[self.XY[x][key][i] for x in self.X])
				ax.plot(xax,[self.XY[x][key][i] for x in self.X], sym, **kw)
			except IndexError:
				break
		if key == 'LogY':
			ax.set_ylim(allY.min(),allY.max())
	def PlotDJDV(self,ax):
		xax = list(self.DJDV.keys())
		xax.sort()
		ax.set_title("Derivitive of Initial Data")
		ax.set_xlabel("Potential (V)")
		ax.set_ylabel(r'Normalized $\mathregular{\frac{dJ}{dV}}$')
		#ax.set_ylim(-0.05,0.2)
		i = -1
		while True:
			i += 1
			try:
				if i not in self.ohmic:
					ax.plot(xax, [self.DJDV[x][i] for x in xax], "-", lw=2)
				else:
					ax.plot(xax, [self.DJDV[x][i] for x in xax], "-", lw=0.5, color='grey')
			except IndexError:
				break

	def PlotHist(self,ax):
		ax.set_title("Gaussian Fit and Raw Data")
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
		fig.savefig(self.opts.outfile+"_fig.png", format="png")
