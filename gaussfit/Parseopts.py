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
from getopt import gnu_getopt, GetoptError
from gaussfit.colors import *

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)

def ShowUsage():
	print('''
		%(path)s opts -o <outfile> <infiles>

		%(y)sThis program expects all X values in one column and all Y values in another.
		Ideally, feed it *_data.txt files and it will take V and J. It can extract X- and 
		Y-values from any two columns. Setting a compliance limit excludes Y > compliance
		from gaussian fits, but it rarely influences fits and can truncate data if set too low.
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
		-M	--min		Compute Vtrans from the min(Y) instead of the derivitve of the cubic spline
		-s	--skip		Skip plots with negative d2J/dV2 values at vcutoff for Vtrans calcuation
		-v	--vcutoff 	Voltage (absolute value) cut-off for dJ/dV skipping routine (default is =Vmin/Vmax)%(rs)s

	''' % {'path':os.path.basename(sys.argv[0]) ,'rs':RS,'y':YELLOW,'b':BLUE,'r':RED,'t':TEAL,'g':GREEN,'w':WHITE})
	sys.exit()
try:
	from scipy.optimize import curve_fit,OptimizeWarning
	from scipy.interpolate import UnivariateSpline
	import numpy as np
	# SciPy throws a useless warning for noisy J/V traces
	warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	ShowUsage()


class Opts:
	''' Parse commandline options and perform some basic checks.'''
	def __init__(self):
		try:
			opts, self.in_files = gnu_getopt(sys.argv[1:], 
					"hb:l:d:o:X:,Y:m:pc:nD:GMsv:", ["help" , "bins", 
					"loglevel=", "delimeter=","output=", "Xcol", 
					"Ycol", "maxr", "plot", "compliance", 
					"nowrite","dir:","GUI",
					"min","skip","vcutoff"])
		except GetoptError:
			print(RED+"Invalid option(s)"+RS)
			ShowUsage()
		
		# Set Defaults
		LOGLEVEL = logging.INFO
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
		self.skipohmic = False
		self.vcutoff = -1 # -1 will default to Vmin/Vmax
		# # #

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
				#print(RED+"\n\t\t> > > This feature is broken and -c is ignored! < < <"+RS)
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
                                self.smooth = False
			if opt in ('-s', '--skip'):
				self.skipohmic = True
			if opt in ('-v', '--vcutoff'):
				try:
					self.vcutoff = abs(float(arg))
				except ValueError:
					print(RED+"\n\t\t> > > vcutoff must be a number! < < <"+RS)
					ShowUsage()

		# Set up terminal logging. Set LOG to a file for debugging.
		LOG = False
		if LOG:
			logging.basicConfig(level=LOGLEVEL,format = '%(asctime)s %(process)d %(levelname)s %(message)s', filename=LOG,filemode='a+')
		else:
			logging.basicConfig(level=LOGLEVEL,format = GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE)

		if not len(self.in_files) and not self.GUI:
			print(RED+"\n\t\t> > > No input files! < < < "+RS)
			ShowUsage()
		# Setup CSV parser dialect
		csv.register_dialect('JV', delimiter=self.delim, 
				quoting=csv.QUOTE_MINIMAL)

