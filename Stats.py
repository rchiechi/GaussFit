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
#from gaussfit.Parser import Parse
#from gaussfit.Output import Writer

import sys,os,logging,warnings,csv
from gaussfit.colors import *

#try:
#	from scipy.optimize import curve_fit,OptimizeWarning
#	import scipy.interpolate 
#	from scipy.stats import gmean,kstest,ttest_ind,ttest_rel,ttest_1samp
#	import scipy.misc 
#	import numpy as np
#	# SciPy throws a useless warning for noisy J/V traces
#	warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)
#
#except ImportError as msg:
#	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
#	sys.exit()

#warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*',UserWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in log.*',RuntimeWarning)
#warnings.filterwarnings('ignore','.*invalid value encountered in true_divide.*',RuntimeWarning)
##warnings.filterwarnings('ignore','.*impossible result.*',UserWarning)

from gaussfit.stats import Stats

def Go(opts):
	'''
	Call this function to execute the parsing engine
	i.e., as main()
	'''

	if len(opts.setA) and len(opts.setB):
		s = Stats(opts)

	elif len(opts.setA):
		SigTest(opts)
	else:
		logging.error("No input files to parse!")


opts = Opts
if opts.GUI:
		from gui import statsfilebrowser
		gui = statsfilebrowser.ChooseFiles(opts,Go)
else:
		Go(opts)
