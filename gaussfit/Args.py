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
import argparse
from gaussfit.colors import *

warnings.filterwarnings('ignore','.*divide by zero.*',RuntimeWarning)

try:
	from scipy.optimize import curve_fit,OptimizeWarning
	from scipy.interpolate import UnivariateSpline
	import numpy as np
	# SciPy throws a useless warning for noisy J/V traces
	warnings.filterwarnings('ignore','.*Covariance of the parameters.*',OptimizeWarning)

except ImportError as msg:
	print("\n\t\t%s> > > Error importing numpy/scipy! %s%s%s < < <%s" % (RED,RS,str(msg),RED,RS))
	sys.exit()

desc='''
	This program expects all X values in one column and all Y values in another.
	Ideally, feed it *_data.txt files and it will take V and J. It can extract X- and 
	Y-values from any two columns. Setting a compliance limit excludes Y > compliance
	from gaussian fits, but it rarely influences fits and can truncate data if set too low.
	The maxr parameter is needed to filter out huge values that occur when a junction shorts.
     '''

parser = argparse.ArgumentParser(description=desc,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('in_files', metavar='Files-to-parse', type=str, nargs='*', default=[], 
		help='Datafiles to parse.')
parser.add_argument('-G','--GUI', action='store_true', default=False,
                help="Launch the GUI.")
parser.add_argument('-o','--outfile', metavar="OUTPUT FILE", default="",
                help="Outputfile (taken from first input)")
parser.add_argument('-D','--dir', metavar="OUTPUT DIR", dest='out_dir', default=os.path.join(os.getcwd(),'parsed'),
                help="Output directory (combined with -o).")
parser.add_argument('-l','--loglevel', default='info', choices=('info','warn','error','debug'),
                help="Set the logging level.")
parser.add_argument('-p','--plot', action='store_true', default=False,
                help="Plot data and save to png file.")
parser.add_argument('-n','--nowrite', dest='write', action='store_false', default=True,
                help="Do not write output files (implies -p).")
parser.add_argument('-d','--delim', default='tab', choices=('tab', 'comma', 'space'), 
		help="Delimeter in inputfiles.")
parser.add_argument('-X','--Xcol', type=int, default=1, 
		help="Column to treat as X.")
parser.add_argument('-Y','--Ycol', type=int, default=3, 
		help="Column to treat as Y. Set to 0 to grab all columns except Xcol.")
parser.add_argument('-b','--bins', default=50, type=int,
                help='Number of bins for histograms (except heatmap).')
parser.add_argument('-m','--maxr', type=float, default=2, 
		help="Maximum allowable value of log|R| or (R if -L).")
parser.add_argument('-c','--compliance', default=np.inf, type=float, 
		help="Set compliance limit for gaussian fits.")
parser.add_argument('-M','--minfn', action='store_false', dest='nomin', default=True, 
		help="Compute Vtrans from the min(Y) instead of the derivitve of the cubic spline.")
parser.add_argument('-s','--skip',  action='store_true', dest='skipohmic', default=False, 
		help="Skip plots with negative d2J/dV2 values at vcutoff for Vtrans calcuation and conductance plot.")
parser.add_argument('-S','--smooth', type=float, default=0, 
		help="Cutoff for residuals when smoothing splines for dJ/dV 0 to disable. 1e-4 for artifacts, 1e-12 for smooth plots.")
parser.add_argument('-v','--vcutoff', type=float, default=-1, 
		help="Voltage (absolute value) cut-off for dJ/dV skipping routine (-1 for Vmin/Vmax)")
parser.add_argument('-a','--lower', metavar='LOWER', dest='mlow', type=float, default=-6, 
		help="Lower cutoff value for conductance heat map plot.")
parser.add_argument('-z','--upper', metavar='UPPER', dest='mhi', type=float, default=0, 
		help="Upper cutoff value for conductance heat map plot.")
parser.add_argument('-B','--heatmapbins', default=25, type=int, 
		help="Number of bins for the conductance heatmap plot.")
parser.add_argument('-R','--logr', default=True, action='store_false', 
		help="Compute |R| instead of log|R| for histograms of rectification.")
parser.add_argument('-L','--lorenzian', default=False, action='store_true', 
		help="Fit data to a Lorenzian instead of a Gaussian.")
parser.add_argument('--maxfev', type=int, default=10000, 
		help="Maximum interations for fitting histograms.")
parser.add_argument('--plots', metavar='setA', type=str, nargs='*', default=['J'], 
		help="What to plot: J (default),R.")


Opts=parser.parse_args()

if not Opts.write:
	Opts.plot = True

if len(Opts.in_files) and not Opts.outfile:
	if Opts.in_files[0][-4] == '.':
		Opts.outfile = os.path.basename(Opts.in_files[0])[:-4]
	else:
		Opts.outfile = Opts.in_files[0]
Opts.Xcol -= 1
Opts.Ycol -= 1
			
if not os.path.exists(Opts.out_dir):
	os.mkdir(Opts.out_dir)
	print("Creating %s" % Opts.out_dir)
	#parser.print_help()
	#print(RED+"\n\t\t> > > Output directory "+arg+" does not exist! < < <"+RS)
	#sys.exit()

# Set up terminal logging. Set LOG to a file for debugging.
LOG = False
if LOG:
	logging.basicConfig(level= getattr(logging,Opts.loglevel.upper()),format = '%(asctime)s %(process)d %(levelname)s %(message)s', filename=LOG,filemode='a+')
else:
	logging.basicConfig(level= getattr(logging,Opts.loglevel.upper()),format = GREEN+os.path.basename(sys.argv[0]+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE)

if not len(Opts.in_files) and not Opts.GUI:
	parser.print_help()
	print(RED+"\n\t\t> > > No input files! < < < "+RS)
	sys.exit()


if Opts.delim == 'tab':
	delim='\t'
if Opts.delim == 'comma':
	delim=','
if Opts.delim == 'space':
	delim=' '

# Setup CSV parser dialect
csv.register_dialect('JV', delimiter=delim, quoting=csv.QUOTE_MINIMAL)
