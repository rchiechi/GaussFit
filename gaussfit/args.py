#!/usr/bin/env python3
'''
Copyright (C) 2018 Ryan Chiechi <r.c.chiechi@rug.nl>
Description:
        Parse all of the command line arguments for GaussFit.py
        and Stats.py

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
		help='Datafiles to parse (GaussFit mode).')
parser.add_argument('-A','--setA', metavar='setA', type=str, nargs='*', default=[], 
		help='Datafiles to parse for set A (Stats mode).')
parser.add_argument('-B','--setB',  metavar='setB', type=str, nargs='*', default=[], 
		help='Datafiles to parse for set B (Stats mode).')

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
		help="Maximum allowable value of log|R| or (R if -R).")
parser.add_argument('-c','--compliance', default=np.inf, type=float, 
		help="Set compliance limit for gaussian fits.")
parser.add_argument('--interpolateminfn', action='store_true', default=False, 
		help="Compute Vtrans from the minimum of the derivative of the FN plot.")
parser.add_argument('-s','--skip',  action='store_false', dest='skipohmic', default=True, 
		help="Do NOT skip plots with negative d2J/dV2 values at vcutoff for Vtrans calcuation and conductance plot.")
parser.add_argument('-S','--smooth', type=float, default=0, 
		help="Cutoff for residuals when smoothing splines for dJ/dV 0 to disable. 1e-4 for artifacts, 1e-12 for smooth plots.")
parser.add_argument('-v','--vcutoff', type=float, default=-1, 
		help="Voltage (absolute value) cut-off for dJ/dV skipping routine (-1 for Vmin/Vmax)")
parser.add_argument('--Glower', metavar='GLOWER', dest='mlow', type=float, default=-6, 
		help="Lower cutoff value for conductance heat map plot.")
parser.add_argument('--Gupper', metavar='GUPPER', dest='mhi', type=float, default=0, 
		help="Upper cutoff value for conductance heat map plot.")
parser.add_argument('--NDClower', metavar='NDCLOWER', dest='ndc_mlow', type=float, default=0.05, 
		help="Lower cutoff value for normalized differential conductance heat map plot.")
parser.add_argument('--NDCupper', metavar='NDCUPPER', dest='ndc_mhi', type=float, default=5.0, 
		help="Upper cutoff value for conductance heat map plot.")
parser.add_argument('--heatmapbins', default=25, type=int, 
		help="Number of bins for the conductance heatmap plot.")
parser.add_argument('-R','--logr', default=True, action='store_false', 
		help="Compute |R| instead of log|R| for histograms of rectification.")
parser.add_argument('-L','--lorenzian', default=False, action='store_true', 
		help="Fit data to a Lorenzian instead of a Gaussian.")
parser.add_argument('-N','--nobs', type=int, default=0, 
		help="Number of observations for statistical tests on J (but not Gmean!).")
parser.add_argument('--maxfev', type=int, default=10000, 
		help="Maximum interations for fitting histograms.")
parser.add_argument('--plots', type=str, default='J', choices=['J','R'], 
		help="Log data to plot.")
parser.add_argument('--histplots', type=str, default='NDC', choices=['NDC','G'], 
		help="Derivative data to plot.")
parser.add_argument('-T','--threads', type=int,default=8, 
		help="Use n threads for parsing for marginal speed boost.")
parser.add_argument('--autonobs', default=False, action='store_true', 
		help="Try to find reasonable values of N automatically.")
parser.add_argument('--heatmapd',  type=int, default=1,
		     help="Derivative order of heatmap plot (0, 1, 2) default: 1.")
parser.add_argument('--tracebyfile', default=False, action='store_true', 
		help="Each input file contains one trace.")
parser.add_argument('--lagcutoff', default=0.1, type=float,
		help="Minimum euclidian distance from line to be considered as scatter.")



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
		
if not len(Opts.in_files) and not Opts.GUI and 0 in (len(Opts.setA),len(Opts.setB)):
	parser.print_help()
	print(RED+"\n\t\t> > > No input files! < < < "+RS)
	sys.exit()


# Setup CSV parser dialect and separator for pandas
if Opts.delim == 'tab':
	Opts.delim='\t'
if Opts.delim == 'comma':
	Opts.delim=','
if Opts.delim == 'space':
	Opts.delim=' '
csv.register_dialect('JV', delimiter=Opts.delim, quoting=csv.QUOTE_MINIMAL)
