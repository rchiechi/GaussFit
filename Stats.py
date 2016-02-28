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


from gaussfit.Output import WriteStats,StatPlotter
from gaussfit.StatArgs import Opts
#from gaussfit.Parser import Parse
#from gaussfit.Output import Writer

import sys,os,logging,warnings,csv
from gaussfit.colors import *
from gaussfit.stats import Stats

def Go(opts):
	'''
	Call this function to execute the parsing engine
	i.e., as main()
	'''

	if len(opts.setA) and len(opts.setB):
		statparser = Stats(opts)
	else:
		logging.error("No input files to parse!")

	if opts.plot:
			plotter = StatPlotter(statparser)
			logging.info("Generating plots...")
			try:
					import matplotlib.pyplot as plt
					plotter.DoPlots(plt)
					plt.show()
			except ImportError as msg:
					logging.error("Cannot import matplotlib! %s", str(msg), exc_info=False)


opts = Opts
if opts.GUI:
		from gui import statsfilebrowser
		gui = statsfilebrowser.ChooseFiles(opts,Go)
else:
		Go(opts)
