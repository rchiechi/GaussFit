#!/usr/bin/env python3
'''
Version: 1.0
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

import logging
from gaussfit.Output import Writer,Plotter
#from gaussfit.Parseopts import Opts,ShowUsage 
from gaussfit.Args import Opts
from gaussfit.Parser import Parse

def Go(opts):
	'''
	Call this function to execute the parsing engine
	i.e., as main()
	'''
	parser = Parse(opts)
	parser.ReadFiles(opts.in_files)
	if opts.write:	
			writer = Writer(parser)
			logging.info("Writing files...")
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
	parser.PrintFN()
	if opts.plot:
			plotter = Plotter(parser)
			logging.info("Generating plots...")
			try:
					import matplotlib.pyplot as plt
					plotter.DoPlots(plt)
					plt.show()
			except ImportError as msg:
					logging.error("Cannot import matplotlib! %s", str(msg), exc_info=False)

opts = Opts
if opts.GUI:
		from gui import filebrowser
		gui = filebrowser.ChooseFiles(opts,Go)
else:
		Go(opts)
