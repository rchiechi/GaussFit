#!/usr/bin/env python3
'''
Version: 1.0
Copyright (C) 2016 Ryan Chiechi <r.c.chiechi@rug.nl>
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
from gaussfit import *
from gaussfit.Output import Writer,Plotter
#from gaussfit.args import Opts
#from gaussfit.parse import Parse

#TODO colorama

def Go(opts):
    '''
    Call this function to execute the parsing engine
    i.e., as main()
    '''
    parser = Parse(opts)
    parser.ReadFiles(opts.in_files)
    if opts.write and not parser.error: 
            writer = Writer(parser)
            Parse.doOutput(writer)
    if opts.plot and not parser.error:
            plotter = Plotter(parser)
            try:
                    import matplotlib.pyplot as plt
                    plotter.DoPlots(plt,*opts.plots)
                    plt.show()
            except ImportError as msg:
                    print("Cannot import matplotlib! %s", str(msg), exc_info=False)

opts = Opts
if opts.GUI:
        from gui import filebrowser
        gui = filebrowser.ChooseFiles(opts)
else:
        Go(opts)
