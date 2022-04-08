#!/usr/bin/env python3
'''
Version: 1.0
Copyright (C) 2022 Ryan Chiechi <ryan.chiechi@ncsu.edu>
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
import cProfile
from gaussfit import Parse
from gaussfit.args import Opts
from gaussfit.output import Writer, Plotter
from gaussfit import colors
from gaussfit.output.libwriter import doOutput

try:
    import matplotlib.pyplot as plt
    CAN_PLOT = True
except ImportError as msg:
    print("Cannot import matplotlib! %s", str(msg), exc_info=False)
    CAN_PLOT = False


def do_cprofile(func):
    '''
    Needed for profiling GaussFit performance
    '''
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func


# @do_cprofile
def do_gaussfit():
    '''
    Call this function to execute the parsing engine
    i.e., as main()
    '''
    parser = Parse(Opts)
    parser.readfiles(Opts.in_files)
    if Opts.write and not parser.error:
        writer = Writer(parser)
        doOutput(writer)
    if Opts.plot and not parser.error:
        if CAN_PLOT:
            plotter = Plotter(parser, plt)
            plotter.DoPlots()
            plt.show()
        else:
            print(colors.RED+"Unable to show plots without matplotlib."+colors.RS)


if __name__ == "__main__":
    if Opts.gui:
        from gui import filebrowser
        gui = filebrowser.ChooseFiles()
    else:
        do_gaussfit()
