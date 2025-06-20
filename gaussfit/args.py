'''
Copyright (C) 2025 Ryan Chiechi <ryan.chiechi@ncsu.edu>
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
# pylint: disable=line-too-long

import sys
import os
import shutil
import warnings
import csv
import argparse
import configparser
from gaussfit.colors import RED, YELLOW, RS

warnings.filterwarnings('ignore', '.*divide by zero.*', RuntimeWarning)

VERSION = '1.3.0'

try:
    from scipy.optimize import curve_fit, OptimizeWarning  # pylint: disable=W0611
    from scipy.interpolate import UnivariateSpline  # pylint: disable=W0611
    import pandas as pd
    import numpy as np
    from appdirs import user_config_dir
    # SciPy throws a useless warning for noisy J/V traces
    warnings.filterwarnings('ignore', '.*Covariance of the parameters.*', OptimizeWarning)

except ImportError as msg:
    print("\n\t\t%s> > > Error importing package %s%s%s < < <%s" % (RED, RS, str(msg), RED, RS))
    sys.exit()


def doconfig(config_file):
    '''Parse config file or write a default file.'''
    if not os.path.exists(config_file):
        _pwd = os.path.dirname(os.path.realpath(__file__))
        _tdir = os.path.join(_pwd, '../templates/')
        shutil.copy2(os.path.join(_tdir, 'config.template'), config_file)
        print("Copied default config file.")
    _config = configparser.ConfigParser(allow_no_value=False)
    _config.read(config_file)
    return _config


def get_args():
    """
    Defines, parses, and processes all command-line arguments.

    This function isolates the ArgumentParser object, returning it alongside the
    final, processed options namespace. This ensures the returned 'opts' object
    is pickleable and safe for use with multiprocessing.

    Returns:
        (argparse.ArgumentParser, argparse.Namespace): A tuple containing the
        configured parser object and the parsed, validated options namespace.
    """
    _cachedir = user_config_dir(__package__)
    _configfile = os.path.join(_cachedir, __package__+'.conf')

    desc = '''
            This program expects all X values in one column and all Y values in another.
            Ideally, feed it *_data.txt files and it will take V and J. It can extract X- and
            Y-values from any two columns. Setting a compliance limit excludes Y > compliance
            from gaussian fits, but it rarely influences fits and can truncate data if set too low.
            The maxr parameter is needed to filter out huge values that occur when a junction shorts.
            '''

    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog='%sConfig file:%s%s' % (YELLOW, _configfile, RS))

    # =============================================================================
    # Positional and Core Arguments (Ungrouped)
    # =============================================================================
    parser.add_argument('in_files', metavar='Files-to-parse', type=str, nargs='*', default=[],
                        help='Datafiles to parse (GaussFit mode).')
    parser.add_argument('--loglevel', default='info', choices=('info', 'warn', 'error', 'debug'),
                        help="Set the logging level.")

    # =============================================================================
    # Input & Output Group
    # =============================================================================
    io_group = parser.add_argument_group('Input & Output')
    io_group.add_argument('-o', '--outfile', metavar="OUTPUT FILE", default="",
                        help="Output file name (taken from first input)")
    io_group.add_argument('--outdir', metavar="OUTPUT DIR", dest='out_dir', default=os.path.join(os.getcwd(), 'parsed'),
                        help="Output directory (combined with -o).")
    io_group.add_argument('-n', '--nowrite', dest='write', action='store_false', default=True,
                        help="Do not write output files (implies -p).")
    io_group.add_argument('--delim', default='tab', choices=('tab', 'comma', 'space'),
                        help="Delimeter in inputfiles.")
    io_group.add_argument('--encoding', default='ISO-8859-1', choices=('ISO-8859-1', 'utf-8', 'ascii'),
                        help="Encoding for input files.")
    io_group.add_argument('-X', '--xcol', type=int, default=1,
                        help="Column to treat as X (1-based index).")
    io_group.add_argument('-Y', '--ycol', type=int, default=3,
                        help="Column to treat as Y (1-based index). Set to 0 to grab all columns except Xcol.")

    # =============================================================================
    # Data Processing & Filtering Group
    # =============================================================================
    proc_group = parser.add_argument_group('Data Processing & Filtering')
    proc_group.add_argument('--maxr', type=float, default=3.0,
                        help="Maximum allowable value of log|R| or (R if -R). Filters shorted junctions.")
    proc_group.add_argument('--compliance', default=np.inf, type=float,
                        help="Set compliance limit for gaussian fits.")
    proc_group.add_argument('--smooth', type=float, default=0,
                        help="Cutoff for residuals when smoothing splines for dJ/dV. 0 to disable. 1e-4 for artifacts, 1e-12 for smooth plots.")
    proc_group.add_argument('--noskip', action='store_false', dest='skipohmic', default=True,
                        help="Do NOT skip plots with negative d2J/dV2 values at vcutoff for Vtrans calcuation and conductance plot.")
    proc_group.add_argument('--force', default=False, action='store_true',
                        help="Force GaussFit to attempt to plot/write bad data.")
    proc_group.add_argument('--vcutoff', type=float, default=-1,
                        help="Voltage (absolute value) cut-off for dJ/dV skipping routine (-1 for Vmin/Vmax)")
    proc_group.add_argument('--xrange', type=int, default=0,
                        help="Only parse x-axis to +/- X (useful for mistmatched volage ranges in inputs).")
    proc_group.add_argument('--tracebyfile', default=False, action='store_true',
                        help="Each input file contains one trace (e.g., for AFM data).")
    proc_group.add_argument('--lagcutoff', default=0.1, type=float,
                        help="Minimum euclidian distance from line to be considered as scatter.")
    proc_group.add_argument('--segments', type=int, default=4,
                        help="Number of segments in each J/V trace.")

    # ... (All other argument groups remain the same) ...

    # =============================================================================
    # Plotting & Visualization Group
    # =============================================================================
    plot_group = parser.add_argument_group('Plotting & Visualization')
    plot_group.add_argument('-p', '--plot', action='store_true', default=False,
                        help="Plot data and save to png file.")
    plot_group.add_argument('--plots', type=str, default='J', choices=['J', 'R', 'FN', 'SLM'],
                        help="Log data to plot.")
    plot_group.add_argument('--histplots', type=str, default='NDC', choices=['NDC', 'G'],
                        help="Heatmap data to plot.")
    plot_group.add_argument('--heatmapd', type=int, default=1,
                        help="Derivative order of heatmap plot (0, 1, 2) default: 1. 0 is equivalent to LogJ.")
    plot_group.add_argument('--bins', default=50, type=int,
                        help='Number of bins for histograms (except heatmap).')
    plot_group.add_argument('--heatmapbins', default=25, type=int,
                        help="Number of bins for the conductance heatmap plot.")
    plot_group.add_argument('--Glower', metavar='GLOWER', dest='mlow', type=float, default=-6,
                        help="Lower cutoff value for conductance heat map plot.")
    plot_group.add_argument('--Gupper', metavar='GUPPER', dest='mhi', type=float, default=0,
                        help="Upper cutoff value for conductance heat map plot.")
    plot_group.add_argument('--NDClower', metavar='NDCLOWER', dest='ndc_mlow', type=float, default=0.05,
                        help="Lower cutoff value for normalized differential conductance heat map plot.")
    plot_group.add_argument('--NDCupper', metavar='NDCUPPER', dest='ndc_mhi', type=float, default=5.0,
                        help="Upper cutoff value for conductance heat map plot.")
    plot_group.add_argument('--nologr', dest='logr', default=True, action='store_false',
                        help="Compute |R| instead of log|R| for histograms of rectification.")
    plot_group.add_argument('--nolag', action='store_true', default=False,
                        help="Do NOT compute lag plots (can be slow for some data types).")

    # =============================================================================
    # Statistical & Fitting Parameters Group
    # =============================================================================
    fit_group = parser.add_argument_group('Statistical & Fitting Parameters')
    fit_group.add_argument('--lorenzian', default=False, action='store_true',
                        help="Fit data to a Lorenzian instead of a Gaussian.")
    fit_group.add_argument('--maxfev', type=int, default=10000,
                        help="Maximum interations for fitting histograms.")
    fit_group.add_argument('--nobs', type=int, default=0,
                        help="Number of observations for statistical tests on J (but not Gmean!).")
    # fit_group.add_argument('--autonobs', default=False, action='store_true',
    #                     help="Try to find reasonable values of N automatically.")
    fit_group.add_argument('--alpha', type=float, default=0.025,
                        help="Alpha value for confidence intervals (p-cutoff = 1-Alpha).")
    fit_group.add_argument('--degfree', type=int, default=1,
                        help="Degrees of freedom. Set to 1 to infer from number of input files.")
    fit_group.add_argument('--minr', type=float, default=0.75,
                        help="Smallest R2-value to tolerate for linear fits.")
    fit_group.add_argument('--maxG', type=float, default=100,
                        help="Largest G-value to consider physically realistic.")
    fit_group.add_argument('--oldfn', action='store_true', default=False,
                        help="Use the old FN function.")
    # fit_group.add_argument('--interpolateminfn', action='store_true', default=False,
    #                     help="DEPRECATED: Compute Vtrans from the minimum of the derivative of the FN plot.")

    # =============================================================================
    # SLM & Clustering Group
    # =============================================================================
    adv_group = parser.add_argument_group('SLM & Clustering')
    adv_group.add_argument('--SLM', action='store_true', default=False,
                        help="Compute SLM fits of each trace.")
    adv_group.add_argument('--nmolecules', type=float, default=2.45e6,
                        help="Number of molecules per-junction to use for SLM calculations.")
    adv_group.add_argument('--cluster', action='store_true', default=False,
                        help="Perform clustering analysis on J/V curves.")
    adv_group.add_argument('--cluster-method', dest='cluster_estimation_method',
                        default='elbow', choices=['elbow', 'silhouette', 'gap'],
                        help="Method for determining optimal number of clusters.")
    adv_group.add_argument('--cluster-resolution', dest='cluster_resolution', type=int, default=100,
                        help="Resolution for clustering feature space. Higher values create more detailed features but require more memory/computation.")
    adv_group.add_argument('--cluster-as-egain', dest='cluster_as_egain', action='store_true', default=False,
                        help="When clustering, use complete EGaIn traces instead of individual voltage sweeps.")


    # =============================================================================
    # Config File Parsing and Post-Processing
    # =============================================================================
    config = doconfig(_configfile)

    _defaults = {}
    try:
        if 'BOOLEANS' in config:
            for _key in config['BOOLEANS']:
                _defaults[_key] = config['BOOLEANS'].getboolean(_key)
        if 'INTEGERS' in config:
            for _key in config['INTEGERS']:
                _defaults[_key] = config['INTEGERS'].getint(_key)
        if 'FLOATS' in config:
            for _key in config['FLOATS']:
                _defaults[_key] = config['FLOATS'].getfloat(_key)
        if 'STRINGS' in config:
            for _key in config['STRINGS']:
                _defaults[_key] = config['STRINGS'][_key]
    except KeyError as msg:
        print("%s%s is malformatted (missing: %s).%s" % (
            RED, _configfile, str(msg), RS))

    parser.set_defaults(**_defaults)
    opts = parser.parse_args()

    # The parser object is NOT attached to opts.
    # Attach configfile path, which is a simple string and safe to pickle.
    setattr(opts, 'configfile', _configfile)

    if not opts.write:
        opts.plot = True

    if len(opts.in_files) and not opts.outfile:
        if opts.in_files[0][-4] == '.':
            opts.outfile = os.path.basename(opts.in_files[0])[:-4]
        else:
            opts.outfile = opts.in_files[0]

    opts.slm_dir = os.path.join(opts.out_dir, 'SLM')
    opts.cluster_dir = os.path.join(opts.out_dir, 'Clustering')

    if opts.xcol == opts.ycol:
        print(RED+"Xcol and Ycol must be different."+RS)
        sys.exit()

    if opts.xcol < 1:
        print(RED+"X column must be greater than 0."+RS)
        sys.exit()

    # TODO: Parsing all columns is broken.
    if opts.ycol < 0:
        print(RED+"Y column must be greater than 0."+RS)
        sys.exit()

    # Convert to 0-based index for internal use
    opts.xcol -= 1
    opts.ycol -= 1

    if opts.degfree < 2 and len(opts.in_files) > 1:
        opts.degfree = len(opts.in_files)-1

    if not 0 < opts.alpha < 1:
        print(RED+"Alpha must be between 0 and 1"+RS)
        sys.exit()

    # Setup CSV parser dialect and separator for pandas
    if opts.delim == 'tab':
        opts.delim = '\t'
    elif opts.delim == 'comma':
        opts.delim = ','
    elif opts.delim == 'space':
        opts.delim = ' '
    csv.register_dialect('JV', delimiter=opts.delim, quoting=csv.QUOTE_MINIMAL)

    return parser, opts


# #!/usr/bin/env python3
# '''
# Copyright (C) 2025 Ryan Chiechi <ryan.chiechi@ncsu.edu>
# Description:
#         Parse all of the command line arguments for GaussFit.py
#         and Stats.py
# 
#         This program is free software: you can redistribute it and/or modify
#         it under the terms of the GNU General Public License as published by
#         the Free Software Foundation, either version 3 of the License, or
#         (at your option) any later version.
# 
#         This program is distributed in the hope that it will be useful,
#         but WITHOUT ANY WARRANTY; without even the implied warranty of
#         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#         GNU General Public License for more details.
# 
#         You should have received a copy of the GNU General Public License
#         along with this program.  If not, see <http://www.gnu.org/licenses/>.
# '''
# # pylint: disable=line-too-long
# 
# import sys
# import os
# import shutil
# import warnings
# import csv
# import argparse
# import configparser
# from gaussfit.colors import RED, YELLOW, RS
# 
# warnings.filterwarnings('ignore', '.*divide by zero.*', RuntimeWarning)
# 
# VERSION = '1.2.0'
# 
# try:
#     from scipy.optimize import curve_fit, OptimizeWarning  # pylint: disable=W0611
#     from scipy.interpolate import UnivariateSpline  # pylint: disable=W0611
#     import pandas as pd
#     import numpy as np
#     from appdirs import user_config_dir
#     # SciPy throws a useless warning for noisy J/V traces
#     warnings.filterwarnings('ignore', '.*Covariance of the parameters.*', OptimizeWarning)
# 
# except ImportError as msg:
#     print("\n\t\t%s> > > Error importing package %s%s%s < < <%s" % (RED, RS, str(msg), RED, RS))
#     sys.exit()
# 
# 
# def doconfig(config_file):
#     '''Parse config file or write a default file.'''
#     if not os.path.exists(config_file):
#         _pwd = os.path.dirname(os.path.realpath(__file__))
#         _tdir = os.path.join(_pwd, '../templates/')
#         shutil.copy2(os.path.join(_tdir, 'config.template'), config_file)
#         print("Copied default config file.")
#     _config = configparser.ConfigParser(allow_no_value=False)
#     _config.read(config_file)
#     # if 'BOOLEANS' not in _config:
#     #     _config = configparser.ConfigParser(allow_no_value=False)
#     #     _config.read(os.path.basename(config_file))
#     return _config
# 
# 
# _cachedir = user_config_dir(__package__)
# _configfile = os.path.join(_cachedir, __package__+'.conf')
# 
# desc = '''
#         This program expects all X values in one column and all Y values in another.
#         Ideally, feed it *_data.txt files and it will take V and J. It can extract X- and
#         Y-values from any two columns. Setting a compliance limit excludes Y > compliance
#         from gaussian fits, but it rarely influences fits and can truncate data if set too low.
#         The maxr parameter is needed to filter out huge values that occur when a junction shorts.
#        '''
# 
# parser = argparse.ArgumentParser(description=desc,
#                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#                                  epilog='%sConfig file:%s%s' % (YELLOW, _configfile, RS))
# 
# parser.add_argument('in_files', metavar='Files-to-parse', type=str, nargs='*', default=[],
#                     help='Datafiles to parse (GaussFit mode).')
# # parser.add_argument('-G', '--gui', action='store_true', default=False,
# #                     help="Launch the GUI.")
# parser.add_argument('-o', '--outfile', metavar="OUTPUT FILE", default="",
#                     help="Output file name (taken from first input)")
# parser.add_argument('--outdir', metavar="OUTPUT DIR", dest='out_dir', default=os.path.join(os.getcwd(), 'parsed'),
#                     help="Output directory (combined with -o).")
# parser.add_argument('--loglevel', default='info', choices=('info', 'warn', 'error', 'debug'),
#                     help="Set the logging level.")
# parser.add_argument('-p', '--plot', action='store_true', default=False,
#                     help="Plot data and save to png file.")
# parser.add_argument('-n', '--nowrite', dest='write', action='store_false', default=True,
#                     help="Do not write output files (implies -p).")
# parser.add_argument('--delim', default='tab', choices=('tab', 'comma', 'space'),
#                     help="Delimeter in inputfiles.")
# parser.add_argument('--encoding', default='ISO-8859-1', choices=('ISO-8859-1', 'utf-8', 'ascii'),
#                     help="Delimeter in inputfiles.")
# parser.add_argument('-X', '--xcol', type=int, default=1,
#                     help="Column to treat as X.")
# parser.add_argument('-Y', '--ycol', type=int, default=3,
#                     help="Column to treat as Y. Set to 0 to grab all columns except Xcol.")
# parser.add_argument('--xrange', type=int, default=0,
#                     help="Only parse x-axis to +/- X (useful for mistmatched volage ranges in inputs).")
# parser.add_argument('--bins', default=50, type=int,
#                     help='Number of bins for histograms (except heatmap).')
# parser.add_argument('--maxr', type=float, default=3.0,
#                     help="Maximum allowable value of log|R| or (R if -R).")
# parser.add_argument('--compliance', default=np.inf, type=float,
#                     help="Set compliance limit for gaussian fits.")
# parser.add_argument('--interpolateminfn', action='store_true', default=False,
#                     help="DEPRECATED: Compute Vtrans from the minimum of the derivative of the FN plot.")
# parser.add_argument('--oldfn', action='store_true', default=False,
#                     help="Use the old FN function.")
# parser.add_argument('--noskip', action='store_false', dest='skipohmic', default=True,
#                     help="Do NOT skip plots with negative d2J/dV2 values at vcutoff for Vtrans calcuation and conductance plot.")
# parser.add_argument('--smooth', type=float, default=0,
#                     help="Cutoff for residuals when smoothing splines for dJ/dV 0 to disable. 1e-4 for artifacts, 1e-12 for smooth plots.")
# parser.add_argument('--vcutoff', type=float, default=-1,
#                     help="Voltage (absolute value) cut-off for dJ/dV skipping routine (-1 for Vmin/Vmax)")
# parser.add_argument('--Glower', metavar='GLOWER', dest='mlow', type=float, default=-6,
#                     help="Lower cutoff value for conductance heat map plot.")
# parser.add_argument('--Gupper', metavar='GUPPER', dest='mhi', type=float, default=0,
#                     help="Upper cutoff value for conductance heat map plot.")
# parser.add_argument('--NDClower', metavar='NDCLOWER', dest='ndc_mlow', type=float, default=0.05,
#                     help="Lower cutoff value for normalized differential conductance heat map plot.")
# parser.add_argument('--NDCupper', metavar='NDCUPPER', dest='ndc_mhi', type=float, default=5.0,
#                     help="Upper cutoff value for conductance heat map plot.")
# parser.add_argument('--heatmapbins', default=25, type=int,
#                     help="Number of bins for the conductance heatmap plot.")
# parser.add_argument('--nologr', dest='logr', default=True, action='store_false',
#                     help="Compute |R| instead of log|R| for histograms of rectification.")
# parser.add_argument('--lorenzian', default=False, action='store_true',
#                     help="Fit data to a Lorenzian instead of a Gaussian.")
# parser.add_argument('--nolag', action='store_true', default=False,
#                     help="Do NOT compute lag plots.")
# parser.add_argument('--nobs', type=int, default=0,
#                     help="Number of observations for statistical tests on J (but not Gmean!).")
# parser.add_argument('--maxfev', type=int, default=10000,
#                     help="Maximum interations for fitting histograms.")
# parser.add_argument('--plots', type=str, default='J', choices=['J', 'R', 'FN', 'SLM'],
#                     help="Log data to plot.")
# parser.add_argument('--histplots', type=str, default='NDC', choices=['NDC', 'G'],
#                     help="Heatmap data to plot.")
# parser.add_argument('--threads', type=int, default=8,
#                     help="Use n threads for parsing for marginal speed boost.")
# parser.add_argument('--autonobs', default=False, action='store_true',
#                     help="Try to find reasonable values of N automatically.")
# parser.add_argument('--heatmapd', type=int, default=1,
#                     help="Derivative order of heatmap plot (0, 1, 2) default: 1. 0 is equivalent to LogJ.")
# parser.add_argument('--tracebyfile', default=False, action='store_true',
#                     help="Each input file contains one trace.")
# parser.add_argument('--lagcutoff', default=0.1, type=float,
#                     help="Minimum euclidian distance from line to be considered as scatter.")
# parser.add_argument('--alpha', type=float, default=0.025,
#                     help="Alpha value to use for computing confidence intervals (p-cutoff = 1-Alpha).")
# parser.add_argument('--segments', type=int, default=4,
#                     help="Number of segments in each J/V trace.")
# parser.add_argument('--degfree', type=int, default=1,
#                     help="Number of degrees of freedom (useful with --ycol=0). Set to 1 to infer from number of input files.")
# parser.add_argument('--minr', type=int, default=0.75,
#                     help="Smallest R2-value to tolerate for linear fits.")
# parser.add_argument('--maxG', type=int, default=100,
#                     help="Largest G-value to consider physically realistic.")
# parser.add_argument('--SLM', action='store_true', default=False,
#                     help="Compute SLM fits of each trace.")
# parser.add_argument('--nmolecules', type=int, default=2.45e6,
#                     help="Number of molecules per-junction to use for SLM calculations.")
# parser.add_argument('--cluster', action='store_true', default=False,
#                     help="Perform clustering analysis on J/V curves.")
# parser.add_argument('--cluster-method', dest='cluster_estimation_method', 
#                     default='elbow', choices=['elbow', 'silhouette', 'gap'],
#                     help="Method for determining optimal number of clusters.")
# parser.add_argument('--cluster-resolution', dest='cluster_resolution', type=int, default=100,
#                     help="Resolution for clustering feature space. Higher values create more detailed features but require more memory/computation. Default 100 (~10k features), modern systems can handle 200+ (40k+ features).")
# parser.add_argument('--cluster-as-egain', dest='cluster_as_egain', action='store_true', default=False,
#                     help="When clustering, use complete EGaIn traces (0→min→0→max→0) instead of individual voltage sweeps. Requires --cluster.")
# parser.add_argument('--force', default=False, action='store_true',
#                     help="Force GaussFit to attempt to plot/write bad data.")
# 
# _cachedir = user_config_dir(__package__)
# if not os.path.exists(_cachedir):
#     os.makedirs(_cachedir)
# _configfile = os.path.join(_cachedir, __package__+'.conf')
# config = doconfig(_configfile)
# 
# _defaults = {}
# try:
#     for _key in config['BOOLEANS']:
#         _defaults[_key] = config['BOOLEANS'].getboolean(_key)
#     for _key in config['INTEGERS']:
#         _defaults[_key] = config['INTEGERS'].getint(_key)
#     for _key in config['FLOATS']:
#         _defaults[_key] = config['FLOATS'].getfloat(_key)
#     for _key in config['STRINGS']:
#         _defaults[_key] = config['STRINGS'][_key]
# except KeyError as msg:
#     print("%s%s is malformatted (missing: %s).%s" % (
#         RED, _configfile, str(msg), RS))
# 
# parser.set_defaults(**_defaults)
# Opts = parser.parse_args()
# 
# setattr(Opts, 'configfile', _configfile)
# 
# if not Opts.write:
#     Opts.plot = True
# 
# if len(Opts.in_files) and not Opts.outfile:
#     if Opts.in_files[0][-4] == '.':
#         Opts.outfile = os.path.basename(Opts.in_files[0])[:-4]
#     else:
#         Opts.outfile = Opts.in_files[0]
# 
# Opts.slm_dir = os.path.join(Opts.out_dir, 'SLM')
# Opts.cluster_dir = os.path.join(Opts.out_dir, 'Clustering')
# 
# if Opts.xcol == Opts.ycol:
#     print(RED+"Xcol and Ycol must be different."+RS)
#     sys.exit()
# 
# if Opts.xcol < 1:
#     print(RED+"Xcolum must be greater than 0."+RS)
#     sys.exit()
# 
# # TODO: Parsing all columns is broken.
# if Opts.ycol < 0:
#     print(RED+"Ycolum must be greater than 0."+RS)
#     sys.exit()
# 
# Opts.xcol -= 1
# Opts.ycol -= 1
# 
# if Opts.degfree < 2 and len(Opts.in_files) > 1:
#     Opts.degfree = len(Opts.in_files)-1
# 
# if not 0 < Opts.alpha < 1:
#     print(RED+"Alpha must be between 0 and 1"+RS)
#     sys.exit()
# 
# 
# # Setup CSV parser dialect and separator for pandas
# if Opts.delim == 'tab':
#     Opts.delim = '\t'
# if Opts.delim == 'comma':
#     Opts.delim = ','
# if Opts.delim == 'space':
#     Opts.delim = ' '
# csv.register_dialect('JV', delimiter=Opts.delim, quoting=csv.QUOTE_MINIMAL)
