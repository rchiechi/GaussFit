#!/usr/bin/env python3

import os
import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import csv
import argparse
from colorama import Fore, Style


delim_choices = {'comma':',', 'tab':'\t', 'semicolon':';'}
desc = '''
       Fit EGaIn data to the SLM.
       '''

parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 )

parser.add_argument('in_files', type=str, nargs='*', default=[],
                    help='Datafiles to parse.')
parser.add_argument('--unlog', action='store_true', default=False,
                    help='Input data are log10 and need to be un-logged.')
parser.add_argument('--correctsign', action='store_true', default=False,
                    help='Input data were un-logged without correcting for sign.')
parser.add_argument('-c', '--column', type=int, default=2,
                    help='Column holding current (y-axis) data.')
parser.add_argument('-d', '--delim', default='comma', choices=list(delim_choices.keys()),
                    help='CSV file delimeter.')
parser.add_argument('--params', type=str, default='',
                    help='A csv file with columns <filename> <G> <ε> <γ> matching the input file names.')
parser.add_argument('--plot', action='store_true', default=False,
                    help='Show previews of the fits.')
parser.add_argument('--save', action='store_true', default=False,
                    help='Save previews of the fits as png files.')
parser.add_argument('--scale', type=float, default=1,
                    help='Scale current (y) by this factor.')
parser.add_argument('--traces', action='store_true', default=False,
                    help='Input data are Traces files from GaussFit.')


def parseparams(fn, **kwargs):
    paramfiles = {}
    delim = kwargs.get('delim', delim_choices['comma'])
    paths = kwargs.get('paths', [])
    with open(fn, newline='') as fh:
        _paths = paths + [os.path.split(fn)[0]]
        _csv = csv.reader(fh, delimiter=delim)
        print(Fore.YELLOW, end='')
        for row in _csv:
            csv_fn = None
            _csv_fn, G, eh, y = row[0:4]
            for _path in _paths + [os.path.split(_csv_fn)[0]]:
                _fn = os.path.join(_path, os.path.basename(_csv_fn))
                if os.path.exists(_fn):
                    csv_fn = os.path.join(_fn)
                    break
                elif os.path.exists(f'{_fn}.csv'):
                    csv_fn = f'{_fn}.csv'
                    break
                else:
                    csv_fn = None

            if not _csv_fn or csv_fn is None:
                print(f'{_csv_fn} does not exist in {_paths}')
                continue
            else:
                print(f'{Fore.GREEN}Found:{csv_fn}{Fore.YELLOW}')
            try:
                _params = {'G': float(G),
                           'eh': float(eh),
                           'y': float(y)}
            except ValueError:
                print(f'{Fore.RED}Cannot parse row: {row}', end=f'{Fore.YELLOW}\r')
                continue
            paramfiles[os.path.abspath(csv_fn)] = _params
    print(Style.RESET_ALL)
    return paramfiles

def parse_traces(in_file, **kwargs):
    delim = kwargs.get('delim', delim_choices['comma'])
    scale = kwargs.get('scale', 1.0)
    traces = []
    with open(in_file, newline='') as fh:
        _csv = csv.reader(fh, delimiter=delim)
        _last_trace = None
        V = []
        J = []
        for row in _csv:
            try:
                n_trace, _V, _J, _FN = list(map(float, row))
            except ValueError:
                print(f"parse_traces: {Fore.RED}Could not unpack row {row}{Style.RESET_ALL}")
                continue
            if _last_trace is None:
                _last_trace = n_trace
            elif _last_trace != n_trace:
                _last_trace = n_trace
                traces.append([np.array(V),np.array(J)])
                J, V = [], []
            V.append(_V)
            J.append(_J * scale)
    return traces

def getpaths(in_files):
    paths = []
    for _fn in in_files:
        _p = os.path.split(_fn)[0]
        if _p not in paths:
            paths.append(_p)
    return paths

def SLM(V, G, eh, y):
    # I = GV (eh)2 / (eh - yeV)2 - (eV/2)2
    _a = G*V*np.power(eh, 2)
    _b1 = np.power(eh - (y*V), 2)
    _b2 = np.power(V/2, 2)
    return _a/(_b1-_b2)


def get_data(fn, y, **kwargs):
    V, I = [], []
    y = int(y)
    delim = kwargs.get('delim', delim_choices['comma'])
    unlog = kwargs.get('unlog', False)
    correctsign = kwargs.get('correctsign', False)
    scale = kwargs.get('scale', 1.0)
    # print(Fore.RED, end='')
    with open(fn, newline='') as fh:
        _csv = csv.reader(fh, delimiter=delim)
        for row in _csv:
            try:
                _V = (float(row[0]))
                _I = (float(row[y]))
                if unlog:
                    _I = np.power(10, _I)
                if correctsign or unlog:
                    if _V < 0 and _I > 0:
                        _I = -1 * _I
                    elif _V == 0:
                        _I = 0
                _I = _I * scale
            except ValueError:
                # print(f'No data in:{row}', end='\r')
                continue
            V.append(_V)
            I.append(_I)
    return np.array(V), np.array(I)


def dofit(x, y, **kwargs):
    p0 = kwargs.get('p0', [1e-9, 1, 0.01])
    bounds = kwargs.get('bounds', ([0, 0, -5], [1, 10, 5]))
    popt, pcov = curve_fit(SLM, x, y, p0=p0, bounds=bounds)

    print(f"G: {popt[0]} ε: {popt[1]} γ: {popt[2]}")
    if kwargs.get('plot', False) or kwargs.get('save', None) is not None:
        plt.plot(x, y, 'k', linewidth=2, label='EGaIn Data')
        plt.plot(x, SLM(x, *popt), 'g--',
                 label=f'Fit (G: {popt[0]:0.2E} ε: {popt[1]:0.2f} γ: {popt[2]:0.4f})')
        if p0 is not None:
            plt.plot(x, SLM(x, *p0), 'b--',
                     label=f'Guess (G: {p0[0]:0.2E} ε: {p0[1]:0.2f} γ: {p0[2]:0.4f})')
        if kwargs.get('traces'):
            plt.ylabel('Current Density (A $\mathrm{cm^{-2}}$)')
        else:
            plt.ylabel('Current (A)')
        plt.xlabel('Bias (V)')
        plt.legend()
        plt.tight_layout()
        if kwargs.get('plot', False):
            plt.show()
        if kwargs.get('save', None) is not None:
            plt.savefig(f"{kwargs['save'][0:-4]}.png")
        plt.clf()
    return popt, pcov

def write_output(fn, x, y, p0, popt, pcov, **kwargs):
    delim = kwargs.get('delim', delim_choices['comma'])
    if p0 is None:
        p0 = [0, 0, 0]
    guess = SLM(x, *p0)
    fit = SLM(x, *popt)
    print(f"{Fore.BLUE}{Style.BRIGHT}Writing {fn}{Style.RESET_ALL}")
    with open(fn, 'wt', newline='') as fh:
        _csv = csv.writer(fh, delimiter=delim)
        _csv.writerow(['Bias (V)',
                       'Current (I)',
                       f'Fit: G={popt[0]} ε={popt[1]} γ={popt[2]}',
                       f'Guess: G={p0[0]} ε={p0[1]} γ={p0[2]}']
                      )
        for i in range(len(x)):
            _csv.writerow([x[i], y[i], fit[i], guess[i]])

def parse_files(opts):
    to_parse = {}
    paramfiles = {}
    if opts.params:
        paramfiles = parseparams(os.path.abspath(opts.params),
                                 delim=delim_choices[opts.delim],
                                 paths=getpaths(opts.in_files))
    for _fn in opts.in_files:
        fn = os.path.abspath(_fn)
        if opts.traces:
            _traces = parse_traces(
                fn,
                delim=delim_choices[opts.delim],
                scale=opts.scale)
            for _n, _trace in enumerate(_traces):
                to_parse[f'{_n}_{fn}'] = {'data':_trace, 'p0':None}
                if fn in paramfiles:
                    paramfiles[f'{_n}_{fn}'] = paramfiles[fn]
        else:
            to_parse[fn] = {'data':get_data(
                            fn,
                            opts.column,
                            delim=delim_choices[opts.delim],
                            unlog=opts.unlog,
                            correctsign=opts.correctsign,
                            scale=opts.scale),
                            'p0':None}
        if fn in paramfiles:
            print(f"{Fore.GREEN}Found params for {os.path.basename(fn)}.{Style.RESET_ALL}")
            to_parse[fn]['p0'] = (paramfiles[fn]['G'],
                                  paramfiles[fn]['eh'],
                                  paramfiles[fn]['y'],)
        elif opts.params:
            print(paramfiles.keys())
            print(f"{Style.BRIGHT}{Fore.YELLOW}{os.path.basename(fn)} not found in params.{Style.RESET_ALL}")
    _save = None
    for fn in to_parse:
        if opts.save:
            _save = os.path.basename(fn)
        popt, pcov = dofit(to_parse[fn]['data'][0],
                           to_parse[fn]['data'][1],
                           p0=to_parse[fn]['p0'],
                           plot=opts.plot,
                           traces=opts.traces,
                           save=_save)
        if opts.save:
            write_output(f"{''.join(os.path.basename(fn).split('.')[:-1])}_fit.csv",
                         to_parse[fn]['data'][0],
                         to_parse[fn]['data'][1],
                         to_parse[fn]['p0'],
                         popt, pcov)


if __name__ == '__main__':
    opts = parser.parse_args()
    opts.column = opts.column - 1  # Use human indices not python for column
    if not len(opts.in_files):
        print(f"{Fore.RED}No input files?{Style.RESET_ALL}")
        parser.print_help()
        sys.exit()
    parse_files(opts)
