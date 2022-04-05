#! /usr/bin/env python3
"""A python template for plotting JV curves generated by GaussFit.
   It expects two directories, each containing files to be plotted together."""

import os
import sys
import csv
import numpy as np
from matplotlib import pyplot as plt

######################################################################
# If you want to plot R or NDC or whatever, change these labels
YLABEL = r'$log|J A cm^{-2}|$'
XLABEL = r'Potential (V)'
######################################################################


def parse(in_file):
    data = {'X': [], 'Y': [], 'Yerror': [], 'Yupper': [], 'Ylower': []}
    with open(in_file, newline="") as csvfile:
        print("Parsing %s" % in_file)
        reader = csv.reader(csvfile, dialect='JV')
        for __line in reader:
            try:
                data['X'].append(float(__line[0]))
                data['Y'].append(float(__line[1]))
                data['Yerror'].append(float(__line[2]))
                data['Yupper'].append(data['Y'][-1] + data['Yerror'][-1])
                data['Ylower'].append(data['Y'][-1] - data['Yerror'][-1])

            except ValueError:
                print("Found labels: %s" % __line)

    for key in ('X', 'Y', 'Yerror'):
        data[key] = np.array(data[key])
    return data


csv.register_dialect('JV', delimiter='\t', quoting=csv.QUOTE_MINIMAL)
datasets_top = []
datasets_bottom = []

plotcolors = ('red',
              'indigo',
              'green',
              'cornflowerblue',
              'darkorange',
              'fuchsia',
              'blue'
              )
plotpoints = ('o', 's', 'd', 'v', '^', '>', '<')
label_style = {'size': 'x-large'}
plot_style = {'color': 'tab:purple', 'lw': 2.0}
major_tick_style = {'axis': 'both',
                    'which': 'major',
                    'labelsize': '14',
                    'left': True,
                    'right': True,
                    'top': True,
                    'bottom': True,
                    'direction': 'inout'}

for in_file in os.listdir(sys.argv[1]):
    datasets_bottom.append(parse(os.path.join(sys.argv[1], in_file)))
for in_file in os.listdir(sys.argv[2]):
    datasets_top.append(parse(os.path.join(sys.argv[2], in_file)))

fig = plt.figure()
gs = fig.add_gridspec(2, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)

i = 0
for data in datasets_bottom:
    axs[0].plot('X', 'Y', data=data, marker=plotpoints[i], color=plotcolors[i], lw=1.5)
    axs[0].fill_between('X', 'Ylower', 'Yupper', data=data, interpolate=True, color=plotcolors[i], alpha=0.1)
    i += 1

i = 0
for data in datasets_top:
    axs[1].plot('X', 'Y', data=data, marker=plotpoints[i], color=plotcolors[i], lw=1.5)
    axs[1].fill_between('X', 'Ylower', 'Yupper', data=data, interpolate=True, color=plotcolors[i], alpha=0.1)
    i += 1

for ax in axs:
    ax.label_outer()
    ax.legend(loc=3)
    ax.set_ylabel(YLABEL, **label_style)
    ax.set_xlabel(XLABEL, **label_style)
    ax.tick_params(**major_tick_style)


plt.rcParams["figure.figsize"] = [10.0, 6.0]
plt.tight_layout()
plt.show()
