#! /usr/bin/env python3
"""A python template for plotting JV curves generated by GaussFit"""

import sys
import csv
import numpy as np
from matplotlib import pyplot as plt

csv.register_dialect('JV', delimiter='\t', quoting=csv.QUOTE_MINIMAL)
datasets = []
for in_file in sys.argv[1:]:
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
    datasets.append(data)


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
                     'direction': 'inout'}




i = 0
for data in datasets:
    plt.plot('X', 'Y', data=data, marker=plotpoints[i], color=plotcolors[i], lw=2.0)
    plt.fill_between('X', 'Ylower', 'Yupper', data=data, interpolate=True, color=plotcolors[i], alpha=0.2)
    i += 1



ax = plt.gca()

ax.legend(loc=3)
ax.set_ylabel(r'$log|J A cm^{-2}|$', **label_style)
ax.set_xlabel(r'Potential (V)', **label_style)

plt.rcParams["figure.figsize"] = [10.0, 6.0]
plt.tick_params(**major_tick_style)
plt.tight_layout()
plt.show()
