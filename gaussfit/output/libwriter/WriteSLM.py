import os
import csv
import logging

logger = logging.getLogger('output')


def WriteSLM(self):
    '''Write the transition voltage data plotted against potential (not 1/V).'''
    _fn = os.path.join(self.opts.out_dir, self.opts.outfile+"_SLM_inputs.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Trace", "G", "Vtpos", "Vtneg"])
        traces = list(self.SLM['G'].keys())
        for trace in traces:
            writer.writerow([f'{trace}',
                             f'{self.SLM["G"][trace]}',
                             f'{self.SLM["Vtpos"][trace]}',
                             f'{self.SLM["Vtneg"][trace]}'])
