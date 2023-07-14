import os
import csv
import logging
import numpy as np
from SLM.fit import dofit

logger = logging.getLogger('output')


def WriteSLM(self):
    '''Write the SLM inputs and parameters per-trace.'''

    if not self.opts.SLM:
        return

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_SLM_inputs.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["Trace", "G", "Vtpos", "Vtneg"])
        traces = list(self.SLM['G'].keys())
        for trace in traces:
            writer.writerow([f'{trace}',
                             f'{self.SLM["G"][trace]}',
                             f'{self.SLM["Vtpos"][trace]}',
                             f'{self.SLM["Vtneg"][trace]}'])

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_SLM_params.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        writer.writerow(["G", "ε", "γ", "Γ"])
        traces = list(self.SLM['epsillon'].keys())
        for trace in traces:
            writer.writerow([f'{self.SLM["G"][trace]}',
                             f'{self.SLM["epsillon"][trace]}',
                             f'{self.SLM["gamma"][trace]}',
                             f'{self.SLM["big_gamma"][trace]}'])

    _fn = os.path.join(self.opts.out_dir, self.opts.outfile + "_SLM_calc_avg.txt")
    with open(_fn, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='JV')
        _G = self.SLM['G_avg']
        _eh = self.SLM['epsillon_avg']
        _gamma = self.SLM['gamma_avg']
        _big_gamma = self.SLM['big_gamma_avg']
        writer.writerow(["V", f"J_calc (G={_G:0.4E} ε={_eh:0.4f} γ={_gamma:0.4f} Γ={_big_gamma:0.4f})"])
        _v, _j = self.SLM['calc_avg']
        for _i in range(len(_v)):
            writer.writerow([f'{_v[_i]}', f'{_j[_i]}'])

    traces = list(self.SLM['calc'].keys())
    for trace in traces:
        _fn = os.path.join(self.opts.slm_dir, f"{self.opts.outfile}_SLM_calc_trace{trace:03}.txt")
        with open(_fn, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='JV')
            _G = self.SLM['G'][trace]
            _eh = self.SLM['epsillon'][trace]
            _gamma = self.SLM['gamma'][trace]
            _big_gamma = self.SLM['big_gamma'][trace]
            writer.writerow(["V", "J", f"J_calc (G={_G:0.4E} ε={_eh:0.4f} γ={_gamma:0.4f} Γ={_big_gamma:0.4f})"])
            _v, _j = self.SLM['calc'][trace]
            for _i in range(len(_v)):
                writer.writerow([f'{_v[_i]}', f'{self.SLM["exp"][trace][1][_i]}',
                                 f'{_j[_i]}',
                                 ])
        _v, _j = self.SLM['full'][trace]
        _fn = f"{self.opts.outfile}_SLM_calc_trace{trace:03}.txt"
        dofit(np.array(_v), np.array(_j), p0=[_G, _eh, _gamma],
              path=self.opts.slm_dir,
              save=_fn)
        _v, _j = self.SLM['full'][trace]
        _fn = f"{self.opts.outfile}_SLM_fit_trace{trace:03}.txt"
        dofit(np.array(_v), np.array(_j),
              path=self.opts.slm_dir,
              save=_fn)
