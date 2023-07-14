from SLM.util import SLM_func, Gamma_func, slm_param_func
from gaussfit.parse.libparse.util import throwimportwarning

try:
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)


def doslm(self):
    '''
    Compute SLM curves for each trace that has valid params
    '''

    if not self.opts.SLM:
        return

    _eh_vals = []
    _gamma_vals = []
    _big_gamma_vals = []
    for trace in self.avg.index.levels[0]:
        _G = self.SLM["G"][trace]
        _Vtpos = self.SLM["Vtpos"][trace]
        _Vtneg = self.SLM["Vtneg"][trace]
        epsillon, gamma = slm_param_func(_Vtpos, _Vtneg)
        big_gamma = Gamma_func(_G, self.opts.nmolecules, _Vtpos, _Vtneg)
        if False in [abs(gamma) > 0, abs(epsillon) > 0, abs(big_gamma) > 0]:
            continue
        V = []
        J = []
        for _v in self.avg.loc[trace].index.tolist():
            if abs(_v) <= 1.25 * abs(epsillon):
                V.append(_v)  # Plotting past 1.25 x eh makes no sense
                J += self.avg.loc[trace]['J'][self.avg.loc[trace].index == _v].tolist()
        Y = []
        for _v in V:
            Y.append(SLM_func(_v, _G, epsillon, gamma))
        Y = SLM_func(np.array(V), _G, epsillon, gamma)
        self.SLM['epsillon'][trace] = epsillon
        self.SLM['gamma'][trace] = gamma
        self.SLM['big_gamma'][trace] = big_gamma
        self.SLM['exp'][trace] = [V, J]
        self.SLM['calc'][trace] = [V, Y]
        self.SLM['full'][trace] = [self.avg.loc[trace].index.tolist(), self.avg.loc[trace]['J'].tolist()]
        _eh_vals.append(epsillon)
        _gamma_vals.append(gamma)
        _big_gamma_vals.append(big_gamma)
    self.SLM['epsillon_avg'] = np.average(_eh_vals)
    self.SLM['gamma_avg'] = np.average(_gamma_vals)
    self.SLM['big_gamma_avg'] = np.average(_big_gamma_vals)
    V = [_x for _x in self.avg.loc[0].index.tolist() if _x <= 1.25 * abs(self.SLM['epsillon_avg'])]
    Y = []
    for _v in V:
        Y.append(SLM_func(_v, self.SLM['G_avg'], self.SLM['epsillon_avg'], self.SLM['gamma_avg']))
    self.SLM['calc_avg'] = [V, Y]