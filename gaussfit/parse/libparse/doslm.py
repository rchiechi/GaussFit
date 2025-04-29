from SLM.util import SLM_func, Gamma_func, slm_param_func
from gaussfit.parse.libparse.util import throwimportwarning
from scipy.special import stdtrit
from scipy.stats import gmean

try:
    import numpy as np
except ImportError as msg:
    throwimportwarning(msg)


def doslm(self):
    '''
    Compute SLM curves for each trace that has valid params
    '''

    fits = 0
    if not self.opts.SLM:
        return fits

    _eh_vals = []
    _gamma_vals = []
    _big_gamma_vals = []
    _g_vals = []
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
        if len(J) < len(self.avg.loc[trace].index.tolist()) / 3:
            self.logger.warning(f"Tossing SLM fit with only {len(J)} values.")
            continue
        Y = SLM_func(np.array(V), _G, epsillon, gamma)
        self.SLM['epsillon'][trace] = epsillon
        self.SLM['gamma'][trace] = gamma
        self.SLM['big_gamma'][trace] = big_gamma
        self.SLM['exp'][trace] = [V, J]
        self.SLM['calc'][trace] = [V, Y]
        self.SLM['full'][trace] = [self.avg.loc[trace].index.tolist(), self.avg.loc[trace]['J'].tolist()]
        _eh_vals.append(epsillon)
        _gamma_vals.append(gamma)
        _g_vals.append(_G)
        _big_gamma_vals.append(big_gamma)
        fits += 1
    self.SLM['G_avg'] = gmean(_g_vals, nan_policy='omit')
    self.SLM['epsillon_avg'] = gmean(_eh_vals, nan_policy='omit')
    self.SLM['gamma_avg'] = np.mean(np.array(_gamma_vals))  # no gmean on negative numbers
    self.SLM['big_gamma_avg'] = gmean(_big_gamma_vals, nan_policy='omit')
    try:
        self.SLM['epsillon_std'] = np.nanstd(_eh_vals)
        self.SLM['gamma_std'] = np.nanstd(_gamma_vals)
        self.SLM['big_gamma_std'] = np.nanstd(_big_gamma_vals)
        self.SLM['G_std'] = np.nanstd(_g_vals)
    except RuntimeWarning:
        self.logger.warning("Failed to compute any valid SLM values.")
        self.opts.SLM = False
        return
    self.SLM['epsillon_sem'] = self.SLM['epsillon_std'] / np.sqrt(fits - 1 or 1)
    self.SLM['epsillon_ci'] = self.SLM['epsillon_sem'] * stdtrit(fits - 1 or 1, 1 - self.opts.alpha)
    self.SLM['gamma_sem'] = self.SLM['gamma_std'] / np.sqrt(fits - 1 or 1)
    self.SLM['gamma_ci'] = self.SLM['gamma_sem'] * stdtrit(fits - 1 or 1, 1 - self.opts.alpha)
    self.SLM['big_gamma_sem'] = self.SLM['big_gamma_std'] / np.sqrt(fits - 1 or 1)
    self.SLM['big_gamma_ci'] = self.SLM['big_gamma_sem'] * stdtrit(fits - 1 or 1, 1 - self.opts.alpha)
    self.SLM['G_sem'] = self.SLM['G_std'] / np.sqrt(fits - 1 or 1)
    self.SLM['G_ci'] = self.SLM['G_sem'] * stdtrit(fits - 1 or 1, 1 - self.opts.alpha)

    V = [_x for _x in self.avg.loc[0].index.tolist() if _x <= 1.25 * abs(self.SLM['epsillon_avg'])]
    V = np.array(self.avg.loc[0].index.tolist())
    self.SLM['calc_avg'] = [V, SLM_func(V, self.SLM['G_avg'], self.SLM['epsillon_avg'], self.SLM['gamma_avg'])]
    return fits
