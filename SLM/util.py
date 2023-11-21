import numpy as np
import scipy.constants as constants


G0 = (2 * np.power(constants.e, 2)) / constants.h


def SLM_func(V: float, G: float, eh: float, y: float) -> float:
    _a = G * V * np.power(eh, 2)
    _b1 = np.power(eh - (y * V), 2)
    _b2 = np.power(V / 2, 2)
    return _a / (_b1 - _b2)


def slm_param_func(Vtpos: float, Vtneg: float) -> (float, float):
    if False in (abs(Vtpos) > 0, abs(Vtneg) > 0):
        return np.nan, np.nan
    _sum = Vtpos + Vtneg
    _product = Vtpos * Vtneg
    _den = np.sqrt(((3 / 10) * abs(_product)) + np.power(Vtneg, 2) + np.power(Vtpos, 2))
    epsillon = (2 * abs(_sum)) / _den
    gamma = (-0.5 * _sum) / _den
    return epsillon, gamma


def Gamma_func(G: float, N: float, Vtpos: float, Vtneg: float) -> float:
    if False in (abs(G) > 0, abs(N) > 0, abs(Vtpos) > 0, abs(Vtneg) > 0):
        return np.nan
    epsillon = slm_param_func(Vtpos, Vtneg)[0]
    return epsillon * np.sqrt(G / (N * G0))
