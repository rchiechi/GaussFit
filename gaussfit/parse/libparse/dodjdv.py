import numpy as np
from collections import OrderedDict
import scipy.interpolate  # type: ignore
# from scipy.ndimage import gaussian_filter1d as gaussian_filter  # type: ignore
from gaussfit.parse.libparse.dohistogram import dohistogram  # type: ignore
from gaussfit.args import Opts

def dodjdv(self):
    '''
    Fit a spline function to X/Y data,
    compute dY/dX and normalize.
    '''
    self.logger.info("* * * * * * Computing dY/dX  * * * * * * * *")
    self.loghandler.flush()
    linx = np.linspace(self.df.V.min(), self.df.V.max(), 200)
    if Opts.vcutoff > 0:
        self.logger.debug('Using %s cutoff for dj/dv', Opts.vcutoff)
        vfilterneg, vfilterpos = np.linspace(-1*Opts.vcutoff, 0, 200), np.linspace(0, Opts.vcutoff.max(), 200)
    else:
        vfilterneg, vfilterpos = np.linspace(self.df.V.min(), 0, 200), np.linspace(0, self.df.V.max(), 200)
    if Opts.vcutoff > 0:
        vfilterneg, vfilterpos = linx[-1*Opts.vcutoff < linx < 0], linx[0 < linx < Opts.vcutoff]
    else:
        vfilterneg, vfilterpos = linx[linx < 0], linx[linx > 0]

    spls = OrderedDict()
    spls_norm = OrderedDict()
    splhists = OrderedDict()
    spl_normhists = OrderedDict()
    ndc_cut, ndc_tot = 0, 0
    filtered = [('Potential', 'Fit', 'Y')]
    for x in linx:
        spls[x] = []
        splhists[x] = {'spl': [], 'hist': {}}
        spls_norm[x] = []
        spl_normhists[x] = {'spl': [], 'hist': {}}
    for trace in self.avg.index.levels[0]:
        try:
            spl = scipy.interpolate.UnivariateSpline(
                self.avg.loc[trace].index, self.avg.loc[trace]['J'], k=5, s=Opts.smooth)
            dd = scipy.interpolate.UnivariateSpline(
                self.avg.loc[trace].index, self.avg.loc[trace]['J'], k=5, s=None).derivative(2)
        except Exception as msg:  # pylint: disable=broad-except
            # TODO: Figure out how to catch all the various scipy errors here
            self.logger.error('Error in derivative calulation: %s', str(msg))
            continue
        try:
            spldd = dd(vfilterpos)  # Compute d2J/dV2
            spldd += -1*dd(vfilterneg)  # Compute d2J/dV2
        except ValueError as msg:
            self.logger.error('Error computing second derivative: %s', str(msg))
            continue
        if len(spldd[spldd < 0]):
            # record in the index where dY/dX is < 0 within vcutoff range
            self.ohmic.append(trace)
            if Opts.skipohmic:
                continue
        else:
            for row in self.avg.loc[trace].iterrows():
                # filtered is a list containing only "clean" traces
                filtered.append((row[0], spl(row[0]), row[1].J))
        err = None
        for x in spls:
            try:
                d = spl.derivatives(x)
            except (TypeError, ValueError) as msg:
                err = str(msg)
                self.logger.warning('Error computing derivative: %s', str(msg))
                # print(f'Error computing derivative:{str(msg)}')
                continue
            # except TypeError:
            #     self.logger.warning('Applying Guassian filter to prevent spline errors.')
            #     d = spl.derivatives(gaussian_filter(spl, 5))
            if np.isnan(d[Opts.heatmapd]):
                self.logger.warning("Got NaN computing dJ/dV")
                # print("Got NaN computing dJ/dV")
                continue
            spls[x].append(d[Opts.heatmapd])
            splhists[x]['spl'].append(np.log10(abs(d[Opts.heatmapd])))
            ndc = d[1] * (x/spl(x))
            spls_norm[x].append(ndc)
            ndc_tot += 1
            if 0.0 < ndc < 10.0:
                # ndc values beyond this range can safely be considered artifacts of the numerical derivative
                spl_normhists[x]['spl'].append(ndc)
            else:
                # but we keep track of how many data points we toss as a sanity check
                ndc_cut += 1
        if err:
            self.logger.error("Error while computing derivative: %s", str(err))

    self.logger.info("Non-tunneling traces: %s (out of %0d)",
                     len(self.ohmic), len(self.avg.index.levels[0]))
    if (len(self.ohmic) == len(self.avg.index.levels[0])) and Opts.skipohmic:
        self.logger.error("You have elected to skip all traces: disable skip non-ohmic and re-parse!")
        self.DJDV, self.GHists, self.NDC, self.NDCHists, self.filtered = {}, {}, {}, {}, {}
        return
    if ndc_tot:
        self.logger.info("NDC values not between 0 and 10: %s (%0.2f%%)", ndc_cut, (ndc_cut/ndc_tot)*100)
    self.loghandler.flush()
    for x in splhists:
        splhists[x]['hist'] = dohistogram(np.array(splhists[x]['spl']), label='DJDV', que=self.logqueue)
        spl_normhists[x]['hist'] = dohistogram(np.array(spl_normhists[x]['spl']), label='NDC', que=self.logqueue)
    self.logger.info("dJdV complete.")
    self.loghandler.flush()
    self.DJDV, self.GHists, self.NDC, self.NDCHists, self.filtered = spls, splhists, spls_norm, spl_normhists, filtered
