import os
import warnings
import logging
from gaussfit.colors import GREEN, TEAL, YELLOW, WHITE
from matplotlib.cm import get_cmap
# from gaussfit.args import get_args
from SLM.util import SLM_func


logger = logging.getLogger(__package__)
loghandler = logging.StreamHandler()
loghandler.setFormatter(logging.Formatter(
                fmt=GREEN+os.path.basename('%(name)s'+TEAL)+' %(levelname)s '+YELLOW+'%(message)s'+WHITE))
logger.addHandler(loghandler)

try:
    import numpy as np
    from scipy.interpolate import griddata
    from scipy.special import stdtrit
    # import matplotlib
    from matplotlib import container
except ImportError:
    pass  # We catch numpy import errors in Parser.py
warnings.filterwarnings('ignore', '.*comparison.*', FutureWarning)


class Plotter():
    '''
    This is the main Plotter class for generating
    plots using matplotlib.
    '''

    def __init__(self, parser, plt_or_fig):
        self.parser = parser
        self.opts = self.parser.opts
        __figure = getattr(plt_or_fig, 'figure', None)
        if callable(__figure):
            self.fig = plt_or_fig.figure(figsize=(16, 10))
        else:
            self.fig = plt_or_fig

    def __getattr__(self, name):
        try:
            return getattr(self.parser, name)
        except AttributeError as msg:
            raise AttributeError("Plotter object has no attribute '%s'" % name) from msg

    def PlotData(self, key, ax, sym, **kw):
        xax = np.array(list(self.XY))
        if key == "FN":
            xax = 1/xax
            ax.set_title("Fowler Nordheim Plot of Initial Data")
            ax.set_xlabel(r'$V^{-1}$')
            ax.set_ylabel(r'$\mathregular{ln(\frac{J}{V^2})}$')
        if key == 'Y':
            if self.opts.compliance != np.inf:
                ax.set_ylim((-1*self.opts.compliance, self.opts.compliance))
            ax.set_title("Initial Data")
            ax.set_xlabel("Potenial (V)")
            ax.set_ylabel(r'Current Density ($A cm^{-2}$)')
        if key == 'LogY':
            ax.set_title("Semilog Plot of Initial Data")
            ax.set_xlabel("Potenial (V)")
            ax.set_ylabel(r'Current Density $\mathregular{log_{10}|J(\mathrm{A cm^{-2}})|}$')
        i = -1
        allY = np.array([])
        while True:
            i += 1
            try:
                allY = np.append(allY, [self.XY[x][key].iloc[i] for x in self.XY])
                ax.plot(xax, [self.XY[x][key].iloc[i] for x in self.XY], sym, **kw)
            except IndexError:
                break
            except ValueError:
                break
        ax.axis([xax.min(), xax.max(), allY.min(), allY.max()])

    def PlotSLM(self, ax):
        ax.set_title("Fit to SLM")
        ax.set_ylabel(r'log|Current Density $J(\mathrm{A cm}^{-2})$|')
        ax.set_xlabel(r'Potential (V)')
        G = self.SLM["Gauss"]["G"]
        epsillon = self.SLM["Gauss"]["epsillon"]
        gamma = self.SLM["Gauss"]["gamma"]
        x = np.array([x for x in self.XY])
        ax.plot(x, np.log10(abs(SLM_func(x, G, epsillon, gamma))), lw=3.0, color='k', label='From GaussJ')
        traces = list(self.SLM['calc'].keys())
        for trace in traces:
            _v, _j = self.SLM['calc'][trace]
            ax.plot(_v, np.log10(abs(np.array(_j))), ':', lw=0.5, color='r')
        _v, _j = self.SLM['calc_avg']
        ax.plot(_v, np.log10(abs(np.array(_j))), lw=3.0, color='b', label='Gmean of per-Trace')
        ax.legend()

    def PlotClustering(self, ax):
        '''Plot clustering results'''
        if not hasattr(self.parser, 'cluster') or self.parser.cluster['clusterer'] is None:
            ax.text(0.5, 0.5, 'No clustering data available', 
                   horizontalalignment='center', verticalalignment='center', 
                   transform=ax.transAxes)
            ax.set_title("Clustering Analysis")
            return
        
        clusterer = self.parser.cluster['clusterer']
        jv_curves = self.parser.cluster['jv_curves']
        
        try:
            # Plot clustering results
            fig = clusterer.plot_results(jv_curves, figsize=(8, 6), log_scale=False)
            
            # Copy the plot to our axis (simplified version)
            ax.set_title(f"J/V Curve Clustering ({self.parser.cluster['n_clusters']} clusters)")
            ax.set_xlabel("Representative clusters shown in main analysis")
            ax.text(0.5, 0.5, f"Found {self.parser.cluster['n_clusters']} clusters\nSee separate clustering plots", 
                   horizontalalignment='center', verticalalignment='center', 
                   transform=ax.transAxes, fontsize=12, 
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        except Exception as e:
            ax.text(0.5, 0.5, f'Clustering plot error: {str(e)}', 
                   horizontalalignment='center', verticalalignment='center', 
                   transform=ax.transAxes)
            ax.set_title("Clustering Analysis - Error")

    def PlotSegmentedGauss(self, ax):
        '''Plot segmented J/V data'''
        ax.set_title("Semilog Plots of |J| by Trace")
        ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm}^{-2})|$')
        ax.set_xlabel(r'Potential (V)')
        X, Y = {}, {}
        Yerr = {}
        for segment in self.segments:
            for trace in self.segments[segment]:
                # TODO: Fix this hack
                if trace == 'combined':
                    continue
                if trace not in X:
                    X[trace] = []
                    Y[trace] = []
                    Yerr[trace] = []
                for x in self.segments[segment][trace]:
                    X[trace].append(x)
                    _hist = self.segments[segment][trace][x]
                    Y[trace].append(_hist['mean'])
                    _sem = float(_hist['std'])/np.sqrt(self.opts.degfree - 1 or 1)
                    _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
                    Yerr[trace].append(_t_val)
        for trace in Y:
            ax.errorbar(X[trace],
                        Y[trace],
                        yerr=Yerr[trace],
                        marker='o',
                        lw=0,
                        elinewidth=0.25,
                        capsize=0.5,
                        label='Trace %s' % (trace+1))
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
        ax.legend(handles, labels)

    def PlotR(self, ax):
        if self.opts.logr:
            ax.set_title("Semilog Plot of |R|")
            ax.set_ylabel(r'log|R|')
        else:
            ax.set_title("Plot of |R|")
            ax.set_ylabel(r'|R|')
        ax.set_xlabel("Potenial (V)")
        Y, Yerr = [], []
        for x in self.XY:
            Y.append(self.XY[x]['R']["hist"]["mean"])
            _sem = float(self.XY[x]['R']["hist"]["std"])/np.sqrt(self.opts.degfree - 1 or 1)
            _t_val = _sem * stdtrit(self.opts.degfree - 1 or 1, 1 - self.opts.alpha)
            Yerr.append(_t_val)
        ax.errorbar(list(self.XY), Y,  yerr=Yerr, marker='o', lw=0.0, color='k')

    def PlotRFit(self, ax):
        key = list(self.XY)[-1]
        ax.set_title(r'Histogram and fit at $'+str(key)+'$ V')
        ax.set_xlabel(r'$log|R|$')
        ax.set_ylabel('Counts')
        ax.bar(self.XY[key]['R']['hist']['bin'], self.XY[key]['R']['hist']['freq'], width=0.05, color='r')
        ax.plot(self.XY[key]['R']['hist']['bin'], self.XY[key]['R']['hist']['fit'], lw=2.0, color='b', label='Fit')

    def PlotDJDV(self, ax):
        xax = list(self.DJDV.keys())
        xax.sort()
        ax.set_xlabel("Potential (V)")
        ax.set_title("Derivative of Initial Data")
        ax.set_ylabel(r'$\mathregular{\frac{dJ}{dV}}$')
        i = -1
        while True:
            i += 1
            try:
                if i not in self.ohmic:
                    ax.plot(xax, [self.DJDV[x][i] for x in xax], "-", lw=2)
                else:
                    ax.plot(xax, [self.DJDV[x][i] for x in xax], "-", lw=0.5, color='grey')
            except IndexError:
                break

    def PlotG(self, ax):
        ax.set_title("Conductance Plot")
        ax.set_xlabel("Potential (V)")
        if self.opts.heatmapd == 0:
            ax.set_title("Heatmap of Initial Data")
            ax.set_ylabel(r'$log|J A cm^{-2}|$')
        elif self.opts.heatmapd == 1:
            ax.set_title("Derivative of Initial Data")
            ax.set_ylabel(r'$\mathregular{\log|\frac{dJ}{dV}|}$')
        elif self.opts.heatmapd == 2:
            ax.set_title("Second Derivative of Initial Data")
            ax.set_ylabel(r'$\mathregular{\frac{d^2J}{dV^2}}$')
        x, y, z = [], [], []
        for v in self.GHists:
            for i in range(0, len(self.GHists[v]['hist']['bin'])):
                if i in self.ohmic and self.opts.skipohmic:
                    continue
                x.append(v)
                y.append(self.GHists[v]['hist']['bin'][i])
                z.append(self.GHists[v]['hist']['freq'][i])
        x, y, z = np.array(x), np.array(y), np.array(z)
        xmin, xmax = x.min(), x.max()
        ymin, ymax = self.opts.mlow, self.opts.mhi
        x = np.r_[x, xmin, xmax]
        y = np.r_[y, ymin, ymax]
        z = np.r_[z, z[0], z[-1]]
        xi = np.linspace(xmin, xmax, 200)
        yi = np.linspace(ymin, ymax, 200)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((x, y), z, (X, Y), method='nearest')
        ax.axis([xmin, xmax, ymin, ymax])
        ax.pcolormesh(X, Y, Z, cmap=get_cmap('rainbow'))

    def PlotNDC(self, ax):
        ax.set_title("NDC Plot")
        ax.set_xlabel("Potential (V)")
        ax.set_title("Heatmap of NDC")
        ax.set_ylabel('Normalized Differential Condutance')
        x, y, z = [], [], []
        for v in self.NDCHists:
            for i in range(0, len(self.NDCHists[v]['hist']['bin'])):
                if i in self.ohmic and self.opts.skipohmic:
                    continue
                x.append(v)
                y.append(self.NDCHists[v]['hist']['bin'][i])
                z.append(self.NDCHists[v]['hist']['freq'][i])
        if not x or not y or not z:
            logger.warn("Not enough data to plot NDC.")
            return
        x, y, z = np.array(x), np.array(y), np.array(z)
        xmin, xmax = x.min(), x.max()
        ymin, ymax = self.opts.ndc_mlow, self.opts.ndc_mhi
        x = np.r_[x, xmin, xmax]
        y = np.r_[y, ymin, ymax]
        z = np.r_[z, z[0], z[-1]]
        xi = np.linspace(xmin, xmax, 200)
        yi = np.linspace(ymin, ymax, 200)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((x, y), z, (X, Y), method='nearest')
        ax.axis([xmin, xmax, ymin, ymax])
        ax.pcolormesh(X, Y, Z, cmap=get_cmap('rainbow'))

    def PlotHist(self, ax):
        ax.set_title("Gaussian Fit and Raw Data")
        ax.set_xlabel('Potential (V)')
        ax.set_ylabel(r'Current Density $log_{10}|J(\mathrm{A cm}^{-2})|$')
        Y, Yerr = [], []
        for x in self.XY:
            Y.append(self.XY[x]["hist"]["mean"])
            Yerr.append(self.XY[x]["hist"]["std"])
        ax.errorbar(list(self.XY.keys()), Y, yerr=Yerr, lw=3.0, color='k')

    def PlotVtrans(self, ax):
        ax.set_title(r'Histogram and fit of $V_{trans}$')
        ax.set_xlabel(r'$V_{trans}$')
        ax.set_ylabel('Counts')
        for key in ('pos', 'neg'):
            ax.bar(self.FN[key]['bin'], self.FN[key]['freq'], width=0.01, color='g')
            ax.plot(self.FN[key]['bin'], self.FN[key]['fit'], lw=2.0, color='b', label='Fit')

    def PlotFit(self, ax):
        key = list(self.XY.keys())[-1]
        ax.set_title(r'Histogram and fit at $'+str(key)+'$ V')
        ax.set_xlabel(r'$log|J|$')
        ax.set_ylabel('Counts')
        ax.bar(self.XY[key]['hist']['bin'], self.XY[key]['hist']['freq'], width=0.1, color='r')
        ax.plot(self.XY[key]['hist']['bin'], self.XY[key]['hist']['fit'], lw=2.0, color='b', label='Fit')

    def DoPlots(self):
        ax1 = self.fig.add_axes([0.06, 0.55, 0.4, 0.4])
        ax2 = self.fig.add_axes([0.56, 0.55, 0.4, 0.4])
        ax3 = self.fig.add_axes([0.06, 0.05, 0.4, 0.4])
        ax4 = self.fig.add_axes([0.56, 0.05, 0.4, 0.4])
        if self.opts.plots == 'J':
            self.PlotFit(ax1)
            self.PlotData('LogY', ax2, ':', lw=0.25, color='c')
            self.PlotHist(ax2)
            ax4.set_ylim(ax2.get_ylim())
        elif self.opts.plots == 'R':
            self.PlotRFit(ax1)
            self.PlotR(ax2)
        elif self.opts.plots == 'FN':
            self.PlotVtrans(ax1)
            self.PlotData('LogY', ax2, ':', lw=0.25, color='c')
            self.PlotHist(ax2)
            ax4.set_ylim(ax2.get_ylim())
        elif self.opts.plots == 'SLM':
            self.PlotFit(ax1)
            self.PlotSLM(ax2)
        self.PlotSegmentedGauss(ax4)
        if self.opts.histplots == 'NDC':
            self.PlotNDC(ax3)
        elif self.opts.histplots == 'G':
            self.PlotG(ax3)
        
        # Show clustering results after main plots if --cluster is enabled
        if self.opts.cluster and hasattr(self.parser, 'cluster') and self.parser.cluster['clusterer'] is not None:
            try:
                # Create and display separate clustering plots
                clusterer = self.parser.cluster['clusterer']
                jv_curves = self.parser.cluster['jv_curves']
                
                # Show the comprehensive clustering results
                clustering_fig = clusterer.plot_results(jv_curves, figsize=(15, 10), log_scale=False)
                
                # Show 2D histograms if there are enough clusters
                if self.parser.cluster['n_clusters'] <= 5:
                    histogram_fig = clusterer.plot_2d_histograms(jv_curves, n_examples=3)
                
            except Exception as e:
                logger.warning(f"Failed to display clustering plots: {e}")
        
        if self.opts.write:
            # workaround for older versions of matplot lib / python3
            try:
                self.fig.savefig(os.path.join(self.opts.out_dir, self.opts.outfile+"_fig.png"), format="png")
            except AttributeError:
                logger.warning('Unable to save plots as image. Try upgrading python3 and matplotlib.')
