from gaussfit.output import Plotter
try:
    import matplotlib.pyplot as plt
    HASMATPLOTLIB = True
except ImportError:
    HASMATPLOTLIB = False


def GUIPlot(self, parser=None):
    '''An interactive matplotlib pyplot that
       doesn't block the main GUI thread.'''

    if not HASMATPLOTLIB:
        self.logger.warn("Matplot lib is not installed, cannot generate plots!")
        return

    if parser is None:
        __plts = []
        while self.plots:
            __plts.append(self.plots.pop())
        for __plt in __plts:
            if __plt.get_fignums():
                __plt.gcf().canvas.draw_idle()
                __plt.gcf().canvas.start_event_loop(0.1)
                self.plots.append(__plt)
        if self.plots:
            self.master.after('100', self.GUIPlot)
        return
    self.logger.info("Generating plots...")
    plotter = Plotter(parser, plt)
    plotter.DoPlots()
    plt.ion()
    plt.pause(1)
    self.plots.append(plt)
    self.master.after('100', self.GUIPlot)
