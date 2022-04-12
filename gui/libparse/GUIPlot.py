from tkinter import Toplevel, Button, BOTTOM, TOP, BOTH, X
from gaussfit.output import Plotter
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
    HASMATPLOTLIB = True
except ImportError:
    HASMATPLOTLIB = False


class PlotWindow(Toplevel):
    '''A new window to hold the plots'''

    def __init__(self, parent):
        super().__init__(parent)
        self.title("Data Preview")


def GUIPlot(self, parser):
    '''An interactive matplotlib pyplot that
   doesn't block the main GUI thread.'''

    if not HASMATPLOTLIB:
        self.logger.warn("Matplot lib is not installed, cannot generate plots!")
        return

    fig = Figure(figsize=(16, 10), dpi=100)
    plotter = Plotter(parser, fig)
    plotter.DoPlots()
    root = PlotWindow(self.master)
    canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
    canvas.draw()
    if int(matplotlib.__version__.split('.')[0]) > 2:
        toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=BOTTOM, fill=X)
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)


# def GUIPlotMultiprocessing(self, parser=None):
#     '''An interactive matplotlib pyplot that
#        doesn't block the main GUI thread.'''

#     if not HASMATPLOTLIB:
#         self.logger.warn("Matplot lib is not installed, cannot generate plots!")
#         return

#     if parser is None:
#         __plts = []
#         while self.plots:
#             __plts.append(self.plots.pop())
#         for __plt in __plts:
#             if __plt.get_fignums():
#                 __plt.gcf().canvas.draw_idle()
#                 __plt.gcf().canvas.start_event_loop(0.1)
#                 self.plots.append(__plt)
#         if self.plots:
#             self.master.after('100', self.GUIPlot)
#         return
#     self.logger.info("Generating plots in background...")
#     plotter = Plotter(parser, plt)
#     plotter.DoPlots()
#     plt.ion()
#     plt.pause(1)
#     self.plots.append(plt)
#     self.master.after('100', self.GUIPlot)
