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


# In your GUI class

def GUIPlot(self, parser):
    '''An interactive matplotlib pyplot that
   doesn't block the main GUI thread.'''

    if not HASMATPLOTLIB:
        self.logger.warn("Matplotlib is not installed, cannot generate plots!")
        return

    # --- Main Plot Window ---
    fig = Figure(figsize=(16, 10), dpi=100)
    plotter = Plotter(parser, fig)
    
    # Capture the extra figures returned by DoPlots
    extra_figs = plotter.DoPlots()

    root = PlotWindow(self.master)
    root.title("Main Data Plots") # Give it a more specific title
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    
    if int(matplotlib.__version__.split('.')[0]) > 2:
        toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=BOTTOM, fill=X)
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    # --- Create New Windows for Extra Plots ---
    for i, extra_fig in enumerate(extra_figs):
        # Create a new Toplevel window for each extra figure
        extra_root = Toplevel(self.master)
        extra_root.title(f"Clustering Plot {i+1}")
        
        extra_canvas = FigureCanvasTkAgg(extra_fig, master=extra_root)
        extra_canvas.draw()
        
        # Add a toolbar to the new window for interactivity
        extra_toolbar = NavigationToolbar2Tk(extra_canvas, extra_root, pack_toolbar=False)
        extra_toolbar.update()
        extra_toolbar.pack(side=BOTTOM, fill=X)
        
        extra_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)


