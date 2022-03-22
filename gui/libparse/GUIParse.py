from tkinter import DISABLED, NORMAL
from gaussfit.output import Writer, Plotter
from gaussfit import Parse
from gui.libparse import ParseThread
from gaussfit.output.libwriter import doOutput


def GUIParse(self):
    '''Start the parser in a background thread and wait for it to complete.'''

    def preParse():
        '''We need to check a couple of things right before we start parsing'''
        self.degfreedom = self.opts.degfree
        if self.opts.degfree == 0 and len(self.opts.in_files):
            self.opts.degfree = len(self.opts.in_files)-1

    def postParse():
        '''We need to check a couple of things right after we finish parsing'''
        self.opts.degfree = self.degfreedom

    preParse()
    if self.gothreads:
        self.ButtonParse['state'] = DISABLED
        for _t in self.gothreads:
            if _t.is_alive():
                self.ButtonParse.after('500', self.GUIParse)
                return
        self.logger.info("Parse complete!")
        gothread = self.gothreads.pop()
        self.ButtonParse['state'] = NORMAL
        if self.opts.write and not gothread.parser.error:
            writer = Writer(gothread.parser)
            doOutput(writer)
        if self.opts.plot and not gothread.parser.error:
            try:
                import matplotlib.pyplot as plt
                plotter = Plotter(gothread.parser, plt)
                self.logger.info("Generating plots...")
                plotter.DoPlots()
                plt.show(block=False)
            except ImportError as msg:
                self.logger.error("Cannot import matplotlib! %s", str(msg), exc_info=False)
    else:
        if self.opts.in_files:
            parser = Parse(self.opts, handler=self.handler)
            self.gothreads.append(ParseThread(self.opts, parser))
            self.gothreads[-1].start()
            self.ButtonParse.after('500', self.GUIParse)
        else:
            self.logger.warning("No input files!")
    postParse()
    self.handler.flush()
