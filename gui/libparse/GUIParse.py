# from tempfile import NamedTemporaryFile
from tkinter import DISABLED, NORMAL
from gaussfit.output import Writer
from gaussfit import Parse
from gui.libparse import ParseThread
from gaussfit.output.libwriter import doOutput
from multiprocessing import Queue
from gaussfit.logger import DelayedMultiprocessHandler


def GUIParse(self):
    '''Start the parser in a background thread and wait for it to complete.'''

    def preParse():
        '''We need to check a couple of things right before we start parsing'''
        self.logque = Queue(-1)
        self.ButtonParse['state'] = DISABLED
        self.degfreedom = self.opts.degfree
        if self.opts.degfree == 0 and len(self.opts.in_files):
            self.opts.degfree = len(self.opts.in_files)-1

    def postParse():
        '''We need to check a couple of things right after we finish parsing'''
        self.ButtonParse['state'] = NORMAL
        self.opts.degfree = self.degfreedom

    if self.gothreads:
        for _t in self.gothreads:
            if _t.is_alive():
                self.ButtonParse.after('500', self.GUIParse)
                while not self.logque.empty():
                    self.handler.emit(self.logque.get_nowait())
                self.handler.flush()
                return
        self.logger.info("Parse complete!")
        gothread = self.gothreads.pop()
        postParse()
        if self.opts.write and not gothread.parser.error:
            writer = Writer(gothread.parser)
            doOutput(writer)
        if self.opts.plot and not gothread.parser.error:
            self.GUIPlot(parser=gothread.parser)

    else:
        if self.opts.in_files:
            preParse()
            parser = Parse(self.opts, handler=DelayedMultiprocessHandler(self.logque))
            self.gothreads.append(ParseThread(self.opts, parser))
            self.gothreads[-1].start()
            self.ButtonParse.after('500', self.GUIParse)
        else:
            self.logger.warning("No input files!")
    self.handler.flush()
