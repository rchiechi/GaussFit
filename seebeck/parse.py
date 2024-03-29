import os
from tkinter import DISABLED, NORMAL
from queue import Queue
import logging
import threading
from seebeck.util import get_raw_data
from seebeck.stats import linear_fit, GHistograms
from seebeck.writer import doOutput


def GUIParse(self):
    '''Start the parser in a background thread and wait for it to complete.'''

    def preParse():
        '''We need to check a couple of things right before we start parsing'''
        self.logque = Queue(-1)
        self.ButtonParse['state'] = DISABLED
        self.ButtonQuit['state'] = DISABLED

        if not os.path.exists(self.opts.out_dir):
            os.mkdir(self.opts.out_dir)
        # self.degfreedom = self.opts.degfree
        # if self.opts.degfree == 0 and len(self.opts.in_files):
        #     self.opts.degfree = len(self.opts.in_files)-1

    def postParse():
        '''We need to check a couple of things right after we finish parsing'''
        try:
            _linear = linear_fit(self.opts, self.raw_data)
            _histograms = GHistograms(self.opts, self.raw_data)
            if self.opts.write:
                doOutput(self.opts, linear_fit=_linear, histograms=_histograms)
        except ValueError:
            logging.error("Could not generate fits.")
        # print(self.raw_data)
        # self.opts.degfree = self.degfreedom
        self.ButtonParse['state'] = NORMAL
        self.ButtonQuit['state'] = NORMAL

    if self.gothreads:
        while not self.logque.empty():
            logging.info(self.logque.get_nowait())
        for _t in self.gothreads:
            if _t.is_alive():
                self.ButtonParse.after('100', self.GUIParse)
                return
        for _t in self.gothreads:
            _t.join()
        logging.info("Parse complete!")
        self.gothreads.pop()
        postParse()
        # if self.opts.write and not gothread.parser.error:
        #     writer = Writer(gothread.parser)
        #     doOutput(writer)
        # if self.opts.plot and not gothread.parser.error:
        #     self.GUIPlot(parser=gothread.parser)

    else:
        if self.opts.in_files:
            preParse()
            self.gothreads.append(threading.Thread(target=get_raw_data,
                                                   args=[self.opts, self.logque, self.raw_data]))
            self.gothreads[-1].start()
            self.ButtonParse.after('500', self.GUIParse)
        else:
            logging.warning("No input files!")
