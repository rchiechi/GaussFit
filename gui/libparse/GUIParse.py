# from tempfile import NamedTemporaryFile
from tkinter import DISABLED, NORMAL
from gaussfit.output import Writer
from gaussfit.parse import Parse
from gaussfit.parse import readfiles
from .ParseThread import ParseThread
from gaussfit.output.libwriter import doOutput
from multiprocessing import Queue
from gaussfit.logger import DelayedMultiprocessHandler


def GUIParse(self):
    '''Start the parser in a background thread and wait for it to complete.'''
    def preParse():
        '''We need to check a couple of things right before we start parsing'''
        
        self.ButtonParse['state'] = DISABLED
        self.degfreedom = self.opts.degfree
        if self.opts.degfree == 1 and len(self.opts.in_files) > 1:
            self.opts.degfree = len(self.opts.in_files)-1

    def postParse():
        '''We need to check a couple of things right after we finish parsing'''
        self.ButtonParse['state'] = NORMAL
        self.opts.degfree = self.degfreedom
    try:
        if self.gothreads:
            for _t in self.gothreads:
                if _t.is_alive():
                    self.ButtonParse.after('100', self.GUIParse)
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
                df =  readfiles(self.opts.in_files)
                preParse()
                self.gothreads.append(ParseThread(Parse(df, handler=self.handler)))
                self.gothreads[-1].start()
                self.ButtonParse.after('500', self.GUIParse)
            else:
                self.logger.warning("No input files!")
        self.handler.flush()
    except Exception as e:
        print(f"Unhandled exception in GUIParse: {e}")
        import traceback
        traceback.print_exc()
        raise