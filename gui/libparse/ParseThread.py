import threading


class ParseThread(threading.Thread):
    '''A Thread object to run the parse in so it
       doesn't block the main GUI thread.'''

    def __init__(self, opts, parser):
        threading.Thread.__init__(self)
        self.parser = parser
        self.opts = opts

    def run(self):
        self.parser.readfiles(self.opts.in_files)
