import threading


class ParseThread(threading.Thread):
    '''A Thread object to run the parse in so it
       doesn't block the main GUI thread.'''

    def __init__(self, parser):
        super().__init__()
        self.parser = parser

    def run(self):
        self.parser.readfiles(self.parser.opts.in_files)
