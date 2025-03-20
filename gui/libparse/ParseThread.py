import asyncio
import threading
from gaussfit.parse import Parse
from gaussfit.parse import readfiles

class ParseThread(threading.Thread):
    '''A Thread object to run the parse in so it
       doesn't block the main GUI thread.'''

    def __init__(self, parser):
        super().__init__()
        self.parser = parser

    def run(self):
        asyncio.run(self.parser.parse())
