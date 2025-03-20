import asyncio
import threading
from gaussfit.parse import Parse
from gaussfit.parse import readfiles
from gaussfit.logger import DelayedMultiprocessHandler

class ParseThread(threading.Thread):
    '''A Thread object to run the parse in so it
       doesn't block the main GUI thread.'''

    def __init__(self, in_files, handler):
        super().__init__()
        self.in_files = in_files
        self.handler = handler
        self.parser = None

    async def parse(self):
        df = await readfiles(self.in_files)
        self.parser = Parse(df, handler=self.handler)
        await self.parser.parse()
        
    def run(self):
        asyncio.run(self.parse())
