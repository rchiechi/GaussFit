import asyncio
import threading
from gaussfit.parse import Parse
from gaussfit.parse import readfiles
from gaussfit.logger import DelayedMultiprocessHandler

class ParseThread(threading.Thread):
    '''A Thread object to run the parse in so it
       doesn't block the main GUI thread.'''

    def __init__(self, opts, logque):
        super().__init__()
        self.opts = opts
        self.logque = logque
        self.parser = None

    async def parse(self):
        df = await readfiles(self.opts)
        self.parser = Parse(df,
                            opts=self.opts,
                            handler=DelayedMultiprocessHandler(self.logque),
                            logque=self.logque)
        await self.parser.parse()
        
    def run(self):
        asyncio.run(self.parse())
