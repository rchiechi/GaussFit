import asyncio
import threading
from gaussfit.parse import Parse
from gaussfit.parse import readfiles
from gaussfit.logger import DelayedMultiprocessHandler

class ParseThread(threading.Thread):
    """
    A wrapper thread that creates a Parse object and runs its async .parse() method
    in a separate event loop, allowing the GUI to remain responsive.
    """
    def __init__(self, opts, logque):
        super().__init__()
        # Dependencies needed to create the Parse object
        self.opts = opts
        self.logque = logque
        # This will be populated by the thread with the finished parser instance
        self.parser = None
        
        self.loop = None
        self.parse_future = None

    async def _create_and_run_parser(self):
        """
        This is the actual async work for the thread. It creates the Parse
        object and then calls its .parse() method.
        """
        # This logic is what was in your original ParseThread.parse()
        df = await readfiles(self.opts)
        self.parser = Parse(df,
                    opts=self.opts,
                    handler=DelayedMultiprocessHandler(self.logque),
                    logque=self.logque)
        await self.parser.parse() # Call the method ON the instance

    def run(self):
        """
        Manages the event loop and schedules the parser task.
        The loop will stop and the thread will terminate when the task is done.
        """
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        
        try:
            # Create a task from our wrapper method
            task = self.loop.create_task(self._create_and_run_parser())
            
            # Add the callback to stop the loop when the task finishes
            task.add_done_callback(lambda future: self.loop.stop())
            
            self.parse_future = task
            self.loop.run_forever()
        finally:
            self.loop.close()
            print("ParseThread event loop closed and thread is terminating.")

    def stop(self):
        """External command to gracefully stop the thread mid-parse."""
        if self.loop and self.loop.is_running():
            if self.parse_future and not self.parse_future.done():
                self.parse_future.cancel()
            self.loop.call_soon_threadsafe(self.loop.stop)
