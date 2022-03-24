import sys
import logging
from multiprocessing import Queue
from collections import Counter


# TODO Isn't this done with QueueHandler?
class DelayedHandler(logging.Handler):
    '''A log handler that buffers messages and
    folds repeat messages into one line.'''

    buff = []

    def __init__(self):
        logging.Handler.__init__(self)
        self.createLock()
        self._delay = False

    def emit(self, message):  # Overwrites the default handler's emit method
        self.buff.append(message)
        if not self._delay:
            self.flush()

    def _emit(self, message, level):
        self.acquire()
        if level in ('INFO' 'DEBUG'):
            sys.stdout.write(message)
            sys.stdout.write("\n")
            sys.stdout.flush()
        else:
            sys.stderr.write(message)
            sys.stderr.write("\n")
            sys.stderr.flush()
        self.release()

    def flush(self):
        try:
            msgs = Counter(map(self.format, self.buff))
        except TypeError as msg:
            self._emit('Error formatting message buffer: %s' % str(msg), 'ERROR')
            for _buff in self.buff:
                try:
                    _buff.getMessage()
                except TypeError:
                    self._emit(str(_buff), 'ERROR')
            self.buff = []
        emitted = []
        for message in self.buff:
            if not str(message).strip():
                continue
            # FIFO
            fmsg = self.format(message)
            if fmsg not in emitted:
                emitted.append(fmsg)
                i = msgs[fmsg]
                if i > 1:
                    self._emit('%s (repeated %s times)' % (self.format(message), i), message.levelname)
                else:
                    self._emit(self.format(message), message.levelname)
        self.buff = []

    def setDelay(self):
        self._delay = True

    def unsetDelay(self):
        self._delay = False
        self.flush()


class GUIHandler(DelayedHandler):
    '''A log handler that buffers messages and folds repeats
    into a single line. It expects a tkinter widget as input.'''

    from tkinter import NORMAL, DISABLED, END

    def __init__(self, console):
        DelayedHandler.__init__(self)
        self.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.console = console  # Any text widget, you can use the class above or not

    def _emit(self, message, level):
        self.console["state"] = self.NORMAL
        self.console.insert(self.END, message+"\n")  # Inserting the logger message in the widget
        self.console["state"] = self.DISABLED
        self.console.see(self.END)


class DelayedMultiprocessHandler(DelayedHandler):

    def __init__(self, que):
        DelayedHandler.__init__(self)
        self.que = que

    def _emit(self, message, level):
        self.que.put(message)
