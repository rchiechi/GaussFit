import sys
import logging
# from logging.handlers import QueueHandler
import queue
# import hashlib
from multiprocessing import Queue
from collections import Counter
from gaussfit.colors import GREEN, TEAL, WHITE, YELLOW

FMT = GREEN+'%(name)s'+TEAL+' %(levelname)s '+YELLOW+'%(message)s'+WHITE


def emittoconsole(message, level):
    if level in ('INFO' 'DEBUG'):
        sys.stdout.write(message)
        sys.stdout.write("\n")
        sys.stdout.flush()
    else:
        sys.stderr.write(message)
        sys.stderr.write("\n")
        sys.stderr.flush()


class DelayedHandler(logging.Handler):
    '''A log handler that buffers messages and
    folds repeat messages into one line.'''

    def __init__(self, buff=None):
        super().__init__()
        self.buff = buff or queue.Queue(-1)
        self.createLock()
        self._delay = False

    def emit(self, message):  # Overwrites the default handler's emit method
        self.buff.put_nowait(message)
        if not self._delay:
            self.flush()

    def _emit(self, message, level):
        emittoconsole(message, level)

    def flush(self):
        _buff = []
        while not self.buff.empty():
            try:
                __record = self.buff.get(1)
                if __record is not None:
                    _buff.append(__record)
            except EOFError:
                return
            except queue.Empty:
                return
        try:
            msgs = Counter(map(self.format, _buff))
        except TypeError as msg:
            emittoconsole('Error formatting message buffer: %s' % str(msg), 'ERROR')
            for __buff in _buff:
                try:
                    __buff.getMessage()
                except TypeError:
                    self._emit(str(__buff), 'ERROR')

        emitted = []
        for message in _buff:
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
        # self.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.console = console  # Any text widget, you can use the class above or not

    def _emit(self, message, level):
        self.console["state"] = self.NORMAL
        self.console.insert(self.END, message+"\n")  # Inserting the logger message in the widget
        self.console["state"] = self.DISABLED
        self.console.see(self.END)
        emittoconsole(message, level)


class DelayedMultiprocessHandler(logging.Handler):
    '''A dummy-delayed log handler that queues unique log messages for running
       in a background thread.'''

    def __init__(self, que):
        super().__init__()
        self.que = que
        self.buff = []

    def emit(self, record):
        self.buff.append(record.lineno)
        if Counter(self.buff)[record.lineno] == 1:
            self.que.put_nowait(record)
        emittoconsole(logging.Formatter(FMT).format(record), record.levelname)

    def setDelay(self):
        return

    def unsetDelay(self):
        return

    def flush(self):
        return
