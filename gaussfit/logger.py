
#!/usr/bin/env python3
'''
Copyright (C) 2015 Ryan Chiechi <r.c.chiechi@rug.nl>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import sys,os,logging
from  gaussfit.colors import *
from collections import Counter


class DelayedHandler(logging.Handler):
    '''A log handler that buffers messages and 
    folds repeat messages into one line.'''

    buff = []

    def __init__(self,delay=False):
        logging.Handler.__init__(self)
        self._delay = delay
    def emit(self, message): # Overwrites the default handler's emit method
        self.buff.append(self.format(message))
        if not self._delay:
            self._flush()

    def _emit(self,msg):
        print(msg)

    def _flush(self):
        msgs = Counter(self.buff)
        emitted = []
        #for msg in msgs:
        for msg in self.buff:
            # FIFO
            if msg not in emitted:
                emitted.append(msg)
                i = msgs[msg]
                if i > 1:
                    self._emit('%s (repeated %s times)' % (str(msg),i))
                else:
                    self._emit(str(msg))
        self.buff = []
    def setDelay(self):
        self._delay = True
    def unsetDelay(self):
        self._delay = False
        self._flush()
    def flush(self):
        self._flush()

class GUIHandler(DelayedHandler):
    '''A log handler that buffers messages and folds repeats
    into a single line. It expects a tkinter widget as input.'''

    from tkinter import NORMAL,DISABLED,END

    def __init__(self, console, delay=False):
        DelayedHandler.__init__(self,delay)
        self.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.console = console #Any text widget, you can use the class above or not

    def _emit(self,message):
        self.console["state"] = self.NORMAL
        self.console.insert(self.END, message+"\n") #Inserting the logger message in the widget
        self.console["state"] = self.DISABLED
        self.console.see(self.END)

