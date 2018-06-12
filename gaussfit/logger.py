
#!/usr/bin/env python3
'''
Copyright (C) 2018 Ryan Chiechi <r.c.chiechi@rug.nl>

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
from collections import Counter,OrderedDict

#TODO Isn't this done with QueueHandler?
class DelayedHandler(logging.Handler):
    '''A log handler that buffers messages and 
    folds repeat messages into one line.'''

    buff = []

    def __init__(self,delay=False):
        logging.Handler.__init__(self)
        self.createLock()
        self._delay = delay

    def emit(self, message): # Overwrites the default handler's emit method
        self.buff.append(message)
        if not self._delay:
            self.flush()

    def _emit(self,message,level):
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
        msgs = Counter(map(self.format,self.buff))
        emitted = []
        for message in self.buff:
            #FIFO
            fmsg = self.format(message)
            if fmsg not in emitted:
                emitted.append(fmsg)
                i = msgs[fmsg]
                if i > 1:
                    self._emit('%s (repeated %s times)' % (self.format(message),i), message.levelname)
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

    from tkinter import NORMAL,DISABLED,END

    def __init__(self, console, delay=False):
        DelayedHandler.__init__(self,delay)
        self.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.console = console #Any text widget, you can use the class above or not

    def _emit(self,message,level):
        self.console["state"] = self.NORMAL
        self.console.insert(self.END, message+"\n") #Inserting the logger message in the widget
        self.console["state"] = self.DISABLED
        self.console.see(self.END)

