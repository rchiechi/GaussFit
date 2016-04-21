
#!/usr/bin/env python3
'''
Copyright (C) 2015 Ryan Chiechi <r.c.chiechi@rug.nl>
Description:
        This program parses raw current-voltage data obtained from
        molecular tunneling junctions. It is specifically designed
        with EGaIn in mind, but may be applicable to other systems.

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
    """ Used to redirect logging output to the widget passed in parameters """
    def __init__(self,delay=False):
        logging.Handler.__init__(self)
        self.__delay = delay
        self.buff = []
    def emit(self, message): # Overwrites the default handler's emit method
        self.buff.append(self.format(message))
        if not self.__delay:
            self.__flush()
    def __flush(self):
        msgs = Counter(self.buff)
        for msg in msgs:
            i = msgs[msg]
            if i > 1:
                print('%s (repeated %s times)' % (str(msg),i))
            else:
                print(str(msg))
        self.buff = []
    def setDelay(self):
        self.__delay = True
    def unsetDelay(self):
        self.__delay = False
        self.__flush()
    def flush(self):
        self.__flush()

class GUIHandler(logging.Handler):
    """ Used to redirect logging output to the widget passed in parameters """
    from tkinter import NORMAL,DISABLED,END
    def __init__(self, console, delay=False):
        logging.Handler.__init__(self)
        self.__delay = delay
        self.buff = []
        self.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.console = console #Any text widget, you can use the class above or not

    def emit(self, message): # Overwrites the default handler's emit method
        #formattedMessage = self.format(message)  #You can change the format here
        self.buff.append(self.format(message))
        if not self.__delay:
            self.__flush()
    def __emit(self,message):
        self.console["state"] = self.NORMAL
        self.console.insert(self.END, message+"\n") #Inserting the logger message in the widget
        self.console["state"] = self.DISABLED
        self.console.see(self.END)

    def __flush(self):
        msgs = Counter(self.buff)
        for msg in msgs:
            i = msgs[msg]
            if i > 1:
                self.__emit('%s (repeated %s times)' % (str(msg),i))
            else:
                self.__emit(str(msg))
        self.buff = []
       
    def setDelay(self):
        self.__delay = True
    def unsetDelay(self):
        self.__delay = False
        self.__flush()
    def flush(self):
        self.__flush()

