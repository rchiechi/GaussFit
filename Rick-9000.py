#!/usr/bin/env python3
'''
Copyright (C) 2022 Ryan Chiechi <ryan.chiechi@ncsu.edu>
Description:

        This is a GUI front-end for parsing Seebeck data from the Rick-9000.

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

import os
import sys
import platform
import logging
from types import SimpleNamespace
import tkinter.ttk as tk
from tkinter import Tk
# from tkinter import Toplevel
from tkinter import filedialog
import tkinter.scrolledtext as ScrolledText
from tkinter import Text, IntVar, StringVar, Listbox, Label, Toplevel, TclError
from tkinter import N, S, E, W, X, Y  # pylint: disable=unused-import
from tkinter import TOP, BOTTOM, LEFT, RIGHT  # pylint: disable=unused-import
from tkinter import END, BOTH, VERTICAL, HORIZONTAL  # pylint: disable=unused-import
from tkinter import EXTENDED, RAISED, SOLID, DISABLED, NORMAL  # pylint: disable=unused-import
from tkinter import PhotoImage
from tkinter.font import Font

BLACK = "#000000"
YELLOW = "#f4e012"
WHITE = "#ffffff"
RED = "#ff0000"
TEAL = "#78CBFD"
GREEN = "#09f218"
BLUE = "#090df2"
GREY = '#e8e8e8'

absdir = os.path.dirname(os.path.abspath(__file__))


class TextHandler(logging.Handler):
    # This class allows you to log to a Tkinter Text or ScrolledText widget
    # Adapted from Moshe Kaplan: https://gist.github.com/moshekaplan/c425f861de7bbf28ef06

    def __init__(self, text):
        # run the regular Handler __init__
        logging.Handler.__init__(self)
        # Store a reference to the Text it will log to
        self.text = text

    def emit(self, record):
        msg = self.format(record)

        def append():
            self.text.configure(state='normal')
            self.text.insert(END, msg + '\n')
            self.text.configure(state='disabled')
            # Autoscroll to the bottom
            self.text.yview(END)
        # This is necessary because we can't modify the Text from other threads
        self.text.after(0, append)


class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0
        self.text = str()

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 27
        y = y + cy + self.widget.winfo_rooty() + 27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        try:
            # For Mac OS
            tw.tk.call("::tk::unsupported::MacWindowStyle",
                       "style", tw._w,
                       "help", "noActivates")
        except TclError:
            pass
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


def createToolTip(widget, text):
    toolTip = ToolTip(widget)

    def enter(event):
        toolTip.showtip(text)

    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)


class ChooseFiles(tk.Frame):
    '''The main frame for adding/removing files, accessing setttings
       and parsing.'''

    from seebeck.parse import GUIParse

    opts = SimpleNamespace(in_files=[],
                           out_dir=os.path.join(os.getcwd(), 'parsed'),
                           out_file='',
                           plot=True,
                           write=True,
                           truetemp=False,
                           col_to_parse=1,
                           dTn=[],
                           cuttoff_to_toss=100)  # TODO: Make this a user option
    raw_data = {}
    gothreads = []
    plots = []
    outdir = ''
    boolmap = {1: True, 0: False}
    colmap = {1:'Raw Voltage (uV)',
              2:'Corrected Voltage (uV)',
              3:'Top T',
              4:'Bottom T',
              5:'Delta-T (°C)',
              6:'Seebeck (uV/K)'}

    FileListBoxFrameLabelVar = None
    FileListBoxFrameLabel = None
    filelist = None
    FileListBox = None
    checks = []
    OutputFileName = None
    DeltaTn = None
    OptionsColString = None
    OptionCol = None

    def __init__(self, master=None):
        if master is None:
            master = Tk()
        super().__init__(master)
        bgimg = PhotoImage(file=os.path.join(absdir, 'gui', 'RCCLabFluidic.png'))
        limg = Label(self.master, i=bgimg, background=GREY)
        limg.pack(side=TOP)
        try:
            self.last_input_path = os.getcwd()
        except KeyError:
            self.last_input_path = os.path.expanduser('~')
        master.tk_setPalette(background=GREY, activeBackground=GREY)
        master.title("RCCLab Rick-9000 Parser")
        master.geometry('800x850+250+250')
        self.pack(fill=BOTH)
        if len(sys.argv) > 1:
            self.opts.in_files = sys.argv[1:]
        self.__createWidgets()
        self.ToFront()

    def __createWidgets(self):

        self.ButtonFrame = tk.Frame(self)
        self.LoggingFrame = tk.Frame(self)
        self.RightOptionsFrame = tk.Frame(self)
        self.FileListBoxFrame = tk.Frame(self)

        # ScrolledText widget to display logging output
        stlogger = ScrolledText.ScrolledText(self.LoggingFrame, state='disabled')
        stlogger.configure(font='TkFixedFont')
        stlogger.pack(side=LEFT, fill=BOTH, expand=True)
        # Create textLogger
        text_handler = TextHandler(stlogger)

        # Logging configuration
        logging.basicConfig(filename='Rick-9000.log',
                            level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')

        # Add the handler to logger
        logger = logging.getLogger(__package__)
        logger.addHandler(text_handler)

        self.__createButtons()
        self.__createFileListBox()
        self.__createOptions()

        self.ButtonFrame.pack(side=BOTTOM, fill=None)
        self.FileListBoxFrame.pack(side=TOP, fill=BOTH, expand=True)
        self.RightOptionsFrame.pack(side=RIGHT, fill=Y)
        self.LoggingFrame.pack(side=BOTTOM, fill=BOTH)

        # self.logger.setLevel(getattr(logging, self.opts.loglevel.upper()))
        self.updateFileListBox()
        self.checkOptions()
        logging.info('Select some files to get started.')

    def __createButtons(self):

        buttons = [{'name': 'Quit', 'text': 'QUIT', 'command': 'Quit', 'side': BOTTOM},
                   {'name': 'SpawnInputDialog', 'text': 'Add Input Files', 'side': LEFT},
                   {'name': 'RemoveFile', 'text': 'Remove Files', 'side': LEFT},
                   {'name': 'SpawnOutputDialog', 'text': 'Choose Output Directory', 'side': LEFT},
                   {'name': 'Parse', 'text': 'Parse!', 'side': LEFT}
                   ]

        for b in buttons:
            button = tk.Button(self.ButtonFrame)
            button.config(text=b['text'], command=getattr(self, b['name']+'Click'))
            button.pack(side=b['side'])
            setattr(self, 'Button'+b['name'], button)

    def __createFileListBox(self):
        self.FileListBoxFrameLabelVar = StringVar()
        self.FileListBoxFrameLabel = tk.Label(self.FileListBoxFrame,
                                              textvariable=self.FileListBoxFrameLabelVar,
                                              font=Font(size=10, weight="bold"))
        self.FileListBoxFrameLabel.pack(side=TOP, fill=X)

        yScroll = tk.Scrollbar(self.FileListBoxFrame, orient=VERTICAL)
        yScroll.pack(side=RIGHT, fill=Y)
        xScroll = tk.Scrollbar(self.FileListBoxFrame, orient=HORIZONTAL)
        xScroll.pack(side=BOTTOM, fill=X)
        self.filelist = StringVar()
        self.FileListBox = Listbox(self.FileListBoxFrame,
                                   listvariable=self.filelist, selectmode=EXTENDED,
                                   height=20, width=0, relief=RAISED, bd=1,
                                   bg=WHITE,
                                   # font=Font(size=10),
                                   xscrollcommand=xScroll.set,
                                   yscrollcommand=yScroll.set)
        xScroll['command'] = self.FileListBox.xview
        yScroll['command'] = self.FileListBox.yview
        self.FileListBox.pack(side=LEFT, fill=BOTH, expand=True)
        self.UpdateFileListBoxFrameLabel()

    def __createOptions(self):
        '''Create the widgets for options and use setattr to assign them to self.'''
        self.checks = [{'name': 'plot', 'text': 'Plot', 'row': 1,
                        'tooltip': "Show summary plots after parsing."},
                       {'name': 'write', 'text': 'Write', 'row': 2,
                        'tooltip': "Write results to text files after parsing."},
                       {'name': 'truetemp', 'text': 'Use Real T', 'row': 3,
                        'tooltip': "Use the measured ΔT instead the ΔT field values."}
                       ]

        for _c in self.checks:
            setattr(self, _c['name'], IntVar())
            check = tk.Checkbutton(self.RightOptionsFrame, text=_c['text'],
                                   variable=getattr(self, _c['name']),
                                   command=self.checkOptions)
            check.grid(column=0, row=_c['row'], sticky=W)
            createToolTip(check, _c['tooltip'])
            setattr(self, 'Check_'+_c['name'], check)
            if getattr(self.opts, _c['name']):
                getattr(self, _c['name']).set(1)

        rowidx = len(self.checks)+1

        tk.Label(self.RightOptionsFrame, text="Output file base name:").grid(
            column=0, row=rowidx)
        self.OutputFileName = tk.Entry(self.RightOptionsFrame, width=20,
                                       font=Font(size=10, slant='italic'))
        for n in ('<Return>', '<Leave>', '<Enter>'):
            self.OutputFileName.bind(n, self.checkOutputFileName)
        self.OutputFileName.grid(column=0, row=rowidx+1)

        if self.opts.out_file:
            self.OutputFileName.insert(0, self.opts.out_file)

        tk.Label(self.RightOptionsFrame, text="ΔT values:").grid(
            column=0, row=rowidx+2, sticky=W)
        self.DeltaTn = tk.Entry(self.RightOptionsFrame, width=20,
                                font=Font(size=10))
        self.DeltaTn.insert(0, '4,8,12')
        for n in ('<Return>', '<Leave>', '<Enter>'):
            self.DeltaTn.bind(None, self.checkOptions)
        self.DeltaTn.grid(column=0, row=rowidx+3, sticky=W)

        tk.Label(self.RightOptionsFrame, text="Column to plot:").grid(
            column=0, row=rowidx+4, sticky=W)
        self.OptionsColString = StringVar()
        self.OptionsColString.set(self.opts.col_to_parse)
        self.OptionCol = tk.OptionMenu(self.RightOptionsFrame,
                                       self.OptionsColString,
                                       self.colmap[self.opts.col_to_parse],
                                       command=self.checkOptions,
                                       *list(self.colmap.values()))
        # __menu = self.nametowidget(self.OptionCol)
        # __menu.config(font=Font(size=10))  # Set the dropdown menu's font
        self.OptionCol.grid(column=0, row=rowidx+5, sticky=W)

        # tk.Label(self.RightOptionsFrame, text="ΔT values:").grid(
        #     column=0, row=rowidx+2, sticky=W)
        # self.OptionsdTnString = StringVar()
        # self.OptionsdTnString.set('.'.join(self.opts.dTn))
        # self.OptionPlots = tk.OptionMenu(self.RightOptionsFrame,
        #                                  self.OptionsdTnString, self.opts.dTn, 'J', 'R',
        #                                  command=self.checkOptions)
        # self.OptionPlots.grid(column=0, row=rowidx+3, sticky=W)

    def RemoveFileClick(self):
        self.checkOptions()
        selected = [self.FileListBox.get(x) for x in self.FileListBox.curselection()]
        todel = []
        filelist = []
        for i in range(0, len(self.opts.in_files)):
            for _s in selected:
                if self.opts.in_files[i].replace(" ", "_") == _s:
                    todel.append(i)

        for i in range(0, len(self.opts.in_files)):
            if i not in todel:
                filelist.append(self.opts.in_files[i])
        self.opts.in_files = filelist
        self.updateFileListBox()
        self.FileListBox.selection_clear(0, END)

    def SpawnInputDialogClick(self):
        # self.checkOptions()
        self.opts.in_files += filedialog.askopenfilename(
            title="Files to parse",
            multiple=True,
            initialdir=self.last_input_path,
            filetypes=[('LabView Files', '*.lvm'), ('Data files', '*_data.txt'),
                       ('Text files', '*.txt'), ('All files', '*')])
        if self.opts.in_files:
            self.last_input_path = os.path.split(self.opts.in_files[0])[0]
            if not self.outdir:
                self.opts.out_dir = self.last_input_path
            self.updateFileListBox()
            if not self.opts.out_file:
                self.OutputFileName.delete(0, END)
                self.OutputFileName.insert(0, os.path.basename(
                    self.opts.in_files[-1]).lower().replace('.lvm', ''))
                self.checkOutputFileName()
        self.updateFileListBox()

    def ParseClick(self):
        self.checkOptions()
        self.GUIParse()

    def checkOptions(self, m=None):
        for c in self.checks:
            setattr(self.opts, c['name'], self.boolmap[getattr(self, c['name']).get()])
        if self.opts.truetemp:
            self.DeltaTn['state'] = DISABLED
        else:
            self.DeltaTn['state'] = NORMAL
        self.opts.dTn = self.DeltaTn.get().split(',')

        for __key in self.colmap:
            if self.colmap[__key] == self.OptionsColString.get():
                self.opts.col_to_parse = __key

    def updateFileListBox(self):
        self.filelist.set(" ".join([x.replace(" ", "_") for x in self.opts.in_files]))

    def SpawnOutputDialogClick(self):
        outdir = filedialog.askdirectory(title="Select Output File(s)",
                                         initialdir=self.opts.out_dir)
        if not outdir:
            return
        if os.path.exists(outdir):
            self.opts.out_dir = outdir
            self.outdir = self.opts.out_dir  # So we know the user set the output dir
            self.UpdateFileListBoxFrameLabel()

    def UpdateFileListBoxFrameLabel(self):
        self.FileListBoxFrameLabelVar.set(
            f"Output to: {self.opts.out_dir}/{self.opts.out_file}_*.txt")

    def checkOutputFileName(self, event=None):
        self.opts.out_file = self.OutputFileName.get()
        self.UpdateFileListBoxFrameLabel()

    def QuitClick(self):
        self.Quit()

    def Quit(self):
        self.master.destroy()

    def ToFront(self):
        '''Try to bring the main window to the front on different platforms'''
        if platform.system() == "Darwin":
            os.system(
                '''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        else:
            self.master.attributes('-topmost', 1)
            self.master.attributes('-topmost', 0)
        self.master.lift()


if __name__ == '__main__':
    print("Starting GUI")
    root = Tk()
    gui = ChooseFiles(root)
    gui.mainloop()
