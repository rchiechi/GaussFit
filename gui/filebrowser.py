'''
Copyright (C) 2022 Ryan Chiechi <ryan.chiechi@ncsu.edu>
Description:

        This is the GUI front-end for the parsing engine. It mostly works ok,
        but some of the options configurable on the command line may not be
        implemented.

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
import platform
import logging
import threading
import tkinter.ttk as tk
from tkinter import Tk
# from tkinter import Toplevel
from tkinter import filedialog
from tkinter import Text, IntVar, StringVar, Listbox, Label
from tkinter import N, S, E, W, X, Y  # pylint: disable=unused-import
from tkinter import TOP, BOTTOM, LEFT, RIGHT  # pylint: disable=unused-import
from tkinter import END, BOTH, VERTICAL, HORIZONTAL  # pylint: disable=unused-import
from tkinter import EXTENDED, RAISED, DISABLED, NORMAL  # pylint: disable=unused-import
from tkinter import PhotoImage
from tkinter.font import Font
from .prefs import PreferencesWindow
from .colors import BLACK, YELLOW, WHITE, RED, TEAL, GREEN, BLUE, GREY  # pylint: disable=unused-import
from .tooltip import createToolTip
from gaussfit.logger import GUIHandler
from gaussfit.args import Opts, VERSION

try:
    import psutil
    HASPSUTIL = True
except ImportError:
    HASPSUTIL = False

absdir = os.path.dirname(os.path.realpath(__file__))
absroot = Tk()


class ChooseFiles(tk.Frame):
    '''The main frame for adding/removing files, accessing settings
       and parsing.'''

    from .libparse import GUIParse
    from .libparse import GUIPlot

    gothreads = []
    plots = []
    outdir = ''
    boolmap = {1: True, 0: False}
    lock = threading.Lock

    def __init__(self):
        self.master = absroot
        super().__init__(self.master)
        # tk.Frame.__init__(self, self.master)
        bgimg = PhotoImage(file=os.path.join(absdir, 'RCCLabFluidic.png'))
        limg = Label(self.master, i=bgimg, background=GREY)
        limg.pack(side=TOP)
        try:
            self.last_input_path = os.getcwd()
        except KeyError:
            self.last_input_path = os.path.expanduser('~')
        self.opts = Opts
        self.degfreedom = {'init': self.opts.degfree, 'user': self.opts.degfree}
        self.master.tk_setPalette(background=GREY, activeBackground=GREY)
        self.master.title(f"RCCLab EGaIn Data Parser v{VERSION}")
        self.master.geometry('800x1000+250+250')
        self.pack(fill=BOTH)
        self.__createWidgets()
        self.ToFront()
        self.master.mainloop()

    def __createWidgets(self):

        self.ButtonFrame = tk.Frame(self)
        # self.LeftOptionsFrame = tk.Frame(self)
        self.LoggingFrame = tk.Frame(self)
        self.RightOptionsFrame = tk.Frame(self)
        self.FileListBoxFrame = tk.Frame(self)

        yScroll = tk.Scrollbar(self.LoggingFrame, orient=VERTICAL)
        self.Logging = Text(self.LoggingFrame, height=40, width=0,
                            bg=BLACK, fg=WHITE, yscrollcommand=yScroll.set)
        yScroll['command'] = self.Logging.yview
        self.Logging.pack(side=LEFT, fill=BOTH, expand=True)
        yScroll.pack(side=RIGHT, fill=Y)

        self.__createButtons()
        self.__createFileListBox()
        self.__createOptions()

        self.ButtonFrame.pack(side=BOTTOM, fill=None)
        self.LabelcpuString = StringVar()
        self.Labelcpu = tk.Label(self.ButtonFrame, textvariable=self.LabelcpuString)
        self.Labelcpu.pack(side=RIGHT)
        self.FileListBoxFrame.pack(side=TOP, fill=BOTH, expand=True)
        self.RightOptionsFrame.pack(side=RIGHT, fill=Y)
        self.LoggingFrame.pack(side=BOTTOM, fill=BOTH)
        self.logger = logging.getLogger(__package__)
        self.handler = GUIHandler(self.Logging)
        self.handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.logger.addHandler(self.handler)
        self.logger.setLevel(getattr(logging, self.opts.loglevel.upper()))
        self.logger.info("Config file:%s", self.opts.configfile)
        self.handler.flush()
        self.updateFileListBox()
        self.UpdateLabelcpu()

    def __createButtons(self):

        buttons = [{'name': 'Quit', 'text': 'QUIT', 'command': 'Quit', 'side': BOTTOM},
                   {'name': 'SpawnInputDialog', 'text': 'Add Input Files', 'side': LEFT},
                   {'name': 'RemoveFile', 'text': 'Remove Files', 'side': LEFT},
                   {'name': 'SpawnOutputDialog', 'text': 'Choose Output Directory', 'side': LEFT},
                   {'name': 'Settings', 'text': 'Settings', 'side': LEFT},
                   {'name': 'Parse', 'text': 'Parse!', 'side': LEFT}
                   ]

        for b in buttons:
            button = tk.Button(self.ButtonFrame)
            button.config(text=b['text'], command=getattr(self, b['name']+'Click'))
            button.pack(side=b['side'])
            setattr(self, 'Button'+b['name'], button)

    def __createOptions(self):
        '''Create the widgets for options and use setattr to assign them to self.'''
        self.checks = [{'name': 'plot', 'text': 'Plot', 'row': 1,
                        'tooltip': "Show summary plots after parsing."},
                       {'name': 'write', 'text': 'Write', 'row': 2,
                        'tooltip': "Write results to text files after parsing."},
                       {'name': 'SLM', 'text': 'Do SLM', 'row': 3,
                        'tooltip': "Compute SLM fits and write them to disk."},
                       {'name': 'skipohmic', 'text': 'Skip bad dJ/dV', 'row': 4,
                        'tooltip': 'Skip plots with d2J/dV2 < 0 between Vcutoff and Vmin/Vmax.'},
                       {'name': 'logr', 'text': 'Use log|R|', 'row': 5,
                        'tooltip': 'Use log |R| instead of |R| when computing histograms.'},
                       {'name': 'lorenzian', 'text': 'Lorenzian', 'row': 6,
                        'tooltip': 'Fit a Lorenzian instead of a Gaussian.'},
                       {'name': 'tracebyfile', 'text': 'AFM Data', 'row': 7,
                        'tooltip': 'Each file contains one (foward/backward) trace.'},
                       {'name': 'nolag', 'text': 'No lag plots.', 'row': 8,
                        'tooltip': 'Some data (often AFM) produces bad lag plots very slowly.'},
                       {'name': 'force', 'text': 'Force parse', 'row': 9,
                        'tooltip': 'Force parsing even when stats cannot be computed.'}
                       ]

        for _c in self.checks:
            setattr(self, _c['name'], IntVar())
            check = tk.Checkbutton(self.RightOptionsFrame, text=_c['text'],
                                   variable=getattr(self, _c['name']), command=self.checkOptions)
            check.grid(column=0, row=_c['row'], sticky=W)
            createToolTip(check, _c['tooltip'])
            setattr(self, 'Check_'+_c['name'], check)
            if getattr(self.opts, _c['name']):
                getattr(self, _c['name']).set(1)

        rowidx = len(self.checks)+1

        tk.Label(self.RightOptionsFrame, text="Output file base name:").grid(
            column=0, row=rowidx)
        self.OutputFileName = tk.Entry(self.RightOptionsFrame, width=20,
                                       font=Font(size=8, slant='italic'))
        for n in ('<Return>', '<Leave>', '<Enter>'):
            self.OutputFileName.bind(n, self.checkOutputFileName)
        self.OutputFileName.grid(column=0, row=rowidx+1)

        if self.opts.outfile:
            self.OutputFileName.insert(0, self.opts.outfile)

        tk.Label(self.RightOptionsFrame, text="Data to plot:").grid(
            column=0, row=rowidx+2, sticky=W)
        self.OptionsPlotsString = StringVar()
        self.OptionsPlotsString.set(self.opts.plots)
        self.OptionPlots = tk.OptionMenu(self.RightOptionsFrame,
                                         self.OptionsPlotsString, self.opts.plots, 'J', 'R', 'FN', 'SLM',
                                         command=self.checkOptions)
        self.OptionPlots.grid(column=0, row=rowidx+3, sticky=W)

        tk.Label(self.RightOptionsFrame, text="Histogram to plot:").grid(
            column=0, row=rowidx+4, sticky=W)
        self.OptionsHistPlotsString = StringVar()
        self.OptionsHistPlotsString.set(self.opts.histplots)
        self.OptionHistPlots = tk.OptionMenu(self.RightOptionsFrame,
                                             self.OptionsHistPlotsString, self.opts.histplots, 'NDC', 'G',
                                             command=self.checkOptions)
        createToolTip(self.OptionHistPlots, '(Normalized) Differential Conductance')
        self.OptionHistPlots.grid(column=0, row=rowidx+5, sticky=W)

        tk.Label(self.RightOptionsFrame, text="Derivative for heatmap:").grid(
            column=0, row=rowidx+6, sticky=W)
        self.OptionsHeatmapdString = StringVar()
        self.OptionsHeatmapdString.set(self.opts.heatmapd)
        self.OptionHeatmapd = tk.OptionMenu(self.RightOptionsFrame,
                                            self.OptionsHeatmapdString, self.opts.heatmapd, '0', '1', '2',
                                            command=self.checkOptions)
        createToolTip(self.OptionHeatmapd, '0 is equivalent to LogJ, 1 is Log|dJ/dV|, etc.')
        self.OptionHeatmapd.grid(column=0, row=rowidx+7, sticky=W)

        self.checkOptions()

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

    def ToFront(self):
        '''Try to bring the main window to the front on different platforms'''
        if platform.system() == "Darwin":
            os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        else:
            self.master.attributes('-topmost', 1)
            self.master.attributes('-topmost', 0)
        self.master.lift()

    def QuitClick(self):
        self.Quit()

    def Quit(self):
        self.master.destroy()

    def ParseClick(self):
        self.checkOptions()
        self.GUIParse()

    def SettingsClick(self):
        self.checkOptions()
        PreferencesWindow(self.master, self.opts)

    def RemoveFileClick(self):
        self.checkOptions()
        selected = [self.FileListBox.get(x) for x in self.FileListBox.curselection()]
        todel = []
        filelist = []
        for i in range(0, len(self.opts.in_files)):
            for s in selected:
                if self.opts.in_files[i].replace(" ", "_") == s:
                    todel.append(i)

        for i in range(0, len(self.opts.in_files)):
            if i not in todel:
                filelist.append(self.opts.in_files[i])
        self.opts.in_files = filelist
        self.updateFileListBox()
        self.FileListBox.selection_clear(0, END)

    def SpawnInputDialogClick(self):
        self.checkOptions()
        self.opts.in_files += filedialog.askopenfilename(
            title="Files to parse",
            multiple=True,
            initialdir=self.last_input_path,
            filetypes=[('Data files', '*_data.txt'),
                       ('Text files', '*.txt'), ('All files', '*')])
        if len(self.opts.in_files):
            self.last_input_path = os.path.split(self.opts.in_files[0])[0]
            if not self.outdir:
                self.opts.out_dir = self.last_input_path
            self.updateFileListBox()
            if not self.opts.outfile:
                self.OutputFileName.delete(0, END)
                self.OutputFileName.insert(0, os.path.basename(
                    self.opts.in_files[-1]).lower().replace('.txt', ''))
                self.checkOutputFileName()

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
        self.FileListBoxFrameLabelVar.set("Output to: %s/%s_*.txt" % (self.opts.out_dir, self.opts.outfile))

    def UpdateLabelcpu(self):
        if HASPSUTIL:
            self.LabelcpuString.set("CPU: %0.1f %%" % psutil.cpu_percent())
            self.Labelcpu.after(500, self.UpdateLabelcpu)

    def checkOptions(self, event=None):
        for c in self.checks:
            setattr(self.opts, c['name'], self.boolmap[getattr(self, c['name']).get()])

        if self.OptionsHistPlotsString.get() == 'NDC':
            self.OptionHeatmapd["state"] = DISABLED
            self.opts.heatmapd = 1
            self.OptionsHeatmapdString.set('1')
        else:
            self.OptionHeatmapd["state"] = NORMAL

        if not self.opts.write:
            self.opts.plot = True
            self.plot.set(1)
            self.Check_plot["state"] = DISABLED
        else:
            self.Check_plot["state"] = NORMAL

        self.opts.plots = self.OptionsPlotsString.get()
        self.opts.histplots = self.OptionsHistPlotsString.get()
        self.opts.heatmapd = int(self.OptionsHeatmapdString.get())

    def checkOutputFileName(self, event=None):
        self.opts.outfile = self.OutputFileName.get()
        self.UpdateFileListBoxFrameLabel()
