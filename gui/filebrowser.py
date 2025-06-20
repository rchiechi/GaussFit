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
import asyncio
import platform
import logging
import threading
import copy
import tkinter.ttk as tk
from tkinter import Toplevel
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
from queue import Queue, Empty
from gaussfit.logger import DelayedMultiprocessHandler
from gaussfit.logger import GUIHandler
from gaussfit.args import get_args, VERSION
# from gaussfit.parse import Parse
from gaussfit.parse import readfiles
from .libparse.ParseThread import ParseThread
from gaussfit.output.libwriter import doOutput
from gaussfit.output import Writer

try:
    import psutil
    HASPSUTIL = True
except ImportError:
    HASPSUTIL = False

absdir = os.path.dirname(os.path.realpath(__file__))

class ChooseFiles(tk.Frame):
    '''The main frame for adding/removing files, accessing settings
       and parsing.'''

    # from .libparse import GUIParse
    from .libparse import GUIPlot

    gothreads = []
    plots = []
    outdir = ''
    boolmap = {1: True, 0: False}
    lock = threading.Lock

    def __init__(self, master=None, **kwargs):
        super().__init__(master)
        if master is None:
            self.master = Toplevel() # Better, but STILL not the best approach
        else:
            self.master = master
        self.loop = kwargs.get('loop', asyncio.get_event_loop())
        self.logque = Queue(-1)
        self.bgimg = PhotoImage(file=os.path.join(absdir, 'RCCLabFluidic.png'))
        limg = Label(self.master, image=self.bgimg, background=GREY)
        limg.pack(side=TOP)
        try:
            self.last_input_path = os.getcwd()
        except KeyError:
            self.last_input_path = os.path.expanduser('~')
        self.opts_parser, self.opts = get_args()
        self.degfreedom = {'init': self.opts.degfree, 'user': self.opts.degfree}
        self.master.tk_setPalette(background=GREY, activeBackground=GREY)
        self.master.title(f"RCCLab EGaIn Data Parser v{VERSION}")
        self.master.geometry('800x1000+250+250')\

        self.pack(fill=BOTH)
        self.__createWidgets()
        self.ToFront()
        # self.master.mainloop()

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

    # def ParseClick(self):        
    #     self.checkOptions()
    #     self.GUIParse()

    def ParseClick(self):
        self.checkOptions()  # Assuming this is synchronous
        self.ButtonParse['state'] = DISABLED  # Disable the button immediately
        if self.gothreads and any(t.is_alive() for t in self.gothreads):
            self.logger.warning("Another parse is already running!") #Or display on GUI
            self.ButtonParse['state'] = NORMAL  # Re-enable if already running
            return

        self.preParse()
        # Schedule _async_GUIParse to run on the asyncio loop
        # self.loop.call_soon_threadsafe(self._async_GUIParse) #Key change
        self.GUIParse()
        # self.master.after(100, self.check_queue)  # Start checking the queue


    def GUIParse(self):
        # async def do_parse(): #Inner async function
        if self.opts.in_files:
            try:
                # df = await readfiles(self.opts.in_files) # Await here!
                # parser = Parse(df, handler=DelayedMultiprocessHandler(self.logque))
                thread = ParseThread(self.opts, self.logque)
                self.gothreads.append(thread)
                thread.start()
                self.master.after(100, self.check_queue)

            except Exception as e:
                logger.error(f"Unhandled exception in _async_GUIParse: {e}")
                import traceback
                traceback.print_exc()
                # Handle the exception (e.g., put an error message in the queue)
                self.master.after(0, self.postParse) #Safely call on main thread

        else:
            self.logger.warning("No input files!")
            self.handler.flush()  # No `await` needed
            self.master.after(0, self.postParse) #Safely call on main thread
        # asyncio.run_coroutine_threadsafe(do_parse(), self.loop) #Run do_parse on event loop

    def check_queue(self):
        if any(t.is_alive() for t in self.gothreads):
            while not self.logque.empty():
                self.handler.emit(self.logque.get_nowait())
            self.handler.flush()
            self.master.after(100, self.check_queue)
        else:
            self.master.after(0, self.postParse) #Safely call on main thread

    def preParse(self):
        '''We need to check a couple of things right before we start parsing'''
        self.degfreedom = self.opts.degfree
        if self.opts.degfree == 1 and len(self.opts.in_files) > 1:
            self.opts.degfree = len(self.opts.in_files) - 1

    def postParse(self):
        '''We need to check a couple of things right after we finish parsing'''
        self.ButtonParse['state'] = NORMAL
        self.opts.degfree = self.degfreedom
        self.logger.info("Parse complete!")
        try:
            gothread = self.gothreads.pop()
            if self.opts.write and not gothread.parser.error:
                writer = Writer(gothread.parser)
                doOutput(writer)
            if self.opts.plot and not gothread.parser.error:
                self.GUIPlot(parser=gothread.parser)
        except (IndexError, ValueError) as e:
            self.logger.error(e)

    def SettingsClick(self):
        """
        Opens a modal preferences window. If the user saves, it updates the
        canonical self.opts object in-place and then syncs the GUI.
        """
        # First, ensure self.opts is up-to-date with any last-minute GUI changes.
        self.checkOptions()
        dialog = PreferencesWindow(self.master, self.opts_parser, self.opts)
        self.master.wait_window(dialog)
    
        if dialog.saved:
            self._update_gui_from_opts()
            self.logger.info("Preferences have been updated and applied.")
        else:
            self.logger.info("Preference changes were discarded.")

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
        """
        Reads all values from the GUI widgets and updates the self.opts namespace.
        This is the primary public callback for all interactive widgets.
        It then applies any dependent UI logic.
        """
        # Part 1: Read all values from the GUI widgets into self.opts
        for c in self.checks:
            setattr(self.opts, c['name'], self.boolmap[getattr(self, c['name']).get()])
    
        self.opts.plots = self.OptionsPlotsString.get()
        self.opts.histplots = self.OptionsHistPlotsString.get()
        self.opts.heatmapd = int(self.OptionsHeatmapdString.get())
        self.opts.outfile = self.OutputFileName.get()
    
        # Part 2: Apply UI logic based on the new, authoritative state in self.opts
        if self.opts.histplots == 'NDC':
            self.OptionHeatmapd["state"] = DISABLED
            # This logic can override user settings for consistency.
            # We must update both the UI and the underlying opts object.
            self.opts.heatmapd = 1
            self.OptionsHeatmapdString.set('1')
        else:
            self.OptionHeatmapd["state"] = NORMAL
    
        if not self.opts.write:
            # If write is off, plot must be on. Enforce this in opts and the UI.
            self.opts.plot = True
            self.plot.set(1) # Also update the widget that gets disabled
            self.Check_plot["state"] = DISABLED
        else:
            self.Check_plot["state"] = NORMAL
    
        # Finally, update any dependent display labels
        self.UpdateFileListBoxFrameLabel()


    def _update_gui_from_opts(self):
        """
        (Internal Helper) Reads all values from self.opts and forces the GUI widgets to match.
        This should be called ONLY after self.opts has been updated from an
        external source, like the preferences dialog.
        """
        # Force all GUI widgets to match the state of self.opts.
        for c in self.checks:
            val = 1 if getattr(self.opts, c['name'], False) else 0
            getattr(self, c['name']).set(val)
    
        self.OptionsPlotsString.set(self.opts.plots)
        self.OptionsHistPlotsString.set(self.opts.histplots)
        self.OptionsHeatmapdString.set(str(self.opts.heatmapd))
        self.OutputFileName.delete(0, 'end')
        self.OutputFileName.insert(0, self.opts.outfile)
    
        # After forcing the GUI to match the new opts, we must call checkOptions()
        # to apply any dependent UI logic (like disabling widgets).
        self.checkOptions()



    def checkOutputFileName(self, event=None):
        self.opts.outfile = self.OutputFileName.get()
        self.UpdateFileListBoxFrameLabel()


    def run_async_task(self):
        async def wrapper():
            result = await self.other_class.async_task(self.update_gui)
            self.update_gui(f"Result: {result}")

        asyncio.run(wrapper())

    def update_gui(self, message):
        #Safely update GUI from within async task
        self.root.after(0, lambda: self.label.config(text=message))