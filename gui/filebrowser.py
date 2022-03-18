'''
Copyright (C) 2018 Ryan Chiechi <r.c.chiechi@rug.nl>
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
import threading
import logging

import tkinter.ttk as tk
from tkinter import filedialog
from tkinter import Text,IntVar,StringVar,Listbox
from tkinter import N,S,E,W,X,Y # pylint: disable=unused-import
from tkinter import TOP,BOTTOM,LEFT,RIGHT # pylint: disable=unused-import
from tkinter import END,BOTH,VERTICAL,HORIZONTAL # pylint: disable=unused-import
from tkinter import EXTENDED,RAISED,DISABLED,NORMAL # pylint: disable=unused-import
from tkinter.font import Font

from gui.prefs import PreferencesWindow
from gui.colors import BLACK,YELLOW,WHITE,RED,TEAL,GREEN,BLUE,GREY # pylint: disable=unused-import
from gui.tooltip import createToolTip
from gaussfit import Parse
from gaussfit.output import Writer,Plotter
from gaussfit.logger import GUIHandler

# pylint: disable=attribute-defined-outside-init,missing-function-docstring

class ParseThread(threading.Thread):
    '''A Thread object to run the parse in so it
       doesn't block the main GUI thread.'''
    def __init__(self,opts,parser):
        threading.Thread.__init__(self)
        self.parser = parser
        self.opts = opts
    def run(self):
        self.parser.readfiles(self.opts.in_files)


class ChooseFiles(tk.Frame):
    '''The main frame for adding/removing files, accessing setttings
       and parsing.'''
    def __init__(self, opts, master=None):
        tk.Frame.__init__(self, master)
        self.boolmap = {1:True, 0:False}
        try:
            self.last_input_path = os.getcwd()
        except KeyError:
            self.last_input_path = os.path.expanduser('~')
        self.opts = opts
        self.gothreads = []
        self.checks = []
        self.degfreedom = {'init': self.opts.degfree, 'user': self.opts.degfree}
        self.outdir = ''
        self.master.tk_setPalette(background=GREY,
            activeBackground=GREY)
        self.master.title("RCCLab EGaIn Data Parser")
        self.master.geometry('700x800+250-250')
        self.pack(fill=BOTH)
        self.__createWidgets()
        self.ToFront()
        self.mainloop()

    def __createWidgets(self):

        self.ButtonFrame = tk.Frame(self)
        # self.LeftOptionsFrame = tk.Frame(self)
        self.LoggingFrame = tk.Frame(self)
        self.RightOptionsFrame = tk.Frame(self)


        self.FileListBoxFrame = tk.Frame(self)

        yScroll = tk.Scrollbar(self.LoggingFrame, orient=VERTICAL)
        self.Logging = Text(self.LoggingFrame, height=20, width=0,
                bg=BLACK, fg=WHITE, yscrollcommand=yScroll.set)
        yScroll['command'] = self.Logging.yview
        self.Logging.pack(side=LEFT,fill=BOTH,expand=True)
        yScroll.pack(side=RIGHT,fill=Y)

        self.__createButtons()
        self.__createFileListBox()
        self.__createOptions()

        self.ButtonFrame.pack(side=BOTTOM,fill=None)
        self.FileListBoxFrame.pack(side=TOP, fil=BOTH, expand=True)
        # self.LeftOptionsFrame.pack(side=LEFT,fill=Y)
        self.RightOptionsFrame.pack(side=RIGHT,fill=Y)
        self.LoggingFrame.pack(side=BOTTOM, fill=BOTH)
        self.logger = logging.getLogger('parser.gui')
        self.handler = GUIHandler(self.Logging)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(getattr(logging,self.opts.loglevel.upper()))
        self.logger.info("Config file:%s", self.opts.configfile)
        self.handler.flush()
        self.updateFileListBox()

    def __createButtons(self):

        buttons = [
            {'name':'Quit','text':'QUIT','command':'Quit','side':BOTTOM},
            {'name':'SpawnInputDialog','text':'Add Input Files','side':LEFT},
            {'name':'RemoveFile','text':'Remove Files','side':LEFT},
            {'name':'SpawnOutputDialog','text':'Choose Output Directory','side':LEFT},
            {'name':'Settings','text':'Settings','side':LEFT},
            {'name':'Parse','text':'Parse!','side':LEFT}
              ]

        for b in buttons:
            button = tk.Button(self.ButtonFrame)
            button.config(text=b['text'],command=getattr(self,b['name']+'Click'))
            button.pack(side=b['side'])
            setattr(self,'Button'+b['name'],button)


    def __createOptions(self):
        '''Create the widgets for options and use setattr to assign them to self.'''
        self.checks = [
              {'name':'plot','text':'Plot','row':1,
                'tooltip':"Show summary plots after parsing."},
              {'name':'write','text':'Write','row':2,
                'tooltip':"Write results to text files after parsing."},
              {'name':'skipohmic','text':'Skip bad dJ/dV','row':3,
                'tooltip':'Skip plots with d2J/dV2 < 0 between Vcutoff and Vmin/Vmax.'},
              {'name':'interpolateminfn','text':'Interpolate FN','row':4,
                'tooltip':'Find FN from the minimum of the derivative.'},
              {'name':'logr','text':'Use log|R|','row':5,
                'tooltip':'Use log |R| instead of |R| when computing histograms.'},
              {'name':'lorenzian','text':'Lorenzian','row':6,
                'tooltip':'Fit a Lorenzian instead of a Gaussian.'},
              {'name':'tracebyfile','text':'AFM Data','row':7,
                'tooltip':'Each file contains one (foward/backward) trace.'},
              {'name':'nolag','text':'No lag plots.','row':8,
                'tooltip':'Some data (often AFM) produces bad lag plots very slowly.'}
              ]

        for _c in self.checks:
            setattr(self,_c['name'],IntVar())
            check = tk.Checkbutton(self.RightOptionsFrame, text=_c['text'],
                    variable=getattr(self,_c['name']), command=self.checkOptions)
            check.grid(column=0,row=_c['row'],sticky=W)
            createToolTip(check,_c['tooltip'])
            setattr(self,'Check_'+_c['name'],check)
            if getattr(self.opts,_c['name']):
                getattr(self,_c['name']).set(1)

        rowidx = len(self.checks)+1

        tk.Label(self.RightOptionsFrame, text="Output file base name:").grid(
            column=0,row=rowidx)
        self.OutputFileName = tk.Entry(self.RightOptionsFrame, width=20,
                font=Font(size=8,slant='italic'))
        for n in ('<Return>','<Leave>','<Enter>'):
            self.OutputFileName.bind(n, self.checkOutputFileName)
        self.OutputFileName.grid(column=0,row=rowidx+1)

        if self.opts.outfile:
            self.OutputFileName.insert(0,self.opts.outfile)

        tk.Label(self.RightOptionsFrame, text="Data to plot:").grid(
            column=0,row=rowidx+2,sticky=W)
        self.OptionsPlotsString = StringVar()
        self.OptionsPlotsString.set(self.opts.plots)
        self.OptionPlots = tk.OptionMenu(self.RightOptionsFrame,
            self.OptionsPlotsString,'J','R',
            command=self.checkOptions)
        self.OptionPlots.grid(column=0,row=rowidx+3,sticky=W)

        tk.Label(self.RightOptionsFrame, text="Histogram to plot:").grid(
            column=0,row=rowidx+4,sticky=W)
        self.OptionsHistPlotsString = StringVar()
        self.OptionsHistPlotsString.set(self.opts.histplots)
        self.OptionHistPlots = tk.OptionMenu(self.RightOptionsFrame,
            self.OptionsHistPlotsString,'NDC','G',
            command=self.checkOptions)
        self.OptionHistPlots.grid(column=0,row=rowidx+5,sticky=W)


        tk.Label(self.RightOptionsFrame, text="Derivative for heatmap:").grid(
            column=0,row=rowidx+6,sticky=W)
        self.OptionsHeatmapdString = StringVar()
        self.OptionsHeatmapdString.set(self.opts.heatmapd)
        self.OptionHeatmapd = tk.OptionMenu(self.RightOptionsFrame,
            self.OptionsHeatmapdString,'0','1','2',
            command=self.checkOptions)
        self.OptionHeatmapd.grid(column=0,row=rowidx+7,sticky=W)

        self.checkOptions()

    def __createFileListBox(self):
        self.FileListBoxFrameLabelVar = StringVar()
        self.FileListBoxFrameLabel = tk.Label(self.FileListBoxFrame,\
                textvariable=self.FileListBoxFrameLabelVar,\
                font=Font(size=10,weight="bold"))
        self.FileListBoxFrameLabel.pack(side=TOP,fill=X)



        yScroll = tk.Scrollbar(self.FileListBoxFrame,orient=VERTICAL)
        yScroll.pack(side=RIGHT,fill=Y)
        xScroll = tk.Scrollbar(self.FileListBoxFrame, orient=HORIZONTAL)
        xScroll.pack(side=BOTTOM,fill=X)
        self.filelist = StringVar()
        self.FileListBox = Listbox(self.FileListBoxFrame,
            listvariable=self.filelist, selectmode=EXTENDED,
            height = 20, width = 0, relief=RAISED, bd=1,
            bg=WHITE,
            #font=Font(size=10),
            xscrollcommand=xScroll.set,
            yscrollcommand=yScroll.set)
        xScroll['command'] = self.FileListBox.xview
        yScroll['command'] = self.FileListBox.yview
        self.FileListBox.pack(side=LEFT, fill=BOTH, expand=True)
        self.UpdateFileListBoxFrameLabel()

    def ToFront(self):
        '''Try to bring hte main window to the front on different platforms'''
        if platform.system() == "Darwin":
            os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')#pylint: disable=line-too-long
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
        self.Parse()

    def SettingsClick(self):
        self.checkOptions()
        PreferencesWindow(self.master, self.opts)

    def Parse(self):
        '''Start the parser in a background thread and wait for it to complete.'''
        self.preParse()
        if self.gothreads:
            self.ButtonParse['state']=DISABLED #pylint: disable=E1101
            for _t in self.gothreads:
                if _t.is_alive():
                    self.ButtonParse.after('500',self.Parse) #pylint: disable=E1101
                    return
            self.logger.info("Parse complete!")
            gothread = self.gothreads.pop()
            self.ButtonParse['state']=NORMAL #pylint: disable=E1101
            if self.opts.write and not gothread.parser.error:
                writer = Writer(gothread.parser)
                Parse.doOutput(writer)
            if self.opts.plot and not gothread.parser.error:
                try:
                    import matplotlib.pyplot as plt #pylint: disable=C0415
                    plotter = Plotter(gothread.parser, plt)
                    self.logger.info("Generating plots...")
                    plotter.DoPlots()
                    plt.show(block=False)
                except ImportError as msg:
                    self.logger.error("Cannot import matplotlib! %s", str(msg), exc_info=False)
        else:
            if self.opts.in_files:
                parser = Parse(self.opts,handler=self.handler)
                self.gothreads.append(ParseThread(self.opts,parser))
                self.gothreads[-1].start()
                self.ButtonParse.after('500',self.Parse) #pylint: disable=E1101
            else:
                self.logger.warning("No input files!")
        self.postParse()
        self.handler.flush()

    def RemoveFileClick(self):
        self.checkOptions()
        selected = [self.FileListBox.get(x) for x in self.FileListBox.curselection()]
        todel = []
        filelist = []
        for i in range(0, len(self.opts.in_files)):
            for s in selected:
                if self.opts.in_files[i].replace(" ","_") == s:
                    todel.append(i)

        for i in range(0, len(self.opts.in_files)):
            if i not in todel:
                filelist.append(self.opts.in_files[i])
        self.opts.in_files = filelist
        self.updateFileListBox()
        self.FileListBox.selection_clear(0,END)


    def SpawnInputDialogClick(self):
        self.checkOptions()
        self.opts.in_files += filedialog.askopenfilename(
            title="Files to parse",
            multiple=True,
            initialdir=self.last_input_path,
            filetypes=[('Data files','*_data.txt'),
                        ('Text files','*.txt'),('All files', '*')])
        if len(self.opts.in_files):
            self.last_input_path = os.path.split(self.opts.in_files[0])[0]
            if not self.outdir:
                self.opts.out_dir = self.last_input_path
            self.updateFileListBox()
            if not self.opts.outfile:
                self.OutputFileName.delete(0,END)
                self.OutputFileName.insert(0,\
                        os.path.basename(self.opts.in_files[-1]).lower().replace('.txt',''))
                self.checkOutputFileName()

    def updateFileListBox(self):
        self.filelist.set(" ".join([x.replace(" ","_") for x in self.opts.in_files]))

    def SpawnOutputDialogClick(self):
        outdir = filedialog.askdirectory(title="Select Output File(s)",
                                         initialdir=self.opts.out_dir)
        if os.path.exists(outdir):
            self.opts.out_dir = outdir
            self.outdir = self.opts.out_dir  # So we know the user set the output dir
            self.UpdateFileListBoxFrameLabel()

    def UpdateFileListBoxFrameLabel(self):
        self.FileListBoxFrameLabelVar.set(
            "Output to: %s/%s_*.txt"% (self.opts.out_dir, self.opts.outfile) )

    def checkOptions(self, event=None): #pylint: disable=W0613
        for c in self.checks:
            setattr(self.opts,c['name'],self.boolmap[getattr(self,c['name']).get()])

        if self.OptionsHistPlotsString.get() == 'NDC':
            self.OptionHeatmapd["state"]=DISABLED
            self.opts.heatmapd=1
            self.OptionsHeatmapdString.set('1')
        else:
            self.OptionHeatmapd["state"]=NORMAL

        if not self.opts.write:
            self.opts.plot = True
            self.plot.set(1) #pylint: disable=E1101
            self.Check_plot["state"]=DISABLED #pylint: disable=E1101
        else:
            self.Check_plot["state"]=NORMAL #pylint: disable=E1101

        self.opts.plots = self.OptionsPlotsString.get()
        self.opts.histplots = self.OptionsHistPlotsString.get()
        self.opts.heatmapd = int(self.OptionsHeatmapdString.get())

    def checkOutputFileName(self, event=None): # pylint: disable=unused-argument
        self.opts.outfile = self.OutputFileName.get()
        self.UpdateFileListBoxFrameLabel()

    def preParse(self):
        '''We need to check a couple of things right before we start parsing'''
        self.degfreedom = self.opts.degfree # Store self.opts.degfree before parsing
        if self.opts.degfree == 0 and len(self.opts.in_files):
            self.opts.degfree = len(self.opts.in_files)-1

    def postParse(self):
        '''We need to check a couple of things right after we finish parsing'''
        self.opts.degfree = self.degfreedom # Retore self.opts.degfree after parsing