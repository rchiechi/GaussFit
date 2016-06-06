#!/usr/bin/python3
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
import sys,os,platform,threading,logging
from gaussfit import *
from gaussfit.Output import Writer,Plotter
from gaussfit.logger import GUIHandler
from tkinter import filedialog #some weird bug...
from tkinter import *
from tkinter.ttk import *
from tkinter.font import *
from gui.colors import *
from gui.tooltip import *



class ParseThread(threading.Thread):
    def __init__(self,opts,parser):
        threading.Thread.__init__(self)
        self.parser = parser
        self.opts = opts
    def run(self):
        self.parser.ReadFiles(self.opts.in_files)


class ChooseFiles(Frame):

    def __init__(self, opts, master=None):
        Frame.__init__(self, master)
        self.boolmap = {1:True, 0:False}
        try:
            self.last_input_path = os.getcwd()
        except KeyError:
            self.last_input_path = os.path.expanduser('~')
        self.opts = opts
        self.gothreads = []
        
        self.master.tk_setPalette(background=GREY,
            activeBackground=GREY)
        self.master.title("RCCLab EGaIn Data Parser")
        self.master.geometry('850x750+250-250')
        self.pack(fill=BOTH)
        self.__createWidgets()
        self.ToFront()
        self.mainloop()

    def __createWidgets(self):
            
        self.ButtonFrame = Frame(self)
        self.LeftOptionsFrame = Frame(self)
        self.LoggingFrame = Frame(self)
        self.RightOptionsFrame = Frame(self)


        self.FileListBoxFrame = Frame(self)
    
        yScroll = Scrollbar(self.LoggingFrame, orient=VERTICAL)
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
        self.LeftOptionsFrame.pack(side=LEFT,fill=Y)
        self.RightOptionsFrame.pack(side=RIGHT,fill=Y)
        self.LoggingFrame.pack(side=BOTTOM, fill=BOTH)
        self.logger = logging.getLogger('parser.gui')
        self.handler = GUIHandler(self.Logging)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(getattr(logging,self.opts.loglevel.upper()))
        self.logger.info("Logging...")
        self.handler.flush()
        self.updateFileListBox()

    def __createButtons(self):
            
        buttons = [
               {'name':'Quit','text':'QUIT','command':'Quit','side':BOTTOM},
               {'name':'SpawnInputDialog','text':'Add Input Files','side':LEFT},
               {'name':'RemoveFile','text':'Remove Files','side':LEFT},
               {'name':'SpawnOutputDialog','text':'Choose Output Directory','side':LEFT},
               {'name':'Parse','text':'Parse!','side':LEFT}
               ]

        for b in buttons:
            button = Button(self.ButtonFrame)
            button.config(text=b['text'],command=getattr(self,b['name']+'Click'))
            button.pack(side=b['side'])
            setattr(self,'Button'+b['name'],button)

        
    def __createOptions(self):
        self.checks = [
              {'name':'plot','text':'Plot','row':1,'tooltip':"Show summary plots after parsing."},
              {'name':'write','text':'Write','row':2,'tooltip':"Write results to text files after parsing."},
              {'name':'skipohmic','text':'Skip bad dJ/dV','row':3,'tooltip':'Skip plots with d2J/dV2 < 0 between Vcutoff and Vmin/Vmax.'},
              {'name':'interpolateminfn','text':'Interpolate FN','row':4,'tooltip':'Find FN from the minimum of the derivative.'},
              {'name':'logr','text':'Use log|R|','row':5,'tooltip':'Use log |R| instead of |R| when computing histograms.'},
              {'name':'lorenzian','text':'Lorenzian','row':6,'tooltip':'Fit a Lorenzian instead of a Gaussian.'},
              {'name':'tracebyfile','text':'AFM Data','row':7,'tooltip':'Each file contains one (foward/backward) trace.'}
              ]

        for c in self.checks:
            setattr(self,c['name'],IntVar())
            check = Checkbutton(self.RightOptionsFrame, text=c['text'],
                    variable=getattr(self,c['name']), command=self.checkOptions)
            check.grid(column=0,row=c['row'],sticky=W)
            createToolTip(check,c['tooltip'])
            setattr(self,'Check_'+c['name'],check)
            if getattr(self.opts,c['name']):
                getattr(self,c['name']).set(1)
        
        rowidx = len(self.checks)+1

        Label(self.RightOptionsFrame, text="Output file base name:").grid(column=0,row=rowidx)
        self.OutputFileName = Entry(self.RightOptionsFrame, width=20,
                font=Font(size=8,slant='italic'))
        for n in ('<Return>','<Leave>','<Enter>'):
            self.OutputFileName.bind(n, self.checkOutputFileName)
        self.OutputFileName.grid(column=0,row=rowidx+1)
        
        if self.opts.outfile:
            self.OutputFileName.insert(0,self.opts.outfile)
    
        Label(self.RightOptionsFrame, text="What to plot:").grid(column=0,row=rowidx+2,sticky=W)
        self.OptionsPlotsString = StringVar()
        self.OptionsPlotsString.set(','.join(self.opts.plots))
        self.OptionPlots = OptionMenu(self.RightOptionsFrame, self.OptionsPlotsString,'J','R',
                                 command=self.checkOptions)
        self.OptionPlots.grid(column=0,row=rowidx+3,sticky=W)
    
        Label(self.RightOptionsFrame, text="Derivative for heatmap:").grid(column=0,row=rowidx+4,sticky=W)
        self.OptionsHeatmapdString = StringVar()
        self.OptionsHeatmapdString.set(self.opts.heatmapd)
        self.OptionHeatmapd = OptionMenu(self.RightOptionsFrame, self.OptionsHeatmapdString,'0','1','2',
                                 command=self.checkOptions)
        self.OptionHeatmapd.grid(column=0,row=rowidx+5,sticky=W)

        lbls = [
            {'name': 'Columns', 'text': "Columns to parse:",
             'tooltip': 'Columns from input data to parse as X/Y data.'},
            {'name': 'Vcutoff', 'text': "Cuttoff for d2J/dV2:",
             'tooltip': "Check the values of d2J/dV2 between |vcutoff| and Vmin/Vmax for line-shape filtering. Set to -1 for Vmin/Vmax."},
            {'name': 'Smooth', 'text': "Smoothing parameter:",
             'tooltip': "The cutoff value for the residual squares (the difference between experimental data points and the fit). The default is 1e-12. Set to 0 to disable smoothing."},
            {'name': 'Gminmax', 'text': "Y-scale for histogram:",
                 'tooltip': "Set Ymin,Ymax for the heapmap plot (lower-left of plot output)."},
            {'name': 'Bins', 'text': "Bins for J/R Histograms:",
             'tooltip': "Set binning for histograms of J and R."},
            {'name': 'Heatmapbins', 'text': "Bins for G Histograms:",
             'tooltip': "Set binning for heatmap histograms."}
            ]
        
        i = 0
        for l in lbls:
            Label(self.LeftOptionsFrame, text=l['text']).grid(column=0,row=i)
            entry = Entry(self.LeftOptionsFrame, width=8, 
                    validate='focus', validatecommand=self.checkOptions)
            entry.bind("<Return>", self.checkOptions)
            entry.bind("<Leave>", self.checkOptions)
            entry.bind("<Enter>", self.checkOptions)
            entry.grid(column=0,row=i+1)
            createToolTip(entry, l['tooltip'])
            setattr(self, 'Entry'+l['name'], entry)
            i+=2

        self.checkOptions()

    def __createFileListBox(self):
        self.FileListBoxFrameLabelVar = StringVar()
        self.FileListBoxFrameLabel = Label(self.FileListBoxFrame,\
                textvariable=self.FileListBoxFrameLabelVar,\
                font=Font(size=10,weight="bold"))
        self.FileListBoxFrameLabel.pack(side=TOP,fill=X)



        yScroll = Scrollbar(self.FileListBoxFrame,orient=VERTICAL)
        yScroll.pack(side=RIGHT,fill=Y)
        xScroll = Scrollbar(self.FileListBoxFrame, orient=HORIZONTAL)
        xScroll.pack(side=BOTTOM,fill=X)
        self.filelist = StringVar()
        self.FileListBox = Listbox(self.FileListBoxFrame, listvariable=self.filelist, selectmode=EXTENDED, 
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
        self.Parse()

    def Parse(self):
        if len(self.gothreads):
            self.ButtonParse['state']=DISABLED
            for c in self.gothreads:
                if c.is_alive():
                    self.ButtonParse.after('500',self.Parse)
                    return  
            self.logger.info("Parse complete!")
            gothread = self.gothreads.pop()
            self.ButtonParse['state']=NORMAL
            if self.opts.write and not gothread.parser.error: 
                writer = Writer(gothread.parser)
                Parse.doOutput(writer)
            if self.opts.plot and not gothread.parser.error:
                plotter = Plotter(gothread.parser)
                self.logger.info("Generating plots...")
                try:
                        import matplotlib.pyplot as plt
                        plotter.DoPlots(plt,*self.opts.plots)
                        plt.show(block=False)
                except ImportError as msg:
                        self.logger.error("Cannot import matplotlib! %s", str(msg), exc_info=False)
        else: 
            if len(self.opts.in_files):
                parser = Parse(self.opts,handler=self.handler)
                self.gothreads.append(ParseThread(self.opts,parser))
                self.gothreads[-1].start()
                self.ButtonParse.after('500',self.Parse)
            else:
                self.logger.warn("No input files!")
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
        self.opts.in_files += filedialog.askopenfilename(title="Files to parse", multiple=True, \
                        initialdir=self.last_input_path,\
                         filetypes=[('Data files','*_data.txt'),('Text files','*.txt'),('All files', '*')])
        if len(self.opts.in_files):
            self.last_input_path = os.path.split(self.opts.in_files[0])[0]
            self.updateFileListBox()
            if not self.opts.outfile:
                self.OutputFileName.delete(0,END)
                self.OutputFileName.insert(0,\
                        os.path.basename(self.opts.in_files[-1]).lower().replace('.txt',''))
                self.checkOutputFileName(None)

    def updateFileListBox(self):
        self.filelist.set(" ".join([x.replace(" ","_") for x in self.opts.in_files]))
            
    def SpawnOutputDialogClick(self):
        outdir = filedialog.askdirectory(title="Select Output File(s)", initialdir=self.opts.out_dir)
        if os.path.exists(outdir):
            self.opts.out_dir = outdir
            self.UpdateFileListBoxFrameLabel()

    def UpdateFileListBoxFrameLabel(self):
        self.FileListBoxFrameLabelVar.set("Output to: %s/%s_*.txt"% (self.opts.out_dir, self.opts.outfile) )

    def checkOptions(self, event=None):
        for c in self.checks:
            setattr(self.opts,c['name'],self.boolmap[getattr(self,c['name']).get()])
            
        if not self.opts.write:
            self.opts.plot = True
            self.plot.set(1)
            self.Check_plot["state"]=DISABLED
        else:
            self.Check_plot["state"]=NORMAL
        
        for n in (('Smooth',1e-12),('Bins',50),('Heatmapbins',25)):
            try:
                var = getattr(self,'Entry'+n[0]).get()
                if 'Smooth' in n:
                    var = float(var)
                else:
                    var = int(var)
            except ValueError as msg:
                var = n[1]
            if var == 0:
                var = int(0)
            if 'Smooth' in n:
                if var >= 0:
                    setattr(self.opts,n[0].lower(),var)
            elif var > 0:
                setattr(self.opts,n[0].lower(),var)
            getattr(self,'Entry'+n[0]).delete(0,END)
            getattr(self,'Entry'+n[0]).insert(0,getattr(self.opts,n[0].lower()))
    
        self.opts.plots = self.OptionsPlotsString.get().split(',')
        self.opts.heatmapd = int(self.OptionsHeatmapdString.get())
        self.checkGminmaxEntry(event)
        self.checkOutputFileName(event) 
        self.checkColumnEntry(event)
        self.checkVcutoff(event)

    def checkVcutoff(self,event):
        try:
            vcutoff = float(self.EntryVcutoff.get())
            if vcutoff != -1:
                vcutoff = abs(vcutoff)
            self.opts.vcutoff = vcutoff  
        except ValueError:
            self.opts.vcutoff = -1

        self.EntryVcutoff.delete(0, END)
        if self.opts.vcutoff > 0:
            self.EntryVcutoff.insert(0,self.opts.vcutoff)
        else:
            self.EntryVcutoff.insert(0,"Vmax")


    def checkOutputFileName(self, event):
        self.opts.outfile = self.OutputFileName.get()
        self.UpdateFileListBoxFrameLabel()

    def checkColumnEntry(self, event):
        try:
            x, y = self.EntryColumns.get().split(",")
            self.opts.Xcol, self.opts.Ycol = int(x)-1, int(y)-1
        except ValueError as msg:
            pass
        self.EntryColumns.delete(0, END)
        self.EntryColumns.insert(0, ",".join(( str(self.opts.Xcol+1), str(self.opts.Ycol+1) )))

    def checkGminmaxEntry(self, event):
        try:
            x, y = self.EntryGminmax.get().split(",")
            self.opts.mlow, self.opts.mhi = int(x), int(y)
        except ValueError as msg:
            pass
        self.EntryGminmax.delete(0, END)
        self.EntryGminmax.insert(0, ",".join( (str(self.opts.mlow), str(self.opts.mhi)) ))


