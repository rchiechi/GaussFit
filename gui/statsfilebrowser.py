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
from gaussfit.stats import Stats
from gaussfit.Output import WriteStats,StatPlotter
from tkinter import filedialog #some weird bug...
from tkinter import *
from tkinter.ttk import *
from tkinter.font import *
from gui.tooltip import *
from gui.colors import *

class ParseThread(threading.Thread):
    def __init__(self,statparser):
        threading.Thread.__init__(self)
        self.statparser = statparser 
    def run(self):
        self.statparser.parse()
        return


class ChooseFiles(Frame):

    def __init__(self, opts, master=None):
        Frame.__init__(self, master)
        self.boolmap = {1:True, 0:False}
        try:
            self.defaultdir = os.getcwd()
        except KeyError:
            self.defaultdir = os.path.expanduser('~')
        self.last_input_pathA = self.defaultdir
        self.last_input_pathB = self.defaultdir
        self.opts = opts
        self.gothreads = []
        self.master.tk_setPalette(background=GREY,
                activeBackground=GREY)
        self.master.title("RCCLab EGaIn Statistics")
        self.master.geometry('1100x700+250-250')
        self.pack(fill=BOTH)
        self.__createWidgets()
        self.ToFront()
        self.mainloop()

    def __createWidgets(self):
            
        self.FileListBoxFrame = Frame(self)
        self.LeftOptionsFrame = Frame(self)
        self.LoggingFrame = Frame(self)
        self.RightOptionsFrame = Frame(self)
        self.ButtonFrame = Frame(self)
        
        yScroll = Scrollbar(self.LoggingFrame, orient=VERTICAL)
        self.Logging = Text(self.LoggingFrame, height=20, width=0, 
                bg=BLACK, fg=WHITE, yscrollcommand=yScroll.set)
        yScroll['command'] = self.Logging.yview
        self.Logging.pack(side=LEFT,fill=BOTH,expand=True)
        yScroll.pack(side=RIGHT,fill=Y)


        self.__createButtons()
        self.__createFileListBox()
        self.__createOptions()

        self.ButtonFrame.pack(side=BOTTOM)
        self.FileListBoxFrame.pack(side=TOP, fil=BOTH, expand=1)
        self.LeftOptionsFrame.pack(side=LEFT,fill=Y)
        self.RightOptionsFrame.pack(side=RIGHT,fill=Y)
        self.LoggingFrame.pack(side=BOTTOM, fill=BOTH, expand=1)
        
        self.logger = logging.getLogger('gui')
        self.logger.addHandler(LoggingToGUI(self.Logging))
        self.logger.info("Logging...")

        self.updateFileListBox('A')
        self.updateFileListBox('B')

    def __createButtons(self):
            
        buttons = [
               {'name':'Quit','text':'QUIT','command':'Quit','side':BOTTOM},
               {'name':'SpawnInputDialogA','text':'Add Input Files A','side':LEFT},
               {'name':'SpawnInputDialogB','text':'Add Input Files B','side':LEFT},
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
              {'name':'skipohmic','text':'Skip bad dJ/dV','row':5,'tooltip':'Skip plots with d2J/dV2 < 0 between Vcutoff and Vmin/Vmax.'},
              {'name':'nomin','text':'Use dJ/dV for Vtrans','row':6,'tooltip':'Use dJ/dV plots to find the minimum of F-N plots when computing Vtrans.'},
              {'name':'logr','text':'Use log|R|','row':7,'tooltip':'Use log |R| instead of |R| when computing histograms.'},
              {'name':'lorenzian','text':'Lorenzian','row':8,'tooltip':'Fit a Lorenzian instead of a Gaussian.'},
              {'name':'autonobs','text':'Auto N','row':9,'tooltip':'Determine N (degrees of freedom) from metadata (if provided with the _data.txt files).'},
              {'name':'tracebyfile','text':'AFM Data','row':10,'tooltip':'Each file contains one (foward/backward) trace.'}
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
        

        Label(self.RightOptionsFrame, text="Output file base name:").grid(column=0,row=3)
        self.OutputFileName = Entry(self.RightOptionsFrame, width=20,
                font=Font(size=8,slant='italic'))
        for n in ('<Return>','<Leave>','<Enter>'):
            self.OutputFileName.bind(n, self.checkOutputFileName)
        self.OutputFileName.grid(column=0,row=4)
        
        if self.opts.outfile:
            self.OutputFileName.insert(0,self.opts.outfile)
        

        lbls = [
            {'name': 'Columns', 'text': "Columns to parse:",
              'tooltip': 'Columns from input data to parse as X/Y data.'},
            {'name': 'Vcutoff', 'text': "Cuttoff for d2J/dV2:",
             'tooltip': "Check the values of d2J/dV2 between |vcutoff| and Vmin/Vmax for line-shape filtering. Set to -1 for Vmin/Vmax."},
            {'name': 'Nobs', 'text': "N for p-test of J:",
                 'tooltip': "Compute p-values for J (not Gmean!) using this number of observations."},
            {'name': 'Threads', 'text': "Parser threads:",
             'tooltip': "Number of parser threads to start. Provides a small speed boost."}
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
        self.FileListBoxFrameLabel = Label(self.FileListBoxFrame,
                textvariable=self.FileListBoxFrameLabelVar,
                font=Font(size=10,weight="bold"))
        self.FileListBoxFrameLabel.grid(row=0, column=0, columnspan=4)
        self.UpdateFileListBoxFrameLabel()
        
        ab=('A','B')
        col=(0,2)
        for i in range(0,len(ab)):
            yScroll = Scrollbar(self.FileListBoxFrame, orient=VERTICAL)
            yScroll.grid(row=1, column=col[i]+1, sticky=N+S)
            xScroll = Scrollbar(self.FileListBoxFrame, orient=HORIZONTAL)
            xScroll.grid(row=2, column=col[i], sticky=E+W)
            setattr(self, 'filelist'+ab[i], StringVar())
            setattr(self, 'FileListBox'+ab[i], Listbox(self.FileListBoxFrame, listvariable=getattr(self,'filelist'+ab[i]), selectmode=EXTENDED, 
                            height = 20, width = 0, relief=RAISED, bd=1,
                                bg=WHITE,
                                                                #font=Font(family="Helvetica",size=10,slant='italic'),
                                xscrollcommand=xScroll.set, 
                                yscrollcommand=yScroll.set))
            getattr(self, 'FileListBox'+ab[i]).grid(row=1, column=col[i], sticky=N+S+E+W)
            xScroll['command'] = getattr(self,'FileListBox'+ab[i]).xview
            yScroll['command'] = getattr(self,'FileListBox'+ab[i]).yview
            self.FileListBoxFrame.columnconfigure(col[i],weight=1)
        self.FileListBoxFrame.rowconfigure(1,weight=1)

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
            if self.opts.plot:
                plotter = StatPlotter(gothread.statparser)
                self.logger.info("Generating plots...")
                try:
                        import matplotlib.pyplot as plt
                        plotter.DoPlots(plt)
                        plt.show(block=False)
                except ImportError as msg:
                        self.logger.error("Cannot import matplotlib! %s", str(msg), exc_info=False)
        else: 
            if len(self.opts.setA) and len(self.opts.setB):
                statparser = Stats(self.opts,handler=LoggingToGUI(self.Logging))
                self.gothreads.append(ParseThread(statparser))
                self.gothreads[-1].start()
                self.ButtonParse.after('500',self.Parse)
            else:
                self.logger.warn("No input files!")
    def RemoveFileClick(self):
        for ab in ('A','B'):
            flb = getattr(self,'FileListBox'+ab)
            selected = [ flb.get(x) for x in flb.curselection() ]
            todel = []
            filelist = []
            for i in range(0, len( getattr(self.opts,'set'+ab) )):
                for s in selected:
                    if getattr(self.opts, 'set'+ab)[i].replace(" ","_") == s:
                        todel.append(i)
            for i in range(0, len( getattr(self.opts,'set'+ab)  )):
                if i not in todel:
                    filelist.append(getattr(self.opts, 'set'+ab)[i])
            setattr(self.opts, 'set'+ab, filelist)
            self.updateFileListBox(ab)
            getattr(self,'FileListBox'+ab).selection_clear(0,END)

    def SpawnInputDialogClick(self,ab):
        self.checkOptions()
        ab=ab.upper()
        if ab == 'A': a,b='A','B'
        else: a,b = 'B','A'
        fs = getattr(self.opts, 'set'+ab)
        if getattr(self,'last_input_path'+a) == self.defaultdir:
            setattr(self,'last_input_path'+a, getattr(self,'last_input_path'+b))
        fs += filedialog.askopenfilename(title="Files to parse", multiple=True, \
                                 initialdir=getattr(self,'last_input_path'+ab),\
                                 filetypes=[('Data files','*_data.txt'),('Text files','*.txt'),('All files', '*')])
        setattr(self.opts, 'set'+ab, fs)
        if len( getattr(self.opts, 'set'+ab)):
            setattr(self,'last_input_path'+ab, os.path.split(getattr(self.opts, 'set'+ab)[0])[0])
            self.updateFileListBox(ab)
        if not self.opts.outfile:
            self.OutputFileName.delete(0,END)
            self.OutputFileName.insert(0,\
               os.path.basename(getattr(self.opts, 'set'+ab)[-1]).lower().replace('.txt',''))
            self.checkOutputFileName(None)

    def SpawnInputDialogAClick(self):
        self.SpawnInputDialogClick('A')
    def SpawnInputDialogBClick(self):
        self.SpawnInputDialogClick('B')

    def updateFileListBox(self,ab):
        ab = ab.upper()
        getattr(self, 'filelist'+ab).set(" ".join([x.replace(" ","_") for x in getattr(self.opts, 'set'+ab)]))
        self.Check_autonobs["state"] = NORMAL   
        for ab in ('A','B'):
            for f in getattr(self.opts, 'set'+ab):
                if not os.path.exists(f.replace('_data.txt','')):
                    self.autonobs.set(0)
                    self.Check_autonobs["state"] = DISABLED
        self.checkOptions()

    def SpawnOutputDialogClick(self):
        self.checkOptions()
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

        for n in ('Nobs', 'Threads'):
            try:
                var = int(getattr(self,'Entry'+n).get())
            except ValueError:
                var = -1
            if var > 0:
                setattr(self.opts,n.lower(),var)
            getattr(self,'Entry'+n).delete(0,END)
            getattr(self,'Entry'+n).insert(0,getattr(self.opts,n.lower()))
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



class LoggingToGUI(logging.Handler):
    """ Used to redirect logging output to the widget passed in parameters """
    def __init__(self, console):
        logging.Handler.__init__(self)
        self.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
        self.console = console #Any text widget, you can use the class above or not

    def emit(self, message): # Overwrites the default handler's emit method
        formattedMessage = self.format(message)  #You can change the format here
        self.console["state"] = NORMAL
        self.console.insert(END, formattedMessage+"\n") #Inserting the logger message in the widget
        self.console["state"] = DISABLED
        self.console.see(END)

