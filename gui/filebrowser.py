#!/usr/bin/python3
import sys,os,platform
from tkinter import filedialog #some weird bug...
from tkinter import *
from tkinter.ttk import *
import logging



class ChooseFiles(Frame):

    def __init__(self, opts, Go, master=None):
        Frame.__init__(self, master)
        self.Go = Go
        self.boolmap = {1:True, 0:False}
        self.last_input_path = os.environ['HOME']
        self.opts = opts
        self.master.title("File Browser")
        self.master.geometry('850x600+250-250')
        self.pack(fill=BOTH)
        self.createWidgets()
        self.ToFront()
        self.mainloop()

    def createWidgets(self):
            
        self.ButtonFrame = Frame(self)

        self.OptionsFrame = Frame(self)

        self.ColumnFrame = Frame(self)
        self.ColumnFrameLabel = Label(self.ColumnFrame, text="Columns to parse:" ).pack(side=TOP, ipadx=10)

        self.FileListBoxFrame = Frame(self)
        self.FileListBoxFrameLabelVar = StringVar()
        self.FileListBoxFrameLabel = Label(self.FileListBoxFrame, \
                                              textvariable=self.FileListBoxFrameLabelVar )
        self.FileListBoxFrameLabel.grid(row=0, column=0)
        self.UpdateFileListBoxFrameLabel()

        self.LoggingFrame = Frame(self)
        self.Logging = Text(self.LoggingFrame, height=20)
        self.Logging.pack()

        self.createButtons()
        self.createOutputLabel()
        self.createOptions()
        self.createColumnEntry()

        self.ButtonFrame.pack(side=BOTTOM, fill=Y)
        self.FileListBoxFrame.pack(side=TOP, fil=Y, expand=1)
        self.OptionsFrame.pack(side=RIGHT)
        self.ColumnFrame.pack(side=LEFT)
        self.LoggingFrame.pack(side=BOTTOM)
        self.logger = logging.getLogger(None)
        self.logger.addHandler(LoggingToGUI(self.Logging))
        self.logger.info("Logging...")


        self.updateFileListBox()

    def createButtons(self):
            
        self.ButtonQuit = Button(self.ButtonFrame)
        self.ButtonQuit.config(text="QUIT", command=self.Quit)
        self.ButtonQuit.pack(side=BOTTOM)
        
        self.ButtonSpawnInputDialog = Button(self.ButtonFrame)
        self.ButtonSpawnInputDialog.config(text="Add Input Files", command=self.SpawnInputDialogClick)
        self.ButtonSpawnInputDialog.pack(side=LEFT)

        self.ButtonRemoveFile = Button(self.ButtonFrame)
        self.ButtonRemoveFile.config(text="Remove Input Files",command=self.RemoveFileClick)
        self.ButtonRemoveFile.pack(side=LEFT)
        
        self.ButtonSpawnOutputDialog = Button(self.ButtonFrame)
        self.ButtonSpawnOutputDialog.config(text="Choose Output Directory",command=self.SpawnOutputDialogClick)
        self.ButtonSpawnOutputDialog.pack(side=LEFT)
        
        self.ButtonGo = Button(self.ButtonFrame)
        self.ButtonGo.config(text="Parse!", command=self.Parse)
        self.ButtonGo.pack(side=LEFT)

    def ToFront(self):
        if platform.system() == "Darwin":
            os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        else:
            self.master.attributes('-topmost', 1)
            self.master.attributes('-topmost', 0)
        self.master.lift()

    def Quit(self):
        self.master.destroy()

    def Parse(self):
        if len(self.opts.in_files):
            self.Go(self.opts)
        else:
            self.logger.warn("No input files!")

    def RemoveFileClick(self):
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
        self.opts.in_files += filedialog.askopenfilename(title="Files to parse", multiple=True, \
                                         initialdir=self.last_input_path,\
                                         filetypes=[('Data files','*_data.txt'),('Text files','*.txt'),('All files', '*')])
        if len(self.opts.in_files):
            self.last_input_path = os.path.split(self.opts.in_files[0])[0]
            self.updateFileListBox()


    def updateFileListBox(self):
        self.filelist.set(" ".join([x.replace(" ","_") for x in self.opts.in_files]))
            
    def SpawnOutputDialogClick(self):
        outdir = filedialog.askdirectory(title="Select Output File(s)", initialdir=self.opts.out_dir)
        if os.path.exists(outdir):
            self.opts.out_dir = outdir
            self.UpdateFileListBoxFrameLabel()

    def UpdateFileListBoxFrameLabel(self):
        self.FileListBoxFrameLabelVar.set("Output to: %s/%s_*.txt"% (self.opts.out_dir, self.opts.outfile) )

        
    def createOptions(self):
##        Label( self.OptionsFrame,text="Select Options:").grid(column=1, row=0)

        self.plot = IntVar()
        self.write = IntVar()
        self.Check_plot = Checkbutton(self.OptionsFrame, text="Plot", \
                                         variable=self.plot, command=self.checkOptions)
        self.Check_plot.grid(column=0,row=1,sticky=W)
        
        self.Check_write = Checkbutton(self.OptionsFrame, text="Write", \
                                          variable=self.write, command=self.checkOptions)
        self.Check_write.grid(column=0,row=2,sticky=W)

        Label(self.OptionsFrame, text="Output file base name:").grid(column=0,row=3)

        self.OutputFileName = Entry(self.OptionsFrame, width=16)
        self.OutputFileName.bind("<Return>", self.checkOutputFileName)
        self.OutputFileName.bind("<Leave>", self.checkOutputFileName)
        self.OutputFileName.bind("<Enter>", self.checkOutputFileName)
        self.OutputFileName.grid(column=0,row=4)

        if self.opts.write:
            self.write.set(1)
        if self.opts.plot:
            self.plot.set(1)
        if self.opts.outfile:
            self.OutputFileName.insert(0,self.opts.outfile)

    def checkOptions(self):
        self.opts.plot = self.boolmap[self.plot.get()]
        self.opts.write = self.boolmap[self.write.get()]
           
        if not self.opts.write:
            self.opts.plot = True
            self.plot.set(1)
            self.Check_plot["state"]=DISABLED
        else:
            self.Check_plot["state"]=NORMAL

    
    def createColumnEntry(self):

        self.EntryColumns = Entry(self.ColumnFrame, width=8)
        self.EntryColumns.bind("<Return>", self.checkColumnEntry)
        self.EntryColumns.bind("<Leave>", self.checkColumnEntry)
        self.EntryColumns.bind("<Enter>", self.checkColumnEntry)
        self.EntryColumns.pack()
        self.checkColumnEntry(None)
        
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
        

    def createOutputLabel(self):

        self.yScroll = Scrollbar(self.FileListBoxFrame, orient=VERTICAL)
        self.yScroll.grid(row=1, column=1, sticky=N+S)

        self.xScroll = Scrollbar(self.FileListBoxFrame, orient=HORIZONTAL)
        self.xScroll.grid(row=2, column=0, sticky=E+W)

 
        self.filelist = StringVar()
        self.FileListBox = Listbox(self.FileListBoxFrame, listvariable=self.filelist, selectmode=EXTENDED, 
                        height = 20, width = 100, relief=RAISED, bd=1, 
                                      xscrollcommand=self.xScroll.set, 
                                      yscrollcommand=self.yScroll.set)
        self.FileListBox.grid(row=1, column=0, sticky=N+S+E+W)
        self.xScroll['command'] = self.FileListBox.xview
        self.yScroll['command'] = self.FileListBox.yview




#############


class LoggingToGUI(logging.Handler):
    """ Used to redirect logging output to the widget passed in parameters """
    def __init__(self, console):
        logging.Handler.__init__(self)
        self.console = console #Any text widget, you can use the class above or not

    def emit(self, message): # Overwrites the default handler's emit method
        formattedMessage = self.format(message)  #You can change the format here
        self.console["state"] = NORMAL
        self.console.insert(END, formattedMessage+"\n") #Inserting the logger message in the widget
        self.console["state"] = DISABLED
        self.console.see(END)


if __name__ == "__main__":

    class Opts():
            def __init__(self):
                    self.in_files = []
                    self.Xcol = 0
                    self.Ycol = 2
                    self.plot = True
                    self.write= True
                    self.out_dir = os.environ['PWD']
                    self.outfile = "test"
                    logging.basicConfig(level=logging.INFO,format = os.path.basename(sys.argv[0])+' %(levelname)s %(message)s')
    def Go(arg):
        logging.info("Dummy function")
        print("Not to logger.")
        return
 

    gui = ChooseFiles(Opts(),Go)
 
