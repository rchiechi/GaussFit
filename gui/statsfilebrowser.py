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
import sys,os,platform
from tkinter import filedialog #some weird bug...
from tkinter import *
from tkinter.ttk import *
import logging
from gui.tooltip import *



class ChooseFiles(Frame):

	def __init__(self, opts, Go, master=None):
		Frame.__init__(self, master)
		self.Go = Go
		self.boolmap = {1:True, 0:False}
		try:
			self.defaultdir = os.environ['PWD']
		except KeyError:
			self.defaultdir = os.environ['HOME']
		self.last_input_pathA = self.defaultdir
		self.last_input_pathB = self.defaultdir
		self.opts = opts
		self.master.title("File Browser")
		self.master.geometry('1400x700+250-250')
		self.pack(fill=BOTH)
		self.createWidgets()
		self.ToFront()
		self.mainloop()

	def createWidgets(self):
			
		self.ButtonFrame = Frame(self)

		self.OptionsFrame = Frame(self)
		

		self.ColumnFrame = Frame(self)
		self.ColumnFrameLabel = Label(self.ColumnFrame, text="Columns to parse:" ).pack(side=TOP, ipadx=10)
		self.LeftOptionsFrame = Frame(self.ColumnFrame)
		
		self.FileListBoxFrame = Frame(self)
		self.FileListBoxFrameLabelVar = StringVar()
		self.FileListBoxFrameLabel = Label(self.FileListBoxFrame, textvariable=self.FileListBoxFrameLabelVar )
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
		self.LeftOptionsFrame.pack(side=LEFT)
		self.LoggingFrame.pack(side=BOTTOM)
		self.logger = logging.getLogger(None)
		self.logger.addHandler(LoggingToGUI(self.Logging))
		self.logger.info("Logging...")

		self.updateFileListBox('A')
		self.updateFileListBox('B')

	def createButtons(self):
			
		self.ButtonQuit = Button(self.ButtonFrame)
		self.ButtonQuit.config(text="QUIT", command=self.Quit)
		self.ButtonQuit.pack(side=BOTTOM)
		
		self.ButtonSpawnInputDialogA = Button(self.ButtonFrame)
		self.ButtonSpawnInputDialogA.config(text="Add Input Files A", command=self.SpawnInputDialogAClick)
		self.ButtonSpawnInputDialogA.pack(side=LEFT)
		
		self.ButtonSpawnInputDialogB = Button(self.ButtonFrame)
		self.ButtonSpawnInputDialogB.config(text="Add Input Files B", command=self.SpawnInputDialogBClick)
		self.ButtonSpawnInputDialogB.pack(side=LEFT)

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
		if len(self.opts.setA) and len(self.opts.setB):
			self.Go(self.opts)
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
			setattr(self,'last_input_path'+ab, os.path.split(self.opts.setA[0])[0])
			self.updateFileListBox(ab)

	def SpawnInputDialogAClick(self):
		self.SpawnInputDialogClick('A')
	def SpawnInputDialogBClick(self):
		self.SpawnInputDialogClick('B')

	def updateFileListBox(self,ab):
		ab = ab.upper()
		getattr(self, 'filelist'+ab).set(" ".join([x.replace(" ","_") for x in getattr(self.opts, 'set'+ab)]))

	def SpawnOutputDialogClick(self):
		outdir = filedialog.askdirectory(title="Select Output File(s)", initialdir=self.opts.out_dir)
		if os.path.exists(outdir):
			self.opts.out_dir = outdir
			self.UpdateFileListBoxFrameLabel()

	def UpdateFileListBoxFrameLabel(self):
		self.FileListBoxFrameLabelVar.set("Output to: %s/%s_*.txt"% (self.opts.out_dir, self.opts.outfile) )

		
	def createOptions(self):
##		Label( self.OptionsFrame,text="Select Options:").grid(column=1, row=0)

		self.plot = IntVar()
		self.write = IntVar()
		self.Check_plot = Checkbutton(self.OptionsFrame, text="Plot", \
							 variable=self.plot, command=self.checkOptions)
		self.Check_plot.grid(column=0,row=1,sticky=W)
		createToolTip(self.Check_plot, "Show summary plots after parsing.")
		
		self.Check_write = Checkbutton(self.OptionsFrame, text="Write", \
							  variable=self.write, command=self.checkOptions)
		self.Check_write.grid(column=0,row=2,sticky=W)
		createToolTip(self.Check_write, "Write results to text files after parsing.")

		Label(self.OptionsFrame, text="Output file base name:").grid(column=0,row=3)

		self.OutputFileName = Entry(self.OptionsFrame, width=16)
		self.OutputFileName.bind("<Return>", self.checkOutputFileName)
		self.OutputFileName.bind("<Leave>", self.checkOutputFileName)
		self.OutputFileName.bind("<Enter>", self.checkOutputFileName)
		self.OutputFileName.grid(column=0,row=4)

		self.skip = IntVar()
		self.Check_skip = Checkbutton(self.OptionsFrame, text="Skip bad dJ/dV", \
										 variable=self.skip, command=self.checkOptions)
		self.Check_skip.grid(column=0,row=5,sticky=W)
		createToolTip(self.Check_skip, "Skip plots with d2J/dV2 < 0 between Vcutoff and Vmin/Vmax.")

		self.nomin = IntVar()
		self.Check_smooth = Checkbutton(self.OptionsFrame, text="Use dJ/dV for Vtrans", \
										 variable=self.nomin, command=self.checkOptions)
		self.Check_smooth.grid(column=0,row=6,sticky=W)
		createToolTip(self.Check_smooth, "Use dJ/dV plots to find the minimum of F-N plots when computing Vtrans.")

		self.logr = IntVar()
		self.logr.set(True)
		self.Check_logr = Checkbutton(self.OptionsFrame, text="Use log|R|", \
										 variable=self.logr, command=self.checkOptions)
		self.Check_logr.grid(column=0,row=7,sticky=W)
		createToolTip(self.Check_logr, "Use log|R| when computing histograms.")

		self.lorenzian = IntVar()
		self.Check_lorenzian = Checkbutton(self.OptionsFrame, text="Lorenzian", \
										 variable=self.lorenzian, command=self.checkOptions)
		self.Check_lorenzian.grid(column=0,row=8,sticky=W)
		createToolTip(self.Check_lorenzian, "Fit a Lorenzian instead of a Gaussian.")


		Label(self.LeftOptionsFrame, text="Cuttoff for d2J/dV2:").grid(column=0,row=0)
		self.EntryVcutoff = Entry(self.LeftOptionsFrame, width=8)
		self.EntryVcutoff.bind("<Return>", self.checkOptions)
		self.EntryVcutoff.bind("<Leave>", self.checkOptions)
		self.EntryVcutoff.bind("<Enter>", self.checkOptions)
		self.EntryVcutoff.grid(column=0,row=1)
		createToolTip(self.EntryVcutoff, "Check the values of d2J/dV2 between |vcutoff| and Vmin/Vmax for line-shape filtering. Set to -1 for Vmin/Vmax.")

		#Label(self.LeftOptionsFrame, text="Smoothing parameter:").grid(column=0,row=2)
		#self.Smoothingcutoff = Entry(self.LeftOptionsFrame, width=8)
		#self.Smoothingcutoff.bind("<Return>", self.checkOptions)
		#self.Smoothingcutoff.bind("<Leave>", self.checkOptions)
		#self.Smoothingcutoff.bind("<Enter>", self.checkOptions)
		#self.Smoothingcutoff.grid(column=0,row=3)
		#createToolTip(self.Smoothingcutoff, "The cutoff value for the residual squares (the difference between experimental data points and the fit). The default is 1e-12. Set to 0 to disable smoothing.")
		
		#Label(self.LeftOptionsFrame, text="Bins for J/R Histograms:").grid(column=0,row=6)
		#self.EntryJRBins = Entry(self.LeftOptionsFrame, width=8)
		#self.EntryJRBins.bind("<Return>", self.checkOptions)
		#self.EntryJRBins.bind("<Leave>", self.checkOptions)
		#self.EntryJRBins.bind("<Enter>", self.checkOptions)
		#self.EntryJRBins.grid(column=0,row=7)
		#createToolTip(self.EntryJRBins, "Set binning for histograms of J and R.")
		
		#Label(self.LeftOptionsFrame, text="Bins for G Histograms:").grid(column=0,row=8)
		#self.Entryhmbins = Entry(self.LeftOptionsFrame, width=8)
		#self.Entryhmbins.bind("<Return>", self.checkOptions)
		#self.Entryhmbins.bind("<Leave>", self.checkOptions)
		#self.Entryhmbins.bind("<Enter>", self.checkOptions)
		#self.Entryhmbins.grid(column=0,row=9)
		#createToolTip(self.Entryhmbins, "Set binning for conductance heatmap histograms.")

		# checkGminmaxEntry call must come last
		#Label(self.LeftOptionsFrame, text="Y-scale for conductance:").grid(column=0,row=4)
		#self.EntryGminmax = Entry(self.LeftOptionsFrame, width=8)
		#self.EntryGminmax.bind("<Return>", self.checkGminmaxEntry)
		#self.EntryGminmax.bind("<Leave>", self.checkGminmaxEntry)
		#self.EntryGminmax.bind("<Enter>", self.checkGminmaxEntry)
		#self.EntryGminmax.grid(column=0,row=5)
		#createToolTip(self.EntryGminmax, "Set Ymin,Ymax for the conductance plot (lower-left of plot output).")
		#self.checkGminmaxEntry(None)



		if self.opts.write:
			self.write.set(1)
		if self.opts.plot:
			self.plot.set(1)
		if self.opts.outfile:
			self.OutputFileName.insert(0,self.opts.outfile)
		if self.opts.skipohmic:
			self.skip.set(1)
		if self.opts.nomin:
			self.nomin.set(1)

		self.checkOptions()

	def checkOptions(self, event=None):
		self.opts.plot = self.boolmap[self.plot.get()]
		self.opts.write = self.boolmap[self.write.get()]
		self.opts.skipohmic = self.boolmap[self.skip.get()]
		self.opts.nomin = self.boolmap[self.nomin.get()]
		self.opts.logr = self.boolmap[self.logr.get()]
		self.opts.lorenzian = self.boolmap[self.lorenzian.get()]
		
		if not self.opts.write:
			self.opts.plot = True
			self.plot.set(1)
			self.Check_plot["state"]=DISABLED
		else:
			self.Check_plot["state"]=NORMAL

		#if not self.opts.plot:
		#	self.EntryGminmax["state"]=DISABLED
		#else:
		#	self.EntryGminmax["state"]=NORMAL

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

		#try:
		#	self.opts.smooth = abs(float(self.Smoothingcutoff.get()))	 
		#except ValueError:
		#	self.opts.smooth = 1e-12

		#self.Smoothingcutoff.delete(0, END)
		#self.Smoothingcutoff.insert(0,self.opts.smooth)

		#try:
		#	self.opts.bins = abs(int(self.EntryJRBins.get()))
		#except ValueError:
		#	self.opts.bins = 50

		#self.EntryJRBins.delete(0, END)
		#self.EntryJRBins.insert(0,self.opts.bins)

		#try:
		#	self.opts.heatmapbins = abs(int(self.Entryhmbins.get()))
		#except ValueError:
		#	self.opts.heatmapbins = 25

		#self.Entryhmbins.delete(0, END)
		#self.Entryhmbins.insert(0,self.opts.heatmapbins)
	
	


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

#	def checkGminmaxEntry(self, event):
#		self.checkOptions()
#		try:
#			x, y = self.EntryGminmax.get().split(",")
#			self.opts.mlow, self.opts.mhi = int(x), int(y)
#		except ValueError as msg:
#			pass
#		self.EntryGminmax.delete(0, END)
#		self.EntryGminmax.insert(0, ",".join( (str(self.opts.mlow), str(self.opts.mhi)) ))

	def createOutputLabel(self):
		ab=('A','B')
		for i in range(0,2):
			setattr(self,'yScroll'+ab[i], Scrollbar(self.FileListBoxFrame, orient=VERTICAL))
			getattr(self,'yScroll'+ab[i]).grid(row=1, column=1, sticky=N+S)

			setattr(self,'xScroll'+ab[i], Scrollbar(self.FileListBoxFrame, orient=HORIZONTAL))
			getattr(self,'xScroll'+ab[i]).grid(row=2, column=0, sticky=E+W)
	 
			setattr(self, 'filelist'+ab[i], StringVar())
			setattr(self, 'FileListBox'+ab[i], Listbox(self.FileListBoxFrame, listvariable=getattr(self,'filelist'+ab[i]), selectmode=EXTENDED, 
							height = 20, width = 75, relief=RAISED, bd=1, 
										  xscrollcommand=getattr(self,'xScroll'+ab[i]).set, 
										  yscrollcommand=getattr(self,'yScroll'+ab[i]).set))
			getattr(self, 'FileListBox'+ab[i]).grid(row=1, column=i, sticky=N+S+E+W)
			getattr(self,'xScroll'+ab[i])['command'] = getattr(self,'FileListBox'+ab[i]).xview
			getattr(self,'yScroll'+ab[i])['command'] = getattr(self,'FileListBox'+ab[i]).yview






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


#if __name__ == "__main__":
#
#	class Opts():
#			def __init__(self):
#					self.setA = []
#					self.setB = []
#					self.Xcol = 0
#					self.Ycol = 2
#					self.plot = True
#					self.write= True
#					self.out_dir = os.environ['PWD']
#					self.outfile = "test"
#					logging.basicConfig(level=logging.INFO,format = os.path.basename(sys.argv[0])+' %(levelname)s %(message)s')
#	def Go(arg):
#		logging.info("Dummy function")
#		print("Not to logger.")
#		return
# 
#
#	gui = ChooseFiles(Opts(),Go)
 