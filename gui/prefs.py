'''A separate preferences window for settings Opts'''

# import os
# import platform
#import tkinter as tk
import tkinter.ttk as tk
from tkinter import Toplevel,Y,BOTTOM,LEFT,RIGHT,END,BOTH
from gui.colors import GREY
from gui.tooltip import createToolTip

BOOLMAP = {1:True, 0:False}

class PreferencesWindow(Toplevel):
    '''A new window to hold the preferences pane'''
    def __init__(self, parent, opts):
        super().__init__(parent)
        self.title("Preferences")
        self.prefs = PreferencesFrame(self, opts)


class PreferencesFrame(tk.Frame):
    '''A preferences frame for settings options
       to be called from the main FileBrowser
    '''
    def __init__(self, master, opts):
        tk.Frame.__init__(self, master)
        self.opts = opts
        self.master.tk_setPalette(background=GREY,
            activeBackground=GREY)
        self.master.geometry('800x800+250-250')
        self.pack(fill=BOTH)
        self.__createFrames()
        self.__createButtons()
        self.__createOptions()
        self.mainloop()

    
    
    def __createFrames(self):

        self.ButtonFrame = tk.Frame(self)
        self.LeftOptionsFrame = tk.Frame(self)
        self.RightOptionsFrame = tk.Frame(self)
        self.ButtonFrame.pack(side=BOTTOM,fill=None)
        self.LeftOptionsFrame.pack(side=LEFT,fill=Y)
        self.RightOptionsFrame.pack(side=RIGHT,fill=Y)



    def __createButtons(self):

        buttons = [
               {'name':'Save','text':'Save','side':BOTTOM},
            #    {'name':'Quit','text':'QUIT','command':'Quit','side':BOTTOM},
            #    {'name':'SpawnInputDialog','text':'Add Input Files','side':LEFT},
            #    {'name':'RemoveFile','text':'Remove Files','side':LEFT},
            #    {'name':'SpawnOutputDialog','text':'Choose Output Directory','side':LEFT},
            #    {'name':'Parse','text':'Parse!','side':LEFT}
               ]

        for _b in buttons:
            button = tk.Button(self.ButtonFrame)
            button.config(text=_b['text'],command=getattr(self,_b['name']+'Click'))
            button.pack(side=_b['side'])
            setattr(self,'Button'+_b['name'],button)

    def __createOptions(self):
    
        lbls = [
            {'name': 'Columns', 'text': "Columns to parse:",
             'tooltip': 'Columns from input data to parse as X/Y data.',
             'frame':self.LeftOptionsFrame},
            {'name': 'Vcutoff', 'text': "Cuttoff for d2J/dV2:",
             'tooltip': "Check the values of d2J/dV2 between |vcutoff| and Vmin/Vmax for line-shape filtering. Set to -1 for Vmin/Vmax.",
             'frame':self.LeftOptionsFrame},
            {'name': 'Lagcutoff', 'text': "Cuttoff for lag plot filter:",
             'tooltip': "Throw out J-values whose euclidian distance from a linear fit of the lag plot exceed this value for computing filtered histograms.",
             'frame':self.LeftOptionsFrame},
            {'name': 'Smooth', 'text': "Smoothing parameter:",
             'tooltip': "The cutoff value for the residual squares (the difference between experimental data points and the fit). The default is 1e-12. Set to 0 to disable smoothing.",
             'frame':self.LeftOptionsFrame},
            {'name': 'Gminmax', 'text': "Y-scale for G heatmap:",
                 'tooltip': "Set Ymin,Ymax for the heapmap plot (lower-left of plot output).",
             'frame':self.LeftOptionsFrame},
            {'name': 'NDCminmax', 'text': "Y-scale for NDC heatmap:",
                 'tooltip': "Set Ymin,Ymax for the heapmap plot (lower-left of plot output).",
             'frame':self.LeftOptionsFrame},
            {'name': 'Bins', 'text': "Bins for J/R Histograms:",
             'tooltip': "Set binning for histograms of J and R.",
             'frame':self.LeftOptionsFrame},
            {'name': 'Heatmapbins', 'text': "Bins for G Histograms:",
             'tooltip': "Set binning for heatmap histograms.",
             'frame':self.LeftOptionsFrame},
            {'name': 'Alpha', 'text': "⍺-value for CI:",
              'tooltip': "The p-cutoff is 1-⍺, e.g., ⍺ = 0.05 => 95% cutoff.",
             'frame':self.LeftOptionsFrame},
            {'name': 'MaxR', 'text': "Maximum log|R|/R value:",
              'tooltip': "Maximum allowable value of log|R| or (R if -R). Values above MaxR will excluded COMPLETELY from processing!",
             'frame':self.RightOptionsFrame},
            ]

        i = 0
        for _l in lbls:
            tk.Label(_l['frame'], text=_l['text']).grid(column=0,row=i)
            entry = tk.Entry(_l['frame'], width=8,
                    validate='focus', validatecommand=self.checkOptions)
            entry.bind("<Return>", self.checkOptions)
            entry.bind("<Leave>", self.checkOptions)
            entry.bind("<Enter>", self.checkOptions)
            entry.bind("<Tab>", self.checkOptions)
            entry.grid(column=0,row=i+1)
            createToolTip(entry, _l['tooltip'])
            setattr(self, 'Entry'+_l['name'], entry)
            i+=2

        self.checkOptions()

    def SaveClick(self):
        self.checkOptions()
        self.master.destroy()

    def checkOptions(self, event=None): #pylint: disable=W0613
        for n in (('Smooth',1e-12),('Bins',50),('Heatmapbins',25)):
            try:
                var = getattr(self,'Entry'+n[0]).get()
                if 'Smooth' in n:
                    var = float(var)
                else:
                    var = int(var)
            except ValueError:
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

        self.checkGminmaxEntry()
        self.checkNDCminmaxEntry()
        self.checkColumnEntry()
        self.checkVcutoffEntry()
        self.checkLagcutoffEntry()
        self.checkAlphaEntry()
        self.checkMaxREntry()

    def checkVcutoffEntry(self):
        try:
            vcutoff = float(self.EntryVcutoff.get()) #pylint: disable=E1101
            if vcutoff != -1:
                vcutoff = abs(vcutoff)
            self.opts.vcutoff = vcutoff
        except ValueError:
            self.opts.vcutoff = -1

        self.EntryVcutoff.delete(0, END) #pylint: disable=E1101
        if self.opts.vcutoff > 0:
            self.EntryVcutoff.insert(0,self.opts.vcutoff) #pylint: disable=E1101
        else:
            self.EntryVcutoff.insert(0,"Vmax") #pylint: disable=E1101

    def checkLagcutoffEntry(self):
        try:
            lagcutoff = float(self.EntryLagcutoff.get()) #pylint: disable=E1101
            self.opts.lagcutoff = abs(lagcutoff)
        except ValueError:
            self.opts.lagcutoff = 0.1

        self.EntryLagcutoff.delete(0, END) #pylint: disable=E1101
        self.EntryLagcutoff.insert(0,self.opts.lagcutoff) #pylint: disable=E1101

    def checkColumnEntry(self):
        try:
            x, y = self.EntryColumns.get().split(",") #pylint: disable=E1101
            self.opts.xcol, self.opts.ycol = int(x)-1, int(y)-1
        except ValueError:
            pass
        self.EntryColumns.delete(0, END) #pylint: disable=E1101
        self.EntryColumns.insert(0, ",".join(( str(self.opts.xcol+1), str(self.opts.ycol+1) ))) #pylint: disable=E1101

    def checkGminmaxEntry(self):
        try:
            x, y = self.EntryGminmax.get().split(",") #pylint: disable=E1101
            self.opts.mlow, self.opts.mhi = int(x), int(y)
        except ValueError:
            pass
        self.EntryGminmax.delete(0, END) #pylint: disable=E1101
        self.EntryGminmax.insert(0, ",".join( (str(self.opts.mlow), str(self.opts.mhi)) )) #pylint: disable=E1101

    def checkNDCminmaxEntry(self):
        try:
            x, y = self.EntryNDCminmax.get().split(",") #pylint: disable=E1101
            self.opts.ndc_mlow, self.opts.ndc_mhi = float(x), float(y)
        except ValueError:
            pass
        self.EntryNDCminmax.delete(0, END) #pylint: disable=E1101
        self.EntryNDCminmax.insert(0, ",".join( (str(self.opts.ndc_mlow), str(self.opts.ndc_mhi)) )) #pylint: disable=E1101

    def checkAlphaEntry(self):
        try:
            alpha = float(self.EntryAlpha.get()) #pylint: disable=E1101
            if 0 < alpha < 1:
                self.opts.alpha = alpha
        except ValueError:
            pass
        self.EntryAlpha.delete(0, END) #pylint: disable=E1101
        self.EntryAlpha.insert(0,self.opts.alpha) #pylint: disable=E1101

    def checkMaxREntry(self):
        try:
            self.opts.maxr = float(self.EntryMaxR.get()) #pylint: disable=E1101
        except ValueError:
            pass
        self.EntryMaxR.delete(0, END) #pylint: disable=E1101
        self.EntryMaxR.insert(0,self.opts.maxr) #pylint: disable=E1101
