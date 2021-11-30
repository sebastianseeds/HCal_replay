from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar, Radiobutton
from tkinter.ttk import Style, Entry
import json
import os
from datetime import date, datetime

import Lib_HCal_Event_Display as fn
import ROOT
import sys
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)

read_all = 0
loop = 1
test_evt = 2#31823

data = "/home/skbarcus/Machine_Learning/Tutorials/Data/"

#infile = sys.argv[1]
#infile = "/home/skbarcus/JLab/SBS/HCal/Analysis/rootfiles/fadc_f1tdc_1710.root"
#infile = "/home/skbarcus/JLab/SBS/HCal/Analysis/Cosmics/rootfiles/fadc_f1tdc_820.root"
infile = "/data/daq/rootfiles/fadc_f1tdc_1727.root"

print("Reading from", infile)

inFile = ROOT.TFile.Open(infile," READ ")

tree = inFile.Get("T")

h1 = ROOT.TH1D("h1","h1 test",20,0,50000)

print(tree.GetEntries())

#for entryNum in range (0, loop):
    #tree.GetEntry(entryNum)

#Get the first entry.
tree.GetEntry(0)
#Get the number of fADC samples per single fADC event.
nsamps = getattr(tree,"sbs.hcal.nsamps")
nsamps = int(nsamps[0])
print('The number of fADC samples is '+str(nsamps)+'.')

#Get the number of PMTs.
adc_int = getattr(tree,"sbs.hcal.a")
npmt = int(len(adc_int))
print('The number of PMTs is '+str(npmt)+'.')

#Get the number of columns and rows (should generally be 12 cols and 24 rows for HCal).
nrow = getattr(tree,"sbs.hcal.row")
nrow = int(max(nrow))
ncol = getattr(tree,"sbs.hcal.col")
ncol = int(max(ncol))
print('There are '+str(ncol)+' columns and '+str(nrow)+' rows.')

#Get the info for the first event plots.
samps = getattr(tree,"sbs.hcal.samps")

if read_all == 1:
    loop = tree.GetEntries()

adc_samps_vals = []
adc_vals = []
tdc_vals = []
hit_vals = []

class DVCS_Pulser_Control_GUI:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("800x700")
        primary.title("DVCS Pulser Control GUI")

        #Create a frame to hold the exit button.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        #info_frame.pack(side='left', fill=BOTH, expand=True)
        info_frame.pack(side='left', fill=BOTH)

        #Create a label with the instructions to use the control GUI.
        instructions_label_text = StringVar()
        instructions_label_text.set("This GUI controls the DVCS pulser which\nis used to pulse the LEDs in HCal\nfor gain monitoring purposes. This pulser is\nused in conjunction with a NIM gate generator\nto produce logical pulses compatible with\nthe HCal LED power boxes.\n\n")
        instructions_label = Label(info_frame, width = 40,textvariable=instructions_label_text,font='Helvetica 12 bold')
        #instructions_label.pack(side="top")

        #Container canvas to hold display frame and a vertical scroll bar.
        container_canvas = Canvas(primary, borderwidth=0)
        container_canvas.pack(side='left', fill=BOTH, expand=True)
        vsb = Scrollbar(primary, orient="vertical", command=container_canvas.yview)
        container_canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")

        #Create frame to set the individual channel voltage outputs.
        display_frame = Frame(container_canvas, relief=RAISED, borderwidth=2)
        display_frame.pack(side='top', fill=BOTH, expand=True)

        container_canvas.create_window((0,0), window=display_frame, anchor="nw")
        display_frame.bind("<Configure>", lambda event, canvas=container_canvas: onFrameConfigure2(canvas))

        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            container_canvas.configure(scrollregion=container_canvas.bbox("all"))

        #v = Scrollbar(display_frame, orient='vertical')
        #v.config(command=display_frame.yview) #for vertical scrollbar
        #v.pack(side = RIGHT, fill = Y)

        #myscrollbar=Scrollbar(display_frame,orient="vertical")
        #myscrollbar.pack(side="right",fill="y")
        #myscrollbar.grid(row=0,column=13)

        #Create the fADC figures.
        #Arrays to hold matplotlib figures.
        fig_arr = []
        ax_arr = []
        hist_arr = []
        for pmt in range(0,npmt):
            xaxis = []
            xaxis.clear()
            for i in range(nsamps):
                xaxis.append(i)

            yaxis = []
            yaxis.clear()
            for i in range(nsamps):
                yaxis.append(samps[i])

            fig, ax = plt.subplots(figsize=(1,1))
            hist = ax.hist(xaxis,bins=nsamps,weights=yaxis,histtype='bar',ec='black')
            ax.set_ylim(100,500)
            ax.set_xlim(0,nsamps-1)
            fig_arr.append(fig)
            ax_arr.append(ax)
            hist_arr.append(hist)

        #Array to hold all the canvases for all the PMTs data to be drawn on.
        canvas_arr = []

        #Handler to pass the fig to the function as an argument.
        def handler(event,fig_arr=fig_arr,canvas_arr=canvas_arr):
            return fn.single_plot(event,fig_arr,canvas_arr)

        for pmt in range(0,npmt):
            canvas = FigureCanvasTkAgg(fig_arr[pmt], master=display_frame)
            canvas.get_tk_widget().bind("<Button-1>", handler)
            canvas_arr.append(canvas)
            #canvas.draw()

        #toolbar = NavigationToolbar2Tk(canvas,display_frame)
        #toolbar.update()
        #canvas.get_tk_widget().pack()

        #Create an entry for the event number to be displayed.
        display_entry = Entry(info_frame,width=8)
        display_btn = Button(info_frame, text='Display\nEvent', width=8, height=4, font='Helvetica 8 bold', command=lambda canvas_arr=canvas_arr, canvas=canvas,hist_arr=hist_arr,ax_arr=ax_arr, fig_arr=fig_arr, ax=ax, fig=fig, info_frame=info_frame, display_frame=display_frame, tree=tree, display_entry=display_entry: fn.display_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry), bg = "grey")
        #display_btn.pack(side="bottom")
        display_btn.grid(row=1,column=1)

        #display_entry.pack(side="bottom")
        display_entry.grid(row=0,column=1)
        display_entry.insert(0,0)

        #Create a button to step forward one event.
        next_btn = Button(info_frame, text='Next', width=8, height=4, font='Helvetica 8 bold', bg = "grey", command=lambda canvas_arr=canvas_arr, canvas=canvas,hist_arr=hist_arr,ax_arr=ax_arr, fig_arr=fig_arr, ax=ax, fig=fig, info_frame=info_frame, display_frame=display_frame, tree=tree, display_entry=display_entry: fn.next_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry))
        next_btn.grid(row=1,column=2)

        #Create a button to step back one event.
        previous_btn = Button(info_frame, text='Previous', width=8, height=4, font='Helvetica 8 bold', bg = "grey", command=lambda canvas_arr=canvas_arr, canvas=canvas,hist_arr=hist_arr,ax_arr=ax_arr, fig_arr=fig_arr, ax=ax, fig=fig, info_frame=info_frame, display_frame=display_frame, tree=tree, display_entry=display_entry: fn.previous_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry))
        previous_btn.grid(row=1,column=0)

        def key(event):
            print("pressed", repr(event.char))

        def callback(event):
            print("clicked at", event.x, event.y)

        #canvas_clickable= Canvas(display_frame, width=100, height=100, bd=2, relief=RAISED)
        #canvas_clickable.bind("<Key>", key)
        #canvas_clickable.bind("<Button-1>", callback)
        #canvas_clickable.pack()

root = Tk()
my_gui = DVCS_Pulser_Control_GUI(root)
root.mainloop()
