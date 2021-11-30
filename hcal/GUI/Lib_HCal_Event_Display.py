#Library of functions for HCal_Event_Display.py.
from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar, Radiobutton
from tkinter.ttk import Style, Entry
import json
import os
from datetime import date, datetime

import ROOT
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)

def test():
    print('Test')

def test_command():
    print('Test.')

#fig, ax = plt.subplots(figsize=(4,4))
#canvas = FigureCanvasTkAgg(fig, master=display_frame)

def display_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry):
    #Get the event to display from the entry box.
    event = display_entry.get()
    print('Now displaying event '+str(event)+'.')

    #Find the correct event in the ROOT tree.
    event = int(event)
    tree.GetEntry(event)

    #Check how many fADC bins are in the window.
    nsamps = getattr(tree,"sbs.hcal.nsamps")
    nsamps = int(nsamps[0])
    print('The number of fADC samples is '+str(nsamps)+'.')

    #Get the fADC sample values vector. (Size per event is number of pmts * number of fADC bins in the window per PMT)
    samps = getattr(tree,"sbs.hcal.samps")

    #Get the TDC values vector. (Size per event is number of PMTs)
    tdc = getattr(tree,"sbs.hcal.tdc")

    #Get the number of PMTs.
    adc_int = getattr(tree,"sbs.hcal.a")
    npmt = int(len(adc_int))
    print('The number of PMTs is '+str(npmt)+'.')

    #Get the number of columns and rows (should generally be 12 cols and 24 rows for HCal).
    nrow = getattr(tree,"sbs.hcal.row")
    nrow = int(max(nrow))
    ncol = getattr(tree,"sbs.hcal.col")
    ncol = int(max(ncol))

    xaxis = []
    yaxis = []
    for pmt in range(0,npmt):
        fig_arr[pmt].set_figheight(1)
        fig_arr[pmt].set_figwidth(1)

        xaxis.clear()
        for i in range(nsamps):
            xaxis.append(i)

        yaxis.clear()
        first_bin = pmt*nsamps
        last_bin = nsamps+pmt*nsamps
        for i in range(first_bin,last_bin):
            #print('PMT '+str(pmt)+': first_bin='+str(first_bin)+' last_bin='+str(last_bin)+' samps[i]='+str(samps[i])+'.')
            yaxis.append(samps[i])

        ax_arr[pmt].clear()
        #ax_arr[pmt].set_axis_off()
        #ax_arr[pmt].set_xticks([])
        #ax_arr[pmt].set_yticks([])
        #hist1 = ax_arr[pmt].hist(xaxis,bins=nsamps,weights=yaxis,histtype='bar',ec='black')

        #Change the color of PMTs that have a TDC hit (could use an ADC threshold instead).
        if tdc[pmt]!=0:
            hist_arr[pmt] = ax_arr[pmt].hist(xaxis,bins=nsamps,weights=yaxis,histtype='bar',color='orangered',ec='darkorange',lw=1)
        else:
            hist_arr[pmt] = ax_arr[pmt].hist(xaxis,bins=nsamps,weights=yaxis,histtype='bar',ec='c',lw=1)
        ax_arr[pmt].set_ylim(100,500)
        ax_arr[pmt].set_xlim(0,nsamps-1)

        row = int(pmt/ncol)
        col = pmt%ncol
        canvas_arr[pmt].draw()
        canvas_arr[pmt].get_tk_widget().grid(row=row,column=col)
        #print('Plotting PMT '+str(pmt)+' at grid location (row,col) = ('+str(row)+','+str(col)+').')

    #print(display_frame.winfo_children())

    #Plot all of the PMT fADCs on the canvases.
    #canvas.draw()
    #canvas.get_tk_widget().grid(row=0,column=0)
    #for row in range(0,nrow):
        #for col in range(0,ncol):
            #canvas_arr[row*ncol+col].draw()
            #canvas_arr[row*ncol+col].get_tk_widget().grid(row=row,column=col)

def next_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry):
    #Get the current value in the entry box.
    event = display_entry.get()
    #print('Event was '+str(event))

    #Clear the entry box.
    entry_len = len(event)
    display_entry.delete(0,entry_len)

    #Insert the new event number in the entry box.
    display_entry.insert(0,str(int(event)+1))
    event = display_entry.get()
    #print('Event is currently '+str(event))
    display_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry)

def previous_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry):
    #Get the current value in the entry box.
    event = display_entry.get()
    print('Event was '+str(event))

    #Clear the entry box.
    entry_len = len(event)
    display_entry.delete(0,entry_len)

    #Insert the new event number in the entry box.
    display_entry.insert(0,str(int(event)-1))
    event = display_entry.get()
    print('Event is currently '+str(event))
    display_event(canvas_arr,canvas,hist_arr,ax_arr,fig_arr,ax,fig,info_frame,display_frame,tree,display_entry)

def single_plot(event,fig_arr,canvas_arr):
    single_plot_window = Tk()
    single_plot_window.wm_title("Single fADC Plot")
    #single_plot_label = Label(single_plot_window, text="Input")
    #single_plot_label.grid(row=0, column=0)
    #single_plot_btn = Button(single_plot_window, text="Close", width=8, height=4, command=single_plot_window.destroy)
    #single_plot_btn.grid(row=0,column=0)

    #Determine which canvas was clicked to find the PMT number to plot the correct fADC waveform.
    canvas = event.widget
    canvas_name = str(canvas)
    #print('canvas_name = '+canvas_name)
    #Get 3rd to last character in canvas_name to work out what the PMT number is. The canvases go canvas, canvas2,...,canvas288.
    if canvas_name[-4:-3]=='n':
        pmt=0
    elif canvas_name[-4:-3]=='v':
        pmt = int(canvas_name[-1:])-1
    elif canvas_name[-4:-3]=='a':
        pmt = int(canvas_name[-2:])-1
    elif canvas_name[-4:-3]=='s':
        pmt = int(canvas_name[-3:])-1
    else:
        pmt=0
        print('ERROR: Could not identify PMT number of the fADC waveform canvas clicked. Defaulting to using PMT 0.')
    #if canvas in canvas_arr:
        #print('Found it!')
    #print(canvas)
    #fig.figure(figsize=(8,8))
    #pmt=0
    fig_arr[pmt].set_figheight(6)
    fig_arr[pmt].set_figwidth(6)

    canvas = FigureCanvasTkAgg(fig_arr[pmt], master=single_plot_window)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0,column=0)
