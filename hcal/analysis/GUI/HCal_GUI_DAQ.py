from tkinter import Tk, Label, Button, StringVar, W, E, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
#import tkinter.ttk as ttk
from tkinter.ttk import Style, Entry, Frame
import json
#import tkinter.font as tkFont

nrows = 24
ncols = 12
channels = nrows * ncols
pmt_mods = []    #Holds PMT module numbers.
buttons = []

#Dictionary containing all of the relevant connections for each PMT.
# reading the data from the file 
with open('hcal_connections.json') as f: 
    connections = f.read()
connections = json.loads(connections)
#print(connections)

class MyFirstGUI:
    def __init__(self, primary):
        self.primary = primary
        primary.title("Front-End Electronics")

        self.style = Style()
        self.style.theme_use("default")

        #Create frames to hold the contents of each DAQ rack.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')#, fill=BOTH, expand=True)
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")#, padx=5, pady=5, expand=True, fill=BOTH)

        RR4_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR4_frame.pack(side='left', fill=BOTH, expand=True)

        RR5_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR5_frame.pack(side='left', fill=BOTH, expand=True)

        #Create buttons and labels for electronics in RR4.
        RR4_label = Label(RR4_frame, text='RR4', font='Helvetica 16 bold')
        RR4_label.pack(side="top")

        RR4_tdc_pp_btn = Button(RR4_frame, text='TDC Patch Panels', width=10, height=6, font='Helvetica 14 bold')
        RR4_tdc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR4_tdc_pp_btn.bind("<Button-1>", self.RR4_tdc_pp_connections)

        RR4_tdc_disc_btn = Button(RR4_frame, text='TDC Discriminators', width=10, height=3, font='Helvetica 14 bold')
        RR4_tdc_disc_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR4_tdc_disc_btn.bind("<Button-1>", self.RR4_tdc_disc_connections)

        RR4_adc_pp_btn = Button(RR4_frame, text='ADC Patch Panels', width=10, height=6, font='Helvetica 14 bold')
        RR4_adc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR4_adc_pp_btn.bind("<Button-1>", self.RR4_adc_pp_connections)

        #Create buttons and labels for electronics in RR5.
        RR5_label = Label(RR5_frame, text='RR5', font='Helvetica 16 bold')
        RR5_label.pack(side="top")

        RR5_nim_btn = Button(RR5_frame, text='NIM Electronics', width=10, height=3, font='Helvetica 14 bold')
        RR5_nim_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        #RR5_nim_btn.bind("<Button-1>", self.RR5_tdc_pp_connections)

        RR5_vxs2_btn = Button(RR5_frame, text='VXS 2 (intelsbshcal2)', width=10, height=6, font='Helvetica 14 bold')
        RR5_vxs2_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR5_vxs2_btn.bind("<Button-1>", self.RR5_vxs2_connections)

        RR5_vxs1_btn = Button(RR5_frame, text='VXS 1 (intelsbshcal1)', width=10, height=6, font='Helvetica 14 bold')
        RR5_vxs1_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR5_vxs1_btn.bind("<Button-1>", self.RR5_vxs1_connections)

    def ftest(self, event):
        print('Test successful!', event.widget.cget('text'))

    def RR4_tdc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Patch Panel Connections")
        win.geometry("900x600")

        RR4_tdc_pp_frames = []                  #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 5                                    #Number of amplifiers in crate.
        nrows = 4                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.

        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)#If you use canvas instead of container frame the 'close' button floats in place while scrolling, but it covers part of the buttons.
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR4_tdc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR4_tdc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR4_tdc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR4_tdc_pp_frames[i], width=4, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'A-'+channel
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR4_DAQ_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'B-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR4_DAQ_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 2 and i!=2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'C-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR4_DAQ_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 3 and i!=2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'D-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR4_DAQ_tdc_pp_info)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR4_tdc_pp_frames[i], text='TDC Patch Panel '+str(10-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            if i==2:
                nrows = 2
            else:
                nrows = 4
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR4_DAQ_tdc_pp_info(self, event):
        text = event.widget.cget('text')
        text = text[:len(text)]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for DAQ TDC Patch Panel '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from front-end TDC patch panel '+str(connections[pmt][6])+' and its output goes to DAQ TDC discriminator '+str(connections[pmt][8])+'.')
                print('This signal terminates at TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR4_tdc_disc_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Discriminator Connections")
        win.geometry("1200x500")
        RR4_tdc_disc_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 25                                   #Number of slots in crate.
        ndisc = 18                                     #Number of summing modules.
        nrows = 8                                    #Number of rows of channels in a module.
        ncols = 2                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)#If you use canvas instead of container frame the 'close' button floats in place while scrolling, but it covers part of the buttons.
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="bottom", fill="x")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR4_tdc_disc_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                for j in range(0,ncols):
                    RR4_tdc_disc_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+3):
                    RR4_tdc_disc_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR4_tdc_disc_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if k == 0:
                            channel = str(j*2+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(18-nmods)+'-'+channel+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR4_DAQ_tdc_disc_info_in)
                            buttons[nmods].append(btn)
                        elif k == 1:
                            channel = str(j*2+2)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(18-nmods)+'-'+channel+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR4_DAQ_tdc_disc_info_in)
                            buttons[nmods].append(btn)
                #Make a button for each discriminator output.
                btn1 = Button(RR4_tdc_disc_frames[i], width=1, height=8, font='Helvetica 8')
                btn1['text'] = 'Ribbon '+str(18-nmods)+' \n Out 1'
                btn1.bind("<Button-1>", self.RR4_DAQ_tdc_disc_info_out)
                buttons[nmods].append(btn1)
                btn2 = Button(RR4_tdc_disc_frames[i], width=1, height=8, font='Helvetica 8')
                btn2['text'] = 'Ribbon '+str(18-nmods)+' \n Out 2'
                btn2.bind("<Button-1>", self.RR4_DAQ_tdc_disc_info_out)
                buttons[nmods].append(btn2)

                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                label = Label(RR4_tdc_disc_frames[i], text='TDC Discriminator '+str(ndisc-nmod))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
            else:
                label = Label(RR4_tdc_disc_frames[i], text='Empty')
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                for j in range(1,5):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')

                for j in range(7,11):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1-2)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                        
                #Add the two output buttons.
                buttons[nmods][16].grid(row=5, column=0, rowspan=1, sticky='NSEW')
                buttons[nmods][17].grid(row=6, column=0, rowspan=1, sticky='NSEW')
                nmods = nmods +1


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR4_DAQ_tdc_disc_info_in(self, event):
        text = event.widget.cget('text')
        if len(text)==12:
            text = text[:8]
        elif len(text)==13:
            text = text[:9]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for DAQ TDC Discriminator '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from DAQ TDC patch panel '+str(connections[pmt][7])+' and its output goes to F1TDC '+str(connections[pmt][9])+'.')
                print('This signal terminates at TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR4_DAQ_tdc_disc_info_out(self, event):
        ribbon2tdc = {'1':'TDC1A','2':'TDC1B','3':'TDC2A','4':'TDC2B','5':'TDC3A','6':'TDC3B','7':'TDC4A','8':'TDC4B','9':'TDC5A','10':'TDC5B','11':'TDC1C','12':'TDC1D','13':'TDC2C','14':'TDC2D','15':'TDC3C','16':'TDC3D','17':'TDC4C','18':'TDC4D'}
        text = event.widget.cget('text')
        if len(text)==17:
            ribbon = text[7:9]
        elif len(text)==16:
            ribbon = text[7:8]
        output = text[-1:]
        #print(ribbon,output)
        if output == '1':
            print('Ribbon cable '+ribbon+' goes to F1TDC input '+ribbon2tdc[ribbon]+'.')
        elif output == '2':
            print('Empty')

    def RR4_adc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Patch Panel Connections")
        win.geometry("800x600")

        RR4_adc_pp_frames = []                  #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 5                                    #Number of amplifiers in crate.
        nrows = 4                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)#If you use canvas instead of container frame the 'close' button floats in place while scrolling, but it covers part of the buttons.
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR4_adc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR4_adc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR4_adc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR4_adc_pp_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(5-i)+'A-'+channel
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR4_DAQ_adc_pp_info)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(5-i)+'B-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR4_DAQ_adc_pp_info)
                        buttons[i].append(btn)
                    elif j == 2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(5-i)+'C-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR4_DAQ_adc_pp_info)
                        buttons[i].append(btn)
                    elif j == 3:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(5-i)+'D-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR4_DAQ_adc_pp_info)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR4_adc_pp_frames[i], text='fADC Patch Panel '+str(5-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR4_DAQ_adc_pp_info(self, event):
        text = event.widget.cget('text')
        text = text[:len(text)]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for DAQ fADC Patch Panel '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from front-end ADC patch panel '+str(connections[pmt][1])+' and its output goes to fADC '+str(connections[pmt][3])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR5_vxs2_connections(self, event):
        win = Tk()
        win.wm_title("RR5 Upper VXS Crate Connections")
        win.geometry("1200x600")#1400x735
        RR5_vxs2_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        ti_buttons = []
        f1_buttons = []
        nslots = 21                                   #Number of slots in crate.
        ndisc = 2                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 1                                     #Number of columns of channels in a module.
        ti_rows = 18
        ti_cols = 1
        nf1 = 5
        f1_rows = 3
        f1_cols = 2
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)#If you use canvas instead of container frame the 'close' button floats in place while scrolling, but it covers part of the buttons.
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="bottom", fill="x")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR5_vxs2_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 18 or i == 19:
                for j in range(0,ncols):
                    RR5_vxs2_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR5_vxs2_frames[i].rowconfigure(j, pad=3, weight=1)

        #Define standalone module rows and cols.
        RR5_vxs2_frames[0].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[0].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[0].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[1].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[1].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[1].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[10].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[10].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[10].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[11].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[11].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[11].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[14].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[14].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[14].rowconfigure(1, pad=3, weight=1)

        for i in range(0,nslots):
            if i == 2 or i == 8  or i == 9  or i == 12  or i == 13  or i == 15 or i == 16 or i == 17:
                RR5_vxs2_frames[i].columnconfigure(0, pad=3, weight=1)
                RR5_vxs2_frames[i].rowconfigure(0, pad=3, weight=1)
                RR5_vxs2_frames[i].rowconfigure(1, pad=3, weight=1)

        #Define TI button grid.
        for i in range(0,ti_rows+5):
            RR5_vxs2_frames[20].rowconfigure(i, pad=3, weight=1)

        for i in range(0,ti_cols):
            RR5_vxs2_frames[20].columnconfigure(i, pad=3, weight=1)

        #Define F1 button grid.
        for i in range(0,nslots):
            if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                for j in range(0,f1_rows+1):
                    RR5_vxs2_frames[i].rowconfigure(j, pad=3, weight=1)

                for j in range(0,f1_cols):
                    RR5_vxs2_frames[i].columnconfigure(j, pad=3, weight=1)

        #Create and bind the fADC buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if  i == 18 or i == 19:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR5_vxs2_frames[i], width=1, height=1, font='Helvetica 8')
                        #Give names and function binds to the different columns.
                        channel = str(15-j)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'f'+str(nmods+17)+'-'+channel+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR5_vxs2_adc_info)
                        buttons[nmods].append(btn)
                nmods = nmods + 1

        #Create and bind the F1 buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                f1_buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,f1_rows):
                    for k in range(0,f1_cols):
                        f1_btn = Button(RR5_vxs2_frames[i], width=4, height=10, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if j == 0 and k == 0:
                            text = 'TDC'+str(nmods+1)+'A\n Chs'+str(17)+'-'+str(32)
                        if j == 0 and k == 1:
                            text = 'TDC'+str(nmods+1)+'B\n Chs'+str(1)+'-'+str(16)
                        if j == 1 and k == 0:
                            text = 'F1 SD\n Cable'
                            f1_btn['height'] = 1
                        if j == 1 and k == 1:
                            text = ''
                            f1_btn['height'] = 1
                            f1_btn['state'] = 'disabled'
                        if j == 2 and k == 0:
                            text = 'TDC'+str(nmods+1)+'C\n Chs'+str(49)+'-'+str(64)
                        if j == 2 and k == 1:
                            text = 'TDC'+str(nmods+1)+'D\n Chs'+str(33)+'-'+str(48)
                        f1_btn['text'] = text
                        f1_btn.bind("<Button-1>", self.RR5_vxs2_tdc_info)
                        f1_buttons[nmods].append(f1_btn)
                nmods = nmods + 1

        #Create a standalone buttons.
        roc_btn = Button(RR5_vxs2_frames[0], text='ROC 17', width=4, height=49, font='Helvetica 8')
        roc_btn.grid(row=1,columnspan=ncols)
        roc_btn.bind("<Button-1>", self.ftest)

        scaler_btn = Button(RR5_vxs2_frames[1], text='F1TDC SD', width=4, height=49, font='Helvetica 8')
        scaler_btn.grid(row=1,columnspan=ncols)
        scaler_btn.bind("<Button-1>", self.ftest)

        vtp_btn = Button(RR5_vxs2_frames[10], text='VTP', width=4, height=49, font='Helvetica 8')
        vtp_btn.grid(row=1,columnspan=ncols)
        vtp_btn.bind("<Button-1>", self.ftest)

        fadc_sd_btn = Button(RR5_vxs2_frames[11], text='fADC SD', width=4, height=49, font='Helvetica 8')
        fadc_sd_btn.grid(row=1,columnspan=ncols)
        fadc_sd_btn.bind("<Button-1>", self.ftest)

        dvcs_pulser_btn = Button(RR5_vxs2_frames[14], text='DVCS\n Pulser', width=4, height=49, font='Helvetica 8')
        dvcs_pulser_btn.grid(row=1,columnspan=ncols)
        dvcs_pulser_btn.bind("<Button-1>", self.ftest)

        #Create buttons for empty slots.
        empty_btns = []
        for i in range(0,nslots):
            if i == 2 or i == 8  or i == 9  or i == 12  or i == 13  or i == 15 or i == 16 or i == 17:
                empty_btn = Button(RR5_vxs2_frames[i], text='Empty', width=4, height=49, font='Helvetica 8')
                empty_btn.grid(row=1,columnspan=ncols)
                empty_btn.bind("<Button-1>", self.ftest)
                empty_btn['state'] = 'disabled'
                empty_btns.append(empty_btn)

        #Make TI buttons.
        ti_btn_names = ['Optical\n Fiber','','','TS#6\n In','TS#5\n In','TS#4\n In','CLK\n In','TS#3\n In','TS#2\n In','TS#1\n In','TRG\n In','INH\n In','OT#2\n Out','OT#1\n Out','O#4\n Out','O#3\n Out','O#2\n Out','O#1\n Out','TRG\n Out','Busy\n Out']

        for i in range(0,ti_rows+2):
            for j in range(0,ti_cols):
                ti_btn = Button(RR5_vxs2_frames[20], width=1, height=1, font='Helvetica 8')
                ti_btn['text'] = ti_btn_names[i]
                ti_btn.bind("<Button-1>", self.ftest)

                if i == 1 or i == 2:
                    ti_btn['state'] = 'disabled'

                ti_buttons.append(ti_btn)

        #Place labels above the fADC buttons in the grid.
        labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 18 or i == 19:
                label = Label(RR5_vxs2_frames[i], text='fADC250 '+str(i-1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1

        #Place labels above the F1 buttons in the grid.
        f1_labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
           if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                f1_label = Label(RR5_vxs2_frames[i], text='F1 TDC '+str(i-2))
                f1_labels.append(f1_label)
                f1_labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
              
        #Place labels above the empty buttons in the grid.
        empty_labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 2 or i == 8  or i == 9  or i == 12  or i == 13  or i == 15 or i == 16 or i == 17:
                empty_label = Label(RR5_vxs2_frames[i], text='Empty')
                empty_labels.append(empty_label)
                empty_labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
                

        #Create standalone labels.
        roc_label = Label(RR5_vxs2_frames[0], text='ROC17')
        roc_label.grid(row=0,columnspan=ncols)

        scaler_label = Label(RR5_vxs2_frames[1], text='F1TDC SD')
        scaler_label.grid(row=0,columnspan=ncols)

        vtp_label = Label(RR5_vxs2_frames[10], text='VTP')
        vtp_label.grid(row=0,columnspan=ncols)

        fadc_sd_label = Label(RR5_vxs2_frames[11], text='fADC SD')
        fadc_sd_label.grid(row=0,columnspan=ncols)

        ti_label = Label(RR5_vxs2_frames[20], text='TI')
        ti_label.grid(row=0,columnspan=ncols)

        dvcs_pulser_label = Label(RR5_vxs2_frames[14], text='DVCS Pulser')
        dvcs_pulser_label.grid(row=0,columnspan=ncols)

        #Place the fADC buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 18 or i == 19:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

        #Place the TI buttons in the grid.
        for i in range(1,ti_rows+3):
            for j in range(0,ti_cols):
                ti_buttons[i-1].grid(row=i, column=j, sticky='NSEW')

        #Place the f1 buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
           if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                for j in range(1,f1_rows+1):
                    for k in range(0,f1_cols):
                        f1_buttons[nmods][(j-1)*f1_cols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR5_vxs2_tdc_info(self, event):
        ribbon2tdc = {'1':'TDC1A','2':'TDC1B','3':'TDC2A','4':'TDC2B','5':'TDC3A','6':'TDC3B','7':'TDC4A','8':'TDC4B','9':'TDC5A','10':'TDC5B','11':'TDC1C','12':'TDC1D','13':'TDC2C','14':'TDC2D','15':'TDC3C','16':'TDC3D','17':'TDC4C','18':'TDC4D'}
        # function to return key for any value
        def get_key(val):
            for key, value in ribbon2tdc.items():
                if val == value:
                    return key
            return 'Empty'
        text = event.widget.cget('text')
        text = text[:5]
        ribbon = get_key(text)
        if ribbon == 'Empty':
            print('Empty')
        else:
            print('The input to '+text+' is ribbon cable '+ribbon+' which comes from DAQ TDC discriminator '+ribbon+'.')

    def RR5_vxs2_adc_info(self, event):
        text = event.widget.cget('text')
        #text = text[:len(text)]
        print(len(text))
        if len(text)==9:
            text = text[:5]
        elif len(text)==10:
            text = text[:6]
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for fADC '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from DAQ ADC patch panel '+str(connections[pmt][2])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR5_vxs1_connections(self, event):
        win = Tk()
        win.wm_title("RR5 Lower VXS Crate Connections")
        win.geometry("1200x600")
        RR5_vxs1_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        ti_buttons = []
        nslots = 21                                   #Number of slots in crate.
        ndisc = 16                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 1                                     #Number of columns of channels in a module.
        ti_rows = 18
        ti_cols = 1
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)#If you use canvas instead of container frame the 'close' button floats in place while scrolling, but it covers part of the buttons.
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="bottom", fill="x")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR5_vxs1_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                for j in range(0,ncols):
                    RR5_vxs1_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR5_vxs1_frames[i].rowconfigure(j, pad=3, weight=1)
        #Define standalone module rows and cols.
        RR5_vxs1_frames[0].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[0].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[0].rowconfigure(1, pad=3, weight=1)

        RR5_vxs1_frames[1].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[1].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[1].rowconfigure(1, pad=3, weight=1)

        RR5_vxs1_frames[10].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[10].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[10].rowconfigure(1, pad=3, weight=1)

        RR5_vxs1_frames[11].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[11].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[11].rowconfigure(1, pad=3, weight=1)

        #Define TI button grid.
        for i in range(0,ti_rows+5):
            RR5_vxs1_frames[20].rowconfigure(i, pad=3, weight=1)

        for i in range(0,ti_cols):
            RR5_vxs1_frames[20].columnconfigure(i, pad=3, weight=1)

        #Create and bind the buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR5_vxs1_frames[i], width=1, height=1, font='Helvetica 8')
                        #Give names and function binds to the different columns.
                        channel = str(15-j)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'f'+str(nmods+1)+'-'+channel+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR5_vxs1_adc_info)
                        buttons[nmods].append(btn)
                nmods = nmods + 1
        #Create a standalone buttons.
        roc_btn = Button(RR5_vxs1_frames[0], text='ROC 16', width=4, height=49, font='Helvetica 8')
        roc_btn.grid(row=1,columnspan=ncols)
        roc_btn.bind("<Button-1>", self.ftest)

        scaler_btn = Button(RR5_vxs1_frames[1], text='Scaler', width=4, height=49, font='Helvetica 8')
        scaler_btn.grid(row=1,columnspan=ncols)
        scaler_btn.bind("<Button-1>", self.ftest)

        vtp_btn = Button(RR5_vxs1_frames[10], text='VTP', width=4, height=49, font='Helvetica 8')
        vtp_btn.grid(row=1,columnspan=ncols)
        vtp_btn.bind("<Button-1>", self.ftest)

        fadc_sd_btn = Button(RR5_vxs1_frames[11], text='fADC SD', width=4, height=49, font='Helvetica 8')
        fadc_sd_btn.grid(row=1,columnspan=ncols)
        fadc_sd_btn.bind("<Button-1>", self.ftest)

        #Make TI buttons.
        ti_btn_names = ['Optical\n Fiber','','','TS#6\n In','TS#5\n In','TS#4\n In','CLK\n In','TS#3\n In','TS#2\n In','TS#1\n In','TRG\n In','INH\n In','OT#2\n Out','OT#1\n Out','O#4\n Out','O#3\n Out','O#2\n Out','O#1\n Out','TRG\n Out','Busy\n Out']

        for i in range(0,ti_rows+2):
            for j in range(0,ti_cols):
                ti_btn = Button(RR5_vxs1_frames[20], width=1, height=1, font='Helvetica 8')
                ti_btn['text'] = ti_btn_names[i]
                ti_btn.bind("<Button-1>", self.ftest)

                if i == 1 or i == 2:
                    ti_btn['state'] = 'disabled'

                ti_buttons.append(ti_btn)

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                label = Label(RR5_vxs1_frames[i], text='fADC250 '+str(nmod+1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
                 
        #Create standalone labels.
        roc_label = Label(RR5_vxs1_frames[0], text='ROC16')
        roc_label.grid(row=0,columnspan=ncols)

        scaler_label = Label(RR5_vxs1_frames[1], text='Scaler')
        scaler_label.grid(row=0,columnspan=ncols)

        vtp_label = Label(RR5_vxs1_frames[10], text='VTP')
        vtp_label.grid(row=0,columnspan=ncols)

        fadc_sd_label = Label(RR5_vxs1_frames[11], text='fADC SD')
        fadc_sd_label.grid(row=0,columnspan=ncols)

        ti_label = Label(RR5_vxs1_frames[20], text='TI')
        ti_label.grid(row=0,columnspan=ncols)

        #Place the fADC buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

        #Place the TI buttons in the grid.
        for i in range(1,ti_rows+3):
            for j in range(0,ti_cols):
                ti_buttons[i-1].grid(row=i, column=j, sticky='NSEW')

        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR5_vxs1_adc_info(self, event):
        text = event.widget.cget('text')
        #text = text[:len(text)]
        print(len(text))
        if len(text)==9:
            text = text[:5]
        elif len(text)==10:
            text = text[:6]
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for fADC '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from DAQ ADC patch panel '+str(connections[pmt][2])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
                filled = 1
        if filled==0:
            print('Empty')

root = Tk()
root.geometry("1200x800")
my_gui = MyFirstGUI(root)
root.mainloop()
