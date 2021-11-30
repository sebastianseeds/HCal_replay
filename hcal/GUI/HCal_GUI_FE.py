from tkinter import Tk, Label, Button, StringVar, W, E, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar, Toplevel
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

#Fill PMT module numbers.
for i in range(0,channels):
    pmt_mods.append(i+1)
    #print(pmt_mods[i])

#Fill a list with the names for the buttons.
button_names = []
for i in range(0,channels):
    button_names.append("Button"+str(i+1))
#print(button_names)

class MyFirstGUI:
    def __init__(self, primary):
        self.primary = primary
        primary.title("Front-End Electronics")

        self.style = Style()
        self.style.theme_use("default")

        #Create frames to hold the contents of each FE rack.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        RR1_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR1_frame.pack(side='left', fill=BOTH, expand=True)

        RR2_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR2_frame.pack(side='left', fill=BOTH, expand=True)

        RR3_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR3_frame.pack(side='left', fill=BOTH, expand=True)

        #Create buttons and labels for electronics in RR1
        RR1_label = Label(RR1_frame, text='RR1', font='Helvetica 16 bold')
        RR1_label.pack(side="top")

        RR1_sum_mods_btn = Button(RR1_frame, text='Summing Modules', width=10, height=5, font='Helvetica 14 bold')
        RR1_sum_mods_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR1_sum_mods_btn.bind("<Button-1>", self.RR1_sum_connections)

        RR1_splitter_btn = Button(RR1_frame, text='Splitter Panels', width=10, height=13, font='Helvetica 14 bold')
        RR1_splitter_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR1_splitter_btn.bind("<Button-1>", self.RR1_splitter_connections)

        RR1_amp_btn = Button(RR1_frame, text='Amplifiers', width=10, height=5, font='Helvetica 14 bold')
        RR1_amp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR1_amp_btn.bind("<Button-1>", self.RR1_amp_connections)

        #Create buttons and labels for electronics in RR2
        RR2_label = Label(RR2_frame, text='RR2', font='Helvetica 16 bold')
        RR2_label.pack(side="top")

        RR2_upper_tdc_pp_btn = Button(RR2_frame, text='Upper TDC Patch Panels', width=10, height=4, font='Helvetica 14 bold')
        RR2_upper_tdc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_upper_tdc_pp_btn.bind("<Button-1>", self.RR2_upper_tdc_pp_connections)

        RR2_upper_tdc_disc_btn = Button(RR2_frame, text='TDC Discriminators from RR3', width=10, height=5, font='Helvetica 14 bold')
        RR2_upper_tdc_disc_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_upper_tdc_disc_btn.bind("<Button-1>", self.RR2_Upper_TDC_Disc_connections)

        RR2_lower_tdc_pp_btn = Button(RR2_frame, text='Lower TDC Patch Panels', width=10, height=6, font='Helvetica 14 bold')
        RR2_lower_tdc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_lower_tdc_pp_btn.bind("<Button-1>", self.RR2_lower_tdc_pp_connections)

        RR2_lower_tdc_disc_btn = Button(RR2_frame, text='TDC Discriminators from RR1', width=10, height=5, font='Helvetica 14 bold')
        RR2_lower_tdc_disc_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_lower_tdc_disc_btn.bind("<Button-1>", self.RR2_Lower_TDC_Disc_connections)

        RR2_fadc_pp_btn = Button(RR2_frame, text='fADC Patch Panels', width=10, height=11, font='Helvetica 14 bold')
        RR2_fadc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_fadc_pp_btn.bind("<Button-1>", self.RR2_adc_pp_connections)

        #Create buttons and labels for electronics in RR3
        RR3_label = Label(RR3_frame, text='RR3', font='Helvetica 16 bold')
        RR3_label.pack(side="top")

        RR3_sum_mods_btn = Button(RR3_frame, text='Summing Modules', width=10, height=5, font='Helvetica 14 bold')
        RR3_sum_mods_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR3_sum_mods_btn.bind("<Button-1>", self.RR3_sum_connections)

        RR3_splitter_btn = Button(RR3_frame, text='Splitter Panels', width=10, height=13, font='Helvetica 14 bold')
        RR3_splitter_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR3_splitter_btn.bind("<Button-1>", self.RR3_splitter_connections)

        RR3_amp_btn = Button(RR3_frame, text='Amplifiers', width=10, height=5, font='Helvetica 14 bold')
        RR3_amp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR3_amp_btn.bind("<Button-1>", self.RR3_amp_connections)

    def fconnections(self, event):
        pmt = event.widget.cget('text')
        print('PMT Module ',pmt,': amplifier channel = ', connections[pmt][0],', fADC channel = ', connections[pmt][1], ', HV channel = ',connections[pmt][2])
        #print(event.widget.cget('text'))

    def RR3_amp_connections(self, event):
        win = Toplevel()
        win.wm_title("RR3 Amplifier Connections")
        win.geometry("1400x600")

        RR3_amp_frames = []                           #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 16                                    #Number of rows of channels in an amp.
        ncols = 3                                     #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        exit_frame = Frame(win, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        #exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', bg = "red")
        #exit_btn.bind("<Button-1>",self.quit_out)
        exit_btn.pack(side="bottom")

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR3_amp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR3_amp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR3_amp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR3_amp_frames[i], width=1, height=1, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if k == 0:
                        channel = str(j+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'a'+str(i+1)+'-'+channel+'\n In'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR3_amp_info_in)
                        buttons[i].append(btn)
                    elif k == 1:
                        channel = str(j+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'a'+str(i+1)+'-'+channel+'\n Out 1'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR3_amp_info_out1)
                        buttons[i].append(btn)
                    elif k == 2:
                        channel = str(j+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'a'+str(i+1)+'-'+channel+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR3_amp_info_out2)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the amp buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR3_amp_frames[i], text='Amplifier '+str(i+10))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the amp buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')

    def RR3_amp_info_in(self, event):
        text = event.widget.cget('text')
        text = text[:5]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' In'
                print('******************** Information for Amplifier '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from PMT '+pmt+' ('+connections[pmt][12]+') and its output goes to front-end fADC patch panel '+str(connections[pmt][1])+' and splitter panel '+str(connections[pmt][4])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+', F1TDC '+str(connections[pmt][9])+', and summing module '+str(connections[pmt][10])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR3_amp_info_out1(self, event):
        text = event.widget.cget('text')
        text = text[:5]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 1'
                print('******************** Information for Amplifier '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from PMT '+pmt+' ('+connections[pmt][12]+') and its output goes to splitter panel '+str(connections[pmt][4])+'.')
                print('This signal terminates at F1TDC '+str(connections[pmt][9])+' and summing module '+str(connections[pmt][10])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR3_amp_info_out2(self, event):
        text = event.widget.cget('text')
        text = text[:5]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 2'
                print('******************** Information for Amplifier '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from PMT '+pmt+' ('+connections[pmt][12]+') and its output goes to front-end fADC patch panel '+str(connections[pmt][1])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')

    def ftest(self, event):
        print('Test successful!', event.widget.cget('text'))

    def RR1_amp_connections(self, event):
        win = Tk()
        win.wm_title("RR1 Amplifier Connections")
        win.geometry("1400x600")

        RR1_amp_frames = []                           #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 16                                    #Number of rows of channels in an amp.
        ncols = 3                                     #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        exit_frame = Frame(win, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR1_amp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR1_amp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR1_amp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR1_amp_frames[i], width=1, height=1, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if k == 0:
                        channel = str(j+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'b'+str(9-i)+'-'+channel+'\n In'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR1_amp_info_in)
                        buttons[i].append(btn)
                    elif k == 1:
                        channel = str(j+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'b'+str(9-i)+'-'+channel+'\n Out 1'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR1_amp_info_out1)
                        buttons[i].append(btn)
                    elif k == 2:
                        channel = str(j+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'b'+str(9-i)+'-'+channel+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR1_amp_info_out2)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the amp buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR1_amp_frames[i], text='Amplifier '+str(9-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the amp buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')

    def RR1_amp_info_in(self, event):
        text = event.widget.cget('text')
        text = text[:5]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' In'
                print('******************** Information for Amplifier '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from PMT '+pmt+' ('+connections[pmt][12]+') and its output goes to front-end fADC patch panel '+str(connections[pmt][1])+' and splitter panel '+str(connections[pmt][4])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+', F1TDC '+str(connections[pmt][9])+', and summing module '+str(connections[pmt][10])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR1_amp_info_out1(self, event):
        text = event.widget.cget('text')
        text = text[:5]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 1'
                print('******************** Information for Amplifier '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from PMT '+pmt+' ('+connections[pmt][12]+') and its output goes to splitter panel '+str(connections[pmt][4])+'.')
                print('This signal terminates at F1TDC '+str(connections[pmt][9])+' and summing module '+str(connections[pmt][10])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR1_amp_info_out2(self, event):
        text = event.widget.cget('text')
        text = text[:5]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 2'
                print('******************** Information for Amplifier '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from PMT '+pmt+' ('+connections[pmt][12]+') and its output goes to front-end fADC patch panel '+str(connections[pmt][1])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')

    def RR3_sum_connections(self, event):
        win = Tk()
        win.wm_title("RR3 Summing Module Connections")
        win.geometry("1400x600")

        RR3_sum_frames = []                           #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        nsum_mods = 5                                 #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 2                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        exit_frame = Frame(win, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create 9 different frames to hold the sum channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR3_sum_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(0,ncols):
                    RR3_sum_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+5):
                    RR3_sum_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows+4):
                    for k in range(0,ncols):
                        btn = Button(RR3_sum_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if j<nrows:
                            if k == 0:
                                channel = str(j+1)
                                if len(channel)==1:
                                    channel = channel[:0]+'0'+channel[0:]
                                text = 'Sum'+str(10-nmods)+'B-'+channel+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR3_sum_info_in)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                channel = str(j+1)
                                if len(channel)==1:
                                    channel = channel[:0]+'0'+channel[0:]
                                text = 'Sum'+str(10-nmods)+'A-'+channel+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR3_sum_info_in)
                                buttons[nmods].append(btn)
                        #Create different buttons for the outputs.
                        else:
                            if k == 0:
                                text = 'Sum'+str(10-nmods)+'B \n Out '+str(j-nrows+1)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR3_sum_info_out)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                text = 'Sum'+str(10-nmods)+'A \n Out '+str(j-nrows+1)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR3_sum_info_out)
                                buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                label = Label(RR3_sum_frames[i], text='Summing Module '+str(10-nlabels))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(1,nrows+5):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

    def RR3_sum_info_in(self, event):
        text = event.widget.cget('text')
        #print(len(text))
        if len(text)==12:
            text = text[:8]
        elif len(text)==13:
            text = text[:9]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' In'
                print('******************** Information for Summing Module '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from splitter panel '+connections[pmt][4]+' and its output goes the overlapping region trigger.')
                print('This signal terminates at summing module '+str(connections[pmt][10])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR3_sum_info_out(self, event):
        text = event.widget.cget('text')
        #print(len(text))
        if len(text)==13:
            sum_mod = text[3]
            side = text[4]
            if side=='A':
                cluster = int(sum_mod)*2-2
            elif side == 'B':
                cluster = int(sum_mod)*2-1
        elif len(text)==14:
            sum_mod = text[3:5]
            side = text[5]
            if side=='A':
                cluster = int(sum_mod)*2-2
            elif side == 'B':
                cluster = int(sum_mod)*2-1
        if sum_mod=='10' and side=='B':
            print('Empty')
        else:
            print('Sum'+sum_mod+side+' connects to cluster '+str(cluster)+' on HCal')

    def RR1_sum_connections(self, event):
        win = Tk()
        win.wm_title("RR1 Summing Module Connections")
        win.geometry("1400x600")

        RR1_sum_frames = []                           #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        nsum_mods = 5                                 #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 2                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        exit_frame = Frame(win, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create 9 different frames to hold the sum channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR1_sum_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(0,ncols):
                    RR1_sum_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+5):
                    RR1_sum_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows+4):
                    for k in range(0,ncols):
                        btn = Button(RR1_sum_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if j<nrows:
                            if k == 0:
                                channel = str(j+1)
                                if len(channel)==1:
                                    channel = channel[:0]+'0'+channel[0:]
                                text = 'Sum'+str(nmods+1)+'A-'+channel+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR1_sum_info_in)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                channel = str(j+1)
                                if len(channel)==1:
                                    channel = channel[:0]+'0'+channel[0:]
                                text = 'Sum'+str(nmods+1)+'B-'+channel+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR1_sum_info_in)
                                buttons[nmods].append(btn)
                        #Create different buttons for the outputs.
                        else:
                            if k == 0:
                                text = 'Sum'+str(nmods+1)+'A \n Out'+str(j-nrows)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR1_sum_info_out)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                text = 'Sum'+str(nmods+1)+'B \n Out'+str(j-nrows)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.RR1_sum_info_out)
                                buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                label = Label(RR1_sum_frames[i], text='Summing Module '+str(nlabels+1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(1,nrows+5):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

    def RR1_sum_info_in(self, event):
        text = event.widget.cget('text')
        #print(len(text))
        if len(text)==12:
            text = text[:8]
        elif len(text)==13:
            text = text[:9]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' In'
                print('******************** Information for Summing Module '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from splitter panel '+connections[pmt][4]+' and its output goes the overlapping region trigger.')
                print('This signal terminates at summing module '+str(connections[pmt][10])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR1_sum_info_out(self, event):
        text = event.widget.cget('text')
        #print(len(text))
        if len(text)==12:
            sum_mod = text[3]
            side = text[4]
            if side=='A':
                cluster = int(sum_mod)*2-1
            elif side == 'B':
                cluster = int(sum_mod)*2
        if sum_mod=='5' and side=='B':
            print('Empty')
        else:
            print('Sum'+sum_mod+side+' connects to cluster '+str(cluster)+' on HCal')

    def RR2_Upper_TDC_Disc_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Discriminator Connections")
        win.geometry("1400x517")
        RR2_upper_TDC_disc_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        ndisc = 9                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 3                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
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
            RR2_upper_TDC_disc_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(0,ncols):
                    RR2_upper_TDC_disc_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR2_upper_TDC_disc_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR2_upper_TDC_disc_frames[i], width=4, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if k == 0:
                            channel = str(j+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(18-nmods)+'-'+channel+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR2_FE_upper_tdc_disc_info)
                            buttons[nmods].append(btn)
                        elif k == 1:
                            channel = str(j+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(18-nmods)+'-'+channel+'\n Out 1'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR2_FE_upper_tdc_disc_info)
                            buttons[nmods].append(btn)
                        elif k == 2:
                            channel = str(j+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(18-nmods)+'-'+channel+'\n Out 2'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR2_FE_upper_tdc_disc_info)
                            buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                label = Label(RR2_upper_TDC_disc_frames[i], text='TDC Discriminator '+str(18-nlabels))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR2_FE_upper_tdc_disc_info(self, event):
        text = event.widget.cget('text')
        text = text[:9]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for Front-End Lower TDC Discriminator '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from splitter panel '+str(connections[pmt][4])+' and its output goes to DAQ TDC patch panel '+str(connections[pmt][7])+'.')
                print('This signal terminates at TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR2_Lower_TDC_Disc_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Lower TDC Discriminator Connections")
        win.geometry("1400x517")
        RR2_lower_TDC_disc_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        ndisc = 9                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 3                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
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
            RR2_lower_TDC_disc_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(0,ncols):
                    RR2_lower_TDC_disc_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR2_lower_TDC_disc_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR2_lower_TDC_disc_frames[i], width=3, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if k == 0:
                            channel = str(j+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(nmods+1)+'-'+channel+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR2_FE_lower_tdc_disc_info)
                            buttons[nmods].append(btn)
                        elif k == 1:
                            channel = str(j+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(nmods+1)+'-'+channel+'\n Out 1'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR2_FE_lower_tdc_disc_info)
                            buttons[nmods].append(btn)
                        elif k == 2:
                            channel = str(j+1)
                            if len(channel)==1:
                                channel = channel[:0]+'0'+channel[0:]
                            text = 'Disc'+str(nmods+1)+'-'+channel+'\n Out 2'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.RR2_FE_lower_tdc_disc_info)
                            buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                label = Label(RR2_lower_TDC_disc_frames[i], text='TDC Discriminator '+str(nlabels+1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR2_FE_lower_tdc_disc_info(self, event):
        text = event.widget.cget('text')
        text = text[:8]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for Front-End Lower TDC Discriminator '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from splitter panel '+str(connections[pmt][4])+' and its output goes to DAQ TDC patch panel '+str(connections[pmt][7])+'.')
                print('This signal terminates at TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR3_splitter_connections(self, event):
        win = Tk()
        win.wm_title("RR3 Splitter Panel Connections")
        win.geometry("800x600")

        RR3_splitter_frames = []                      #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 3                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR3_splitter_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR3_splitter_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR3_splitter_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR3_splitter_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'SP'+str(18-i)+'-'+channel+'\n Out 1'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR3_splitter_info_out1)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'SP'+str(18-i)+'-'+channel+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR3_splitter_info_in)
                        buttons[i].append(btn)
                    elif j == 2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'SP'+str(18-i)+'-'+channel+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR3_splitter_info_out2)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR3_splitter_frames[i], text='Splitter Panel '+str(18-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR3_splitter_info_in(self, event):
        text = event.widget.cget('text')
        text = text[:7]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' In'
                print('******************** Information for Splitter '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+connections[pmt][0]+' Out 1 and its output goes to front-end TDC discriminator '+str(connections[pmt][5])+' and summing module '+str(connections[pmt][10])+'.')
                print('This signal terminates at F1TDC '+str(connections[pmt][9])+' and summing module '+str(connections[pmt][10])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR3_splitter_info_out1(self, event):
        text = event.widget.cget('text')
        text = text[:7]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 1'
                print('******************** Information for Splitter '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+connections[pmt][0]+' Out 1 and its output goes to summing module '+str(connections[pmt][10])+'.')
                print('This signal terminates at summing module '+str(connections[pmt][10])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR3_splitter_info_out2(self, event):
        text = event.widget.cget('text')
        text = text[:7]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 2'
                print('******************** Information for Splitter '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+connections[pmt][0]+' Out 1 and its output goes to front-end TDC discriminator '+str(connections[pmt][5])+'.')
                print('This signal terminates at F1TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')

    def RR1_splitter_connections(self, event):
        win = Tk()
        win.wm_title("RR1 Splitter Panel Connections")
        win.geometry("800x600")

        RR1_splitter_frames = []                      #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 3                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=win.quit, bg = "red")
        exit_btn.pack(side="bottom")

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR1_splitter_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR1_splitter_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR1_splitter_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR1_splitter_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'SP'+str(nslots-i)+'-'+channel+'\n Out 1'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR1_splitter_info_out1)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'SP'+str(nslots-i)+'-'+channel+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR1_splitter_info_in)
                        buttons[i].append(btn)
                    elif j == 2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'SP'+str(nslots-i)+'-'+channel+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR1_splitter_info_out2)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR1_splitter_frames[i], text='Splitter Panel '+str(nslots-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR1_splitter_info_in(self, event):
        text = event.widget.cget('text')
        text = text[:6]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' In'
                print('******************** Information for Splitter '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+connections[pmt][0]+' Out 1 and its output goes to front-end TDC discriminator '+str(connections[pmt][5])+' and summing module '+str(connections[pmt][10])+'.')
                print('This signal terminates at F1TDC '+str(connections[pmt][9])+' and summing module '+str(connections[pmt][10])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR1_splitter_info_out1(self, event):
        text = event.widget.cget('text')
        text = text[:6]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 1'
                print('******************** Information for Splitter '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+connections[pmt][0]+' Out 1 and its output goes to summing module '+str(connections[pmt][10])+'.')
                print('This signal terminates at summing module '+str(connections[pmt][10])+'.')
                print('The summing module data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')

    def RR1_splitter_info_out2(self, event):
        text = event.widget.cget('text')
        text = text[:6]
        #print(text)
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                text = text + ' Out 2'
                print('******************** Information for Splitter '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+connections[pmt][0]+' Out 1 and its output goes to front-end TDC discriminator '+str(connections[pmt][5])+'.')
                print('This signal terminates at F1TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')

    def RR2_adc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 fADC Patch Panel Connections")
        win.geometry("800x600")

        RR2_adc_pp_frames = []                        #List to hold the amplifier frames.
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

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
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
            RR2_adc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR2_adc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR2_adc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR2_adc_pp_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i)+'A-'+channel
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR2_FE_adc_pp_info)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i)+'B-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_adc_pp_info)
                        buttons[i].append(btn)
                    elif j == 2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i)+'C-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_adc_pp_info)
                        buttons[i].append(btn)
                    elif j == 3:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i)+'D-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_adc_pp_info)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR2_adc_pp_frames[i], text='fADC Patch Panel '+str(nslots-i))
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

    def RR2_FE_adc_pp_info(self, event):
        text = event.widget.cget('text')
        text = text[:7]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for Front-End fADC Patch Panel '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from amplifier '+str(connections[pmt][0])+' and its output goes to DAQ fADC patch panel '+str(connections[pmt][2])+'.')
                print('This signal terminates at fADC '+str(connections[pmt][3])+'.')
                print('The fADC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') --> amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR2_lower_tdc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Lower TDC Patch Panel Connections")
        win.geometry("800x500")

        RR2_lower_tdc_pp_frames = []                  #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 3                                    #Number of amplifiers in crate.
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

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
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
            RR2_lower_tdc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR2_lower_tdc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR2_lower_tdc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR2_lower_tdc_pp_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i+5)+'A-'+channel
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR2_FE_lower_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i+5)+'B-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_lower_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 2 and i!=0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i+5)+'C-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_lower_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 3 and i!=0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(nslots-i+5)+'D-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_lower_tdc_pp_info)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR2_lower_tdc_pp_frames[i], text='Lower TDC Patch Panel '+str(nslots-i+5))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            if i==0:
                nrows = 2
            else:
                nrows = 4
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR2_FE_lower_tdc_pp_info(self, event):
        text = event.widget.cget('text')
        text = text[:7]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for Front-End Lower TDC Patch Panel '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from front-end TDC discriminator '+str(connections[pmt][5])+' and its output goes to DAQ TDC patch panel '+str(connections[pmt][7])+'.')
                print('This signal terminates at TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                filled = 1
        if filled==0:
            print('Empty')

    def RR2_upper_tdc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Patch Panel Connections")
        win.geometry("900x395")

        RR2_upper_tdc_pp_frames = []                  #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 2                                    #Number of amplifiers in crate.
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

        exit_frame = Frame(container_frame, relief=RAISED, borderwidth=2)
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
            RR2_upper_tdc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR2_upper_tdc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR2_upper_tdc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR2_upper_tdc_pp_frames[i], width=4, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'A-'+channel
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.RR2_FE_upper_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 1:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'B-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_upper_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 2:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'C-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_upper_tdc_pp_info)
                        buttons[i].append(btn)
                    elif j == 3:
                        channel = str(k+1)
                        if len(channel)==1:
                            channel = channel[:0]+'0'+channel[0:]
                        text = 'PP'+str(10-i)+'D-'+channel
                        btn['text'] = text
                        btn.bind("<Button-1>", self.RR2_FE_upper_tdc_pp_info)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR2_upper_tdc_pp_frames[i], text='Upper TDC Patch Panel '+str(10-i))
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

    def RR2_FE_upper_tdc_pp_info(self, event):
        text = event.widget.cget('text')
        text = text[:len(text)]
        #print(text)
        filled = 0
        for pmt in range(1,289):
            pmt = str(pmt)
            if text in connections[pmt]:
                print('******************** Information for Front-End Upper TDC Patch Panel '+text+' ********************')
                print(text+' connects to HV channel '+str(connections[pmt][11])+'.')
                print(text+'\'s input comes from front-end TDC discriminator '+str(connections[pmt][5])+' and its output goes to DAQ TDC patch panel '+str(connections[pmt][7])+'.')
                print('This signal terminates at TDC '+str(connections[pmt][9])+'.')
                print('The TDC data flow follows: PMT '+pmt+' ('+connections[pmt][12]+') -->  amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
                filled = 1
        if filled==0:
            print('Empty')

root = Tk()
root.geometry("1200x800")
my_gui = MyFirstGUI(root)
root.mainloop()
