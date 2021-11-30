from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
from tkinter.ttk import Style, Entry
import json
import os
from datetime import date, datetime

#These ranges are set by configuring jumpers on the mpv904 board. +/- 10 V is set by having jumpers on range selectors 2 and 3 only (see mpv904 manual for more settings).
max_volt = 10.0     #Maximum voltage output from mpv904.
min_volt = -10.0    #Minimum voltage output from mpv904.
max_value = 4095    #Maximum value that can be passed to the 12 bit mpv904 to set voltages. 4095 is the lowest voltage value and 0 the highest.
min_value = 0       #Minimum value that can be passed to the 12 bit mpv904 to set voltages. 4095 is the lowest voltage value and 0 the highest.
volt_range = max_volt-min_volt

class DAC_Control_GUI:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("1000x700")
        primary.title("VME Digital to Analog Converter Control GUI")

        global save_file
        save_file = 'mpv904_DAC_Voltage_Settings.json'

        #Read in json file containing a list with all of the saved voltage channel settings.
        voltage_settings = get_voltage_settings(save_file)
        for ch in range(16):
            print('Output channel',ch,'is set to',voltage_settings[ch],'Volts.')

        #Create a frame to hold the exit button.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        info_frame.pack(side='left', fill=BOTH, expand=True)

        save_btn = Button(info_frame, text='Save Settings', width=10, height=4, font='Helvetica 12 bold', command=self.save_settings, bg = "yellow")
        save_btn.pack(side="bottom")

        init_btn = Button(info_frame, text='Initialize\nmpv904', width=10, height=4, font='Helvetica 12 bold', command=initialize_mpv904, bg = "grey")
        init_btn.pack(side="bottom")

        load_btn = Button(info_frame, text='Load mpv904\nLibrary', width=10, height=4, font='Helvetica 12 bold', command=load_library, bg = "grey")
        load_btn.pack(side="bottom")

        #Create a label with the instructions to use the control GUI.
        instructions_label_text = StringVar()
        instructions_label_text.set("This GUI controls the mpv904 VME digital\nto analog converter. The mpv904 is used to set\nseveral remote thresholds for HCal and BigBite\nduring the SBS experiments.\n\nTo set voltages with the mpv904's 16\noutput channels follow these instructions:\n1. On daq@enpcamsonne CODA must be open\nand configured.\n2. Open a terminal on daq@enpcamsonne.\n3. Type 'python3 DAC_Control_GUI.py'\n4. Press the 'Load mpv904 Library' button.\n5. Press the 'Initialize mpv904' button.\n6. Enter a desired voltage for an output channel\nand press the 'Set Voltage' button.\n\n")
        instructions_label = Label(info_frame, width = 40,textvariable=instructions_label_text,font='Helvetica 12 bold')
        instructions_label.pack(side="top")

        #Create a label with the voltage range and resolution.
        range_label_text = StringVar()
        range_label_text.set("The mpv904 VME DAC is currently set to a\nvoltage range of "+str(min_volt)+" V to "+str(max_volt)+" V.\nThe voltage resolution is "+str(round(((max_volt-min_volt)/max_value)*1000.,3))+" mV.")
        range_label = Label(info_frame, width = 40,textvariable=range_label_text,font='Helvetica 12 bold')
        range_label.pack(side="top")

        #canvas = Canvas(info_frame, width = 60,height = 60)
        #my_oval = canvas.create_oval(5, 5, 30, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
        #canvas.itemconfig(my_oval, fill="red")
        #canvas.create_text(17,18,font='Helvetica 8 bold',text='OFF')
        #canvas.pack(side='top')

        #Create frame to set the individual channel voltage outputs.
        ch_settings_frame = Frame(primary, relief=RAISED, borderwidth=2)
        ch_settings_frame.pack(side='left', fill=BOTH, expand=True)

        self.nchs = 16                         #Number of VME DAC output channels.
        nrows = self.nchs                      #Number of rows of channels settings.
        ncols = 6                              #Number of columns per channel setting row.
        self.ch_labels = []                    #Array to hold channel labels.
        self.readback_labels = []              #Array to hold current voltage setting readback labels.
        self.predicted_threshold_labels = []   #Array to hold predicted VME DAC threshold setting.
        self.voltage_setting_entries = []      #Array to hold entry widgets for entering the DAC voltage settings.
        self.power_buttons = []                #Array to hold the power on/off buttons for each output channel.
        self.indicator_lights = []             #Array to hold the indicator lights to see if a channel is energized.

        entry_names = []

        for ch in range(self.nchs):
            #Create and store the channel labels.
            self.ch_label_text = StringVar()
            self.ch_label_text.set("Channel "+str(ch))
            if ch == 0:
                self.ch_label_text.set("HCal Overlapping\nRegion Threshold")
            #elif ch == 4:
            #    self.ch_label_text.set("BBCal Lo\nDiscriminator 1")
            #elif ch == 5:
            #    self.ch_label_text.set("BBCal Lo\nDiscriminator 2")
            elif ch == 6:
                self.ch_label_text.set("BBCal Hi\nDiscriminator 1")
            elif ch == 7:
                self.ch_label_text.set("BBCal Hi\nDiscriminator 2")
            self.ch_label = Label(ch_settings_frame, width = 9,textvariable=self.ch_label_text)
            self.ch_labels.append(self.ch_label)

            #Create and store the channel voltage readback labels.
            self.readback_label_text = StringVar()
            self.readback_label_text.set(str(voltage_settings[ch]))
            self.readback_label = Label(ch_settings_frame, width = 9,textvariable=self.readback_label_text)
            self.readback_labels.append(self.readback_label)

            #Create and store the predicted discriminator threshold labels.
            self.thr_label_text = StringVar()
            self.thr_label_text.set("NA")
            self.thr_label = Label(ch_settings_frame, width = 9,textvariable=self.thr_label_text)
            self.predicted_threshold_labels.append(self.thr_label)

            #Create and store the entry widgets for setting the channel voltage outputs.
            self.voltage_setting_entry = Entry(ch_settings_frame, width = 5)
            self.voltage_setting_entry.insert(0,str(voltage_settings[ch]))
            self.voltage_setting_entries.append(self.voltage_setting_entry)

            #Create and store the buttons for enabling/disabling the channel voltage outputs.
            self.power_btn = Button(ch_settings_frame, text='Set\n Voltage', width=3, height=1, font='Helvetica 8 bold')
            self.power_btn.configure(command=lambda button=self.power_btn, buttons=self.power_buttons, indicator_lights=self.indicator_lights, voltage_setting_entries=self.voltage_setting_entries, readback_labels=self.readback_labels, predicted_threshold_labels=self.predicted_threshold_labels: toggle_power(button,buttons,indicator_lights,voltage_setting_entries,readback_labels,predicted_threshold_labels))
            self.power_buttons.append(self.power_btn)

            #Create and store indicator 'lights' to show whether a channel is energized or not.
            self.canvas = Canvas(ch_settings_frame, width = 40,height = 40)
            my_oval = self.canvas.create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
            self.canvas.itemconfig(my_oval, fill="red")
            self.canvas.create_text(27,18,font='Helvetica 8 bold',text='OFF')
            self.indicator_lights.append(self.canvas)

        #Create labels for the columns.
        col_labels = []
        col_labels_text = ['Output Channels','Current Voltage\n Setting (Volts)','Predicted\n Discriminator\n Threshold (mV)','Enter New\n Voltage Setting\n (Volts)','Set & Save\n Single Channel\n Voltage','Error\n Indicator']
        for col in range(ncols):
            col_label_text = StringVar()
            col_label_text.set(col_labels_text[col])
            col_label = Label(ch_settings_frame, width = 14,height=3,textvariable=col_label_text,font='Helvetica 10 bold')
            col_labels.append(col_label)
            

        #Define the rows and cols for the button grid.
        for row in range(0,nrows+1):
            ch_settings_frame.rowconfigure(row, pad=3, weight=1)
            
        for col in range(0,ncols):
            ch_settings_frame.columnconfigure(col, pad=3, weight=1)

        #Place widgets on the grid.
        for row in range(0,nrows+1):
            if row!=0:
                self.ch_labels[row-1].grid(row=row, column=0, sticky='NSEW')
                self.readback_labels[row-1].grid(row=row, column=1, sticky='NSEW')
                self.predicted_threshold_labels[row-1].grid(row=row, column=2, sticky='NSEW')
                self.voltage_setting_entries[row-1].grid(row=row, column=3)
                self.power_buttons[row-1].grid(row=row, column=4)
                self.indicator_lights[row-1].grid(row=row, column=5)
            for col in range(0,ncols):
                if row==0:
                    col_labels[col].grid(row=row, column=col, sticky='NSEW')

    def save_settings(self):
        voltage_settings = [] #Holds the voltage settings for each of the 16 channels and reads them from the entry boxes when the save button is clicked.
        save_file = 'mpv904_DAC_Voltage_Settings.json'
        #save_file_set_values = open('mpv904_DAC_Voltage_Set_Values.txt','w')
        #save_file_set_values = open('/adaqfs/home/sbs-onl/hcal/cfg/mpv904/mpv904_DAC_Voltage_Set_Values.txt','w')
        save_file_set_values = open('/adaqfs/home/adaq/hcal/gui/mpv904_DAC_Voltage_Set_Values.txt','w')
        for ch in range(self.nchs):
            voltage_settings.append(self.voltage_setting_entries[ch].get())
            print('Output channel',ch,'is now set to',voltage_settings[ch],'Volts.')
            set_value = round( ((float(self.voltage_setting_entries[ch].get())*-1.-min_volt)/volt_range)*max_value )
            print('set value = ',set_value)
            save_file_set_values.write(str(set_value)+' ')
        print('New voltage settings saved to',save_file+'.')
        #Write list to json file.
        out_file = json.dumps(voltage_settings)
        f = open(save_file,"w")
        f.write(out_file)
        f.close()
        save_file_set_values.close()
        
        #Update log file with current settings.
        save_settings_to_log_file(voltage_settings)

#Function to return voltage settings.
def get_voltage_settings(save_file):
    #Read in json file containing a list with all of the saved voltage channel settings.
    with open(save_file) as f: 
        voltage_settings = f.read()
        voltage_settings = json.loads(voltage_settings)
    return voltage_settings

#Function to set the voltage output for individual channels and save the current settings.
def toggle_power(button,buttons,indicator_lights,voltage_setting_entries,readback_labels,predicted_threshold_labels):
    #Check which button was pressed by searching the array holding the power buttons. 
    if button in buttons:
        ch = buttons.index(button)

    #Get the desired voltage setting from the entry box.
    set_voltage = voltage_setting_entries[ch].get()

    #Check that the set voltage is in the allowed DAC voltage range and change nothing if it isn't.
    in_range = volt_in_range(set_voltage)
    if in_range == True:

        #Calculate what value (0-4095) to give the mpv904 output channel closest to that voltage. 4095 is the lowest voltage value and 0 the highest.
        set_value = round( ((float(set_voltage)*-1.-min_volt)/volt_range)*max_value )
        #Calculate the actual voltage output by the DAC.
        real_volt = round(((set_value/max_value)*volt_range+min_volt)*-1,3)
        #Print the actual output in the current voltage column. (Remember this is not a measurement of real voltage!)
        readback_label = readback_labels[ch]
        real_volt_text = StringVar()
        real_volt_text.set(str(real_volt))
        readback_label['textvariable']=real_volt_text

        #Calcultate the expected threshold at the PS706 discriminator (assuming FP threshold screw at min).
        thr = calc_thr(real_volt)
        predicted_threshold_label = predicted_threshold_labels[ch]
        predicted_threshold_text = StringVar()
        predicted_threshold_text.set(str(thr))
        predicted_threshold_label['textvariable']=predicted_threshold_text

        #Send the new voltage to the mpv904 VME DAC.
        os.system('remex hcalROC17 \'mpv904SetDAC(0,'+str(ch)+','+str(set_value)+')\'')

        #Save the voltage setting for the channel just updated.
        voltage_settings = get_voltage_settings(save_file)
        voltage_settings[ch] = set_voltage
        print('Output channel',ch,'is now set to',voltage_settings[ch],'Volts.')
        print('Voltage settings saved to',save_file+'.')
        #Write list to json file.
        out_file = json.dumps(voltage_settings)
        f = open(save_file,"w")
        f.write(out_file)
        f.close()

        save_file_set_values = open('mpv904_DAC_Voltage_Set_Values.txt','w')
        nchs=16
        for ch in range(nchs):
            voltage_settings.append(voltage_setting_entries[ch].get())
            #print('Output channel',ch,'is now set to',voltage_settings[ch],'Volts.')
            set_value = round( ((float(voltage_setting_entries[ch].get())*-1.-min_volt)/volt_range)*max_value )
            #print('set value = ',set_value)
            save_file_set_values.write(str(set_value)+' ')

        save_file_set_values.close()

        #Update log file with current settings.
        save_settings_to_log_file(voltage_settings)

        #text = button.cget('text')
        #off = 'Enable\n Voltage'
        #on = 'Disable\n Voltage'
    
        #Check which button was pressed by searching the array holding the power buttons. 
        #if button in buttons:
        #ch = buttons.index(button)
        #voltage = voltage_settings[ch]
        #print('Channel',ch,'is set to',voltage,'Volts.')

        #if text==off:
        #    print('Turning channel',ch,'on.')
        #    button['text']=on
        #    my_oval = indicator_lights[ch].create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
        #    indicator_lights[ch].itemconfig(my_oval, fill="green")
        #    indicator_lights[ch].create_text(27,18,font='Helvetica 8 bold',text='ON')
        
        #if text==on:
        #    print('Turning channel',ch,'off.')
        #    button['text']=off
        #    my_oval = indicator_lights[ch].create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
        #    indicator_lights[ch].itemconfig(my_oval, fill="red")
        #    indicator_lights[ch].create_text(27,18,font='Helvetica 8 bold',text='OFF')

#Function to save current settings to log file.
def save_settings_to_log_file(voltage_settings):
    #Update log file.
    log_file_name = 'mpv904_DAC_Voltage_Settings_Log_File.txt'
    log_file = open(log_file_name,"a")#append mode
    now = datetime.now()# datetime object containing current date and time
    dt_string = now.strftime("%B %d, %Y, %H:%M:%S")
    text = dt_string+":\nNew Voltage Settings = "+str(voltage_settings)+'\n'
    log_file.write(text)
    log_file.close()

def load_library():
    #Loads the mpv904 library.
    os.system('remex -d hcalROC17 linuxvme/mpv904/libmpv904.so')

def initialize_mpv904():
    #Initializes the mpv904.
    os.system('remex hcalROC17 \'mpv904Init(0x700000,0,1,0x6)\'')

def calc_thr(set_volt):
    thr = round(97.9853*set_volt-9.8114,1)
    return thr

def volt_in_range(set_volt):
    #Check if the set voltage is in the range of allowed voltages.
    if float(set_volt) <= max_volt and float(set_volt) >= min_volt:
        return True
    else:
        print('ERROR: Voltage setting outside of voltage range ('+str(min_volt)+' to '+str(max_volt)+'). Voltage was not changed.')
        return False

root = Tk()
my_gui = DAC_Control_GUI(root)
root.mainloop()
