from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar, Radiobutton
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

rate = 0
rate_text = ['1','2','5','10','20','50','100','200','500','1000','2000','5000','10000','20000','50000','100000']

class DVCS_Pulser_Control_GUI:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("800x700")
        primary.title("DVCS Pulser Control GUI")

        global save_file
        #save_file = 'mpv904_DAC_Voltage_Settings.json'
        save_file = 'DVCS_Pulser_Settings.json'

        #Read in json file containing a list with all of the saved voltage channel settings.
        #voltage_settings = get_voltage_settings(save_file)
        #for ch in range(16):
            #print('Output channel',ch,'is set to',voltage_settings[ch],'Volts.')
        
        saved_settings = get_saved_settings(save_file)
        print('The DVCS pulser rate was last set to '+str(saved_settings)+' Hz.')

        #Create a frame to hold the exit button.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        info_frame.pack(side='left', fill=BOTH, expand=True)

        save_btn = Button(info_frame, text='Save Settings', width=10, height=4, font='Helvetica 12 bold', command=save_settings, bg = "yellow")
        save_btn.pack(side="bottom")

        enable_chs_btn = Button(info_frame, text='Enable Output\nChannels', width=10, height=4, font='Helvetica 12 bold', command=enable_chs, bg = "grey")
        enable_chs_btn.pack(side="bottom")

        init_btn = Button(info_frame, text='Initialize\nDVCS Pulser', width=10, height=4, font='Helvetica 12 bold', command=initialize, bg = "grey")
        init_btn.pack(side="bottom")

        load_btn = Button(info_frame, text='Load DVCS\nPulser Library', width=10, height=4, font='Helvetica 12 bold', command=load_library, bg = "grey")
        load_btn.pack(side="bottom")

        #Create a label with the instructions to use the control GUI.
        instructions_label_text = StringVar()
        instructions_label_text.set("This GUI controls the DVCS pulser which\nis used to pulse the LEDs in HCal\nfor gain monitoring purposes. This pulser is\nused in conjunction with a NIM gate generator\nto produce logical pulses compatible with\nthe HCal LED power boxes.\n\n")
        instructions_label = Label(info_frame, width = 40,textvariable=instructions_label_text,font='Helvetica 12 bold')
        instructions_label.pack(side="top")

        #Create a label with the voltage range and resolution.
        range_label_text = StringVar()
        range_label_text.set("To operate the DVCS pulser first have CODA\nopen and configured on the machine running\nthis GUI. Then download the library ")
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
        ncols = 1                              #Number of columns per channel setting row.
        self.ch_labels = []                    #Array to hold channel labels.
        self.readback_labels = []              #Array to hold current voltage setting readback labels.
        self.predicted_threshold_labels = []   #Array to hold predicted VME DAC threshold setting.
        self.voltage_setting_entries = []      #Array to hold entry widgets for entering the DAC voltage settings.
        self.power_buttons = []                #Array to hold the power on/off buttons for each output channel.
        self.indicator_lights = []             #Array to hold the indicator lights to see if a channel is energized.
        self.rate_radios = []

        #rate_text = ['1','2','5','10','20','50','100','200','500','1000','2000','5000','10000','20000','50000','100000']
        #rate = IntVar()

        for ch in range(self.nchs):
            #Create and store the channel labels.
            self.ch_label_text = StringVar()
            self.ch_label_text.set(rate_text[ch])
            self.ch_label = Label(ch_settings_frame, width = 9,textvariable=self.ch_label_text)
            self.ch_labels.append(self.ch_label)

            #self.rate_radio = Radiobutton(ch_settings_frame, text=rate_text[ch], padx = 20, variable=rate, command=set_rate, value=ch)
            self.rate_radio = Radiobutton(ch_settings_frame, text=rate_text[ch], padx = 20, variable=rate, value=ch,indicatoron = 0)
            self.rate_radio.configure(command=lambda radio=self.rate_radio, radios=self.rate_radios: set_rate(radio,radios))
            self.rate_radios.append(self.rate_radio)

            #Create and store indicator 'lights' to show whether a channel is energized or not.
            self.canvas = Canvas(ch_settings_frame, width = 40,height = 40)
            my_oval = self.canvas.create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
            self.canvas.itemconfig(my_oval, fill="red")
            self.canvas.create_text(27,18,font='Helvetica 8 bold',text='OFF')
            self.indicator_lights.append(self.canvas)

        #Create labels for the columns.
        col_labels = []
        col_labels_text = ['Set Pulser\nRate (Hz)']
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
                self.rate_radios[row-1].grid(row=row, column=0, sticky='NSEW')
            for col in range(0,ncols):
                if row==0:
                    col_labels[col].grid(row=row, column=col, sticky='NSEW')

#Function to return voltage settings.
def get_voltage_settings(save_file):
    #Read in json file containing a list with all of the saved voltage channel settings.
    with open(save_file) as f: 
        voltage_settings = f.read()
        voltage_settings = json.loads(voltage_settings)
    return voltage_settings

#Function to return saved settings.
def get_saved_settings(save_file):
    #Read in json file containing a list with all of the saved setting.
    with open(save_file) as f: 
        saved_settings = f.read()
        saved_settings = json.loads(saved_settings)
    return saved_settings

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
def save_settings_to_log_file(saved_settings):
    #Update log file.
    log_file_name = 'DVCS_Pulser_Settings_Log_File.txt'
    log_file = open(log_file_name,"a")#append mode
    now = datetime.now()# datetime object containing current date and time
    dt_string = now.strftime("%B %d, %Y, %H:%M:%S")
    text = dt_string+":\nNew DVCS Pulser Settings = "+str(saved_settings)+' Hz.\n'
    log_file.write(text)
    log_file.close()

def load_library():
    #Loads the DVCS pulser library.
    os.system('remex -d hcalROC17 linuxvme/dvcsPulser/libdvcsPulser.so')

def initialize():
    #Initializes the DVCS pulser.
    os.system('remex hcalROC17 \'dvcsPulserInit(0xEF0,0)\'')

def enable_chs():
    #Enables the output channels for the DVCS pulser.
    os.system('remex hcalROC17 \'dvcsPulserEnable(0x1d)\'')

def save_settings():
    print('Save Test')

def calc_thr(set_volt):
    thr = round(97.9853*set_volt-9.8114,1)
    return thr

def volt_in_range(set_volt):
    #Check if the set voltage is in the range of allowed voltages.
    if float(set_volt) < max_volt and float(set_volt) > min_volt:
        return True
    else:
        print('ERROR: Voltage setting outside of voltage range ('+str(min_volt)+' to '+str(max_volt)+'). Voltage was not changed.')
        return False

def set_rate(radio,radios):
    if radio in radios:
        idx = radios.index(radio)

    rate = int(rate_text[idx])
    #print(idx,rate_text[idx])

    os.system('remex hcalROC17 \'dvcsPulserSetFrequency('+str(rate)+')\'')
    print('The DVCS pulser rate is now set to '+str(rate)+' Hz.')

    saved_settings = rate

    #Write list to json file.
    out_file = json.dumps(saved_settings)
    f = open(save_file,"w")
    f.write(out_file)
    f.close()
    
    save_settings_to_log_file(saved_settings)

root = Tk()
my_gui = DVCS_Pulser_Control_GUI(root)
root.mainloop()
