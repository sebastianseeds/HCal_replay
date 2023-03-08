from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
from tkinter.ttk import Style, Entry
import json

nrows = 24
ncols = 12
channels = nrows * ncols
pmt_mods = []    #Holds PMT module numbers.
buttons = []

#Dictionary containing all of the relevant connections for each PMT.
# reading the data from the file 
with open('hcal_connections.json') as f: 
    connections = f.read()
#print(type(connections))
# reconstructing the data as a dictionary 
connections = json.loads(connections)
#print(type(connections))
#print(connections['1'][0])

clusters = [[1, 2, 3, 4, 13, 14, 15, 16, 25, 26, 27, 28, 37, 38, 39, 40], [5, 6, 7, 8, 17, 18, 19, 20, 29, 30, 31, 32, 41, 42, 43, 44], [9, 10, 11, 12, 21, 22, 23, 24, 33, 34, 35, 36, 45, 46, 47, 48], [49, 50, 51, 52, 61, 62, 63, 64, 73, 74, 75, 76, 85, 86, 87, 88], [53, 54, 55, 56, 65, 66, 67, 68, 77, 78, 79, 80, 89, 90, 91, 92], [57, 58, 59, 60, 69, 70, 71, 72, 81, 82, 83, 84, 93, 94, 95, 96], [97, 98, 99, 100, 109, 110, 111, 112, 121, 122, 123, 124, 133, 134, 135, 136], [101, 102, 103, 104, 113, 114, 115, 116, 125, 126, 127, 128, 137, 138, 139, 140], [105, 106, 107, 108, 117, 118, 119, 120, 129, 130, 131, 132, 141, 142, 143, 144], [145, 146, 147, 148, 157, 158, 159, 160, 169, 170, 171, 172, 181, 182, 183, 184], [149, 150, 151, 152, 161, 162, 163, 164, 173, 174, 175, 176, 185, 186, 187, 188], [153, 154, 155, 156, 165, 166, 167, 168, 177, 178, 179, 180, 189, 190, 191, 192], [193, 194, 195, 196, 205, 206, 207, 208, 217, 218, 219, 220, 229, 230, 231, 232], [197, 198, 199, 200, 209, 210, 211, 212, 221, 222, 223, 224, 233, 234, 235, 236], [201, 202, 203, 204, 213, 214, 215, 216, 225, 226, 227, 228, 237, 238, 239, 240], [241, 242, 243, 244, 253, 254, 255, 256, 265, 266, 267, 268, 277, 278, 279, 280], [245, 246, 247, 248, 257, 258, 259, 260, 269, 270, 271, 272, 281, 282, 283, 284], [249, 250, 251, 252, 261, 262, 263, 264, 273, 274, 275, 276, 285, 286, 287, 288]]


#Fill PMT module numbers.
for i in range(0,channels):
    pmt_mods.append(i+1)
    #print(pmt_mods[i])

class MyFirstGUI:
    LABEL_TEXT = [
        "This is our first GUI!",
        "Actually, this is our second GUI.",
        "We made it more interesting...",
        "...by making this label interactive.",
        "Go on, click on it again.",
    ]
    def __init__(self, primary):
        self.primary = primary
        primary.title("Hadron Calorimeter")
        primary.geometry("770x600")

        #Style().configure("TButton", padding=(0, 5, 0, 5),font='serif 10')

        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        canvas = Canvas(primary, borderwidth=0)

        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(primary, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

        for i in range(0,ncols):
            #primary.columnconfigure(i, pad=3)
            container_frame.columnconfigure(i, pad=3)

        for i in range(0,nrows+1):
            #primary.rowconfigure(i, pad=3)
            container_frame.rowconfigure(i, pad=3)

        #Fill a list with the actual buttons.
        #buttons = []
        for i in range(0,channels):
            for cluster in range(0,len(clusters)):
                for pmt in range(0,len(clusters[cluster])):
                    if clusters[cluster][pmt] == i+1:
                        if cluster%2==0:
                            bg_color = 'Gray66'
                        else:
                            bg_color = 'Gray88'
            self.button = Button(container_frame, text=str(pmt_mods[i])+'\n'+connections[str(i+1)][12], width=4, height=1,bg=bg_color)
            #button["command"] = self.fconnections
            self.button.bind("<Button-1>", self.fconnections)
            buttons.append(self.button)

        for i in range(0,nrows):
            for j in range(0,ncols):
                buttons[i*ncols+j].grid(row=i, column=j)
                #self.buttons[i*ncols+j].pack()

        #exit_btn = Button(primary, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.destroy, bg = "red") #Don't use destroy! It hangs the terminal for a long time!
        exit_btn = Button(container_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.grid(row=24, column=5, columnspan=2)

        #test_button = Button(primary, text='Test', command=self.fconnections)
        #self.test_button.pack()

        #self.label_index = 0
        #self.label_text = StringVar()
        #self.label_text.set(self.LABEL_TEXT[self.label_index])
        #self.label = Label(primary, textvariable=self.label_text)
        #self.label.bind("<Button-1>", self.cycle_label_text)
        #self.label.pack()

        #self.greet_button = Button(primary, text="Greet", command=self.greet)
        #self.greet_button.pack()

        #self.close_button = Button(primary, text="Close", command=primary.quit)
        #self.close_button.pack()

    def greet(self):
        print("Greetings!")

    def cycle_label_text(self, event):
        print(self,'  ',event)
        print(event.widget)
        self.label_index += 1
        self.label_index %= len(self.LABEL_TEXT) # wrap around
        self.label_text.set(self.LABEL_TEXT[self.label_index])

    def old_fconnections(self, event):
        pmt = event.widget.cget('text')
        print('PMT Module ',pmt,': amplifier channel = ', connections[pmt][0],', fADC channel = ', connections[pmt][1], ', HV channel = ',connections[pmt][2])
        #print(event.widget.cget('text'))

    def fconnections(self, event):
        pmt = event.widget.cget('text')
        #print(len(pmt))
        if pmt[2:3]== 'R':
            pmt = str(pmt[:1])
        elif pmt[2:3]== '\n':
            pmt = str(pmt[:2])
        else:
            pmt = str(pmt[:3])
        print('******************** Information for PMT Module '+pmt+' ('+connections[pmt][12]+') ********************')
        print('PMT module '+pmt+' ('+connections[pmt][12]+') is powered by HV channel '+str(connections[pmt][11])+'.')
        print('PMT module '+pmt+'\'s ('+connections[pmt][12]+') output goes to amplifier '+str(connections[pmt][0])+'.')
        print('This signal terminates at fADC '+str(connections[pmt][3])+', F1TDC '+str(connections[pmt][9])+', and summing module '+str(connections[pmt][10])+'.')
        print('The fADC data flow follows: amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
        print('The TDC data flow follows: amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
        print('The summing module data flow follows: amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')
        #print(event.widget.cget('text'))

root = Tk()
#print(root.call("info", "patchlevel")) #Prints tkinter version.
my_gui = MyFirstGUI(root)
root.mainloop()
