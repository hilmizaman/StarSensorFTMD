from tkinter import *
from PIL import ImageTk,Image
import cv2
import nested_function as nf
from functools import partial
import numpy as np
from contextlib import suppress


def create_star_image_after_button_pressed(ra,de,roll,f,myu,canvas):
    """[create star image command when button is pressed]

    Args:
        ra ([int]): [right ascension from the spinner]
        de ([int]): [declination from the spinner]
        roll ([int]): [roll from the spinner]
        f ([float]): [focal length from the entry in mm]
        myu ([float]): [length per pixel from the entry in nm]
        canvas ([tkinter canvas]): [the canvas to draw the image in]
    """
    ra_calc = (ra.get())
    de_calc = (de.get())
    roll_calc = (roll.get())
    f_calc = float((f.get()))/1000
    myu_calc = float((myu.get()))*(10**-6)
    star_image = nf.create_star_image(ra_calc,de_calc,roll_calc,f_calc,myu_calc).astype(np.uint8)
    dsize = (int(star_image.shape[1] * 25/100),int(star_image.shape[0] * 25/100))
    rescaled_image = cv2.resize(star_image,dsize)
    tkinter_image = ImageTk.PhotoImage(image=Image.fromarray(rescaled_image))
    image_id = canvas.create_image(20,20,anchor=NW,image=tkinter_image)
    canvas.configure(image=tkinter_image)    
    return


mainWindow = Tk()
mainWindow.title("Star Simulator")
mainWindow.geometry("1000x780+8+8")
mainWindow.configure(bg='skyblue')

created_label = Label(mainWindow,bg='steelblue',text="Created by Brian Mohammed Catraguna, Astronautics Laboratory, Faculty of Mechanical and Aerospace Engineering, Bandung Institue of Technology")
created_label.grid(row=2,column=0,sticky='nsew',columnspan=3)
created_label.config(font=("Courier", 7))

mainWindow.columnconfigure(0,weight=1)

#INPUT FRAME
inputframe = Frame(mainWindow,width=1000)
inputframe.grid(row=1,column=0,sticky='n')
inputframe.config(relief='sunken',borderwidth=3)

#Attitude Sub - Frame
attitudeframe = LabelFrame(inputframe,text="Attitude")
attitudeframe.config(font=("system",9))
attitudeframe.grid(row=0,column=0)

ra = IntVar()
de = IntVar()
roll = IntVar()

f = StringVar()
myu = StringVar()

#Setting default values
f.set(3.04)
myu.set(1.12)

ra_input_label = Label(attitudeframe,text="Right Ascension α (degrees): ")
ra_input_label.grid(row=0,column=0)
ra_input = Spinbox(attitudeframe,width=10,values=tuple(range(-180,180)),textvariable=ra)
ra_input.grid(row=0,column=1)

de_input_label = Label(attitudeframe,text="Declination δ (degrees): ")
de_input_label.grid(row=1,column=0)
de_input = Spinbox(attitudeframe,width=10,values=tuple(range(-90,90)),textvariable=de)
de_input.grid(row=1,column=1)

roll_input_label = Label(attitudeframe,text="Roll φ (degrees): ")
roll_input_label.grid(row=2,column=0)
roll_input = Spinbox(attitudeframe,width=10,values=tuple(range(0,360)),textvariable=roll)
roll_input.grid(row=2,column=1)

#Sensor Settings Sub - Frame
settingsframe = LabelFrame(inputframe,text="Sensor Settings")
settingsframe.config(font=("system",9))
settingsframe.grid(row=0,column=1)

focal_length_label = Label(settingsframe,text="Focal Length f (mm): ")
focal_length_label.grid(row=0,column=0)
focal_length = Entry(settingsframe,textvariable=f)
focal_length.grid(row=0,column=1)

miu_label = Label(settingsframe,text="Length per Pixel μ (μm): ")
miu_label.grid(row=1,column=0)
miu = Entry(settingsframe,textvariable=myu)
miu.grid(row=1,column=1)


#OUTPUT FRAME
outputframe = Frame(mainWindow,width=1000)
outputframe.grid(row=0,column=0,sticky='s')
outputframe.config(relief='ridge',borderwidth=3)

#Creating Canvas for Showing Image
canvas = Canvas(outputframe,width=850,height=664)
canvas.grid(row=0,column=0,sticky='nsew')

#Generate Star Image Button
create_star_image_after_button_pressed = partial(create_star_image_after_button_pressed,ra,de,roll,f,myu,canvas)
generate_button = Button(inputframe,text="Generate Star Image",command=create_star_image_after_button_pressed)
generate_button.grid(row=0,column=2)

mainWindow.mainloop()