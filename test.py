#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 16:15:14 2020

Code to test initialising a GUI 
Will contain explanations of some functions

@author: ppxee
"""

import tkinter as tk
### Import requirements for embedding matplotlib ###
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
from PIL import Image, ImageTk # for resizing images

def clicked():
    ''' This is a function that defines what happens when a button is clicked
    '''
    lbl.configure(text="Button was clicked !!"+ txt.get()) # changes text of the welcome label
    
    ### test image ###
    ID = txt.get()
    IDstr = str(ID)
    ob_image =Image.open('all_obs/'+IDstr+'_all_band_images.png')
    ob_image = ob_image.resize((400,400), Image.ANTIALIAS)
    ob_image = ImageTk.PhotoImage(ob_image)
    img = tk.Label(window, image=ob_image) # Add the logo
    img.image = ob_image 
    img.grid(column=0, row=4, padx=1, pady=3)
    
window = tk.Tk() # define the main window of the GUI

window.geometry('800x800') # define the shape/size of the window

window.title("The UDS Database") # Give it a title that will appear in the bar

lbl = tk.Label(window, text="Welcome", font=("Arial Bold", 50)) # label adds text to the window

''' grid itself allows you to define exactly where the writing will appear'''
lbl.grid(column=0, row=0) # the label will not sgow until something like grid is called

''' add a box where you can enter text '''
txt = tk.Entry(window,width=10) # 

txt.grid(column=1, row=0)

''' you can add a button and change its colours as shown'''
btn = tk.Button(window, text="Click Me", bg="orange", fg="red", command=clicked) 

btn.grid(column=2, row=0)

### test plot ###
fig = Figure(figsize=(5, 4), dpi=100)
t = np.arange(0, 3, .01)
fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))

canvas = FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().grid(column=0, row=3)

#### test image ###
#IDstr = str(3)
#ob_image =Image.open('all_obs/'+IDstr+'_all_band_images.png')
#ob_image = ob_image.resize((400,400), Image.ANTIALIAS)
#ob_image = ImageTk.PhotoImage(ob_image)
#img = tk.Label(window, image=ob_image) # Add the logo
#img.grid(column=0, row=4, padx=1, pady=3)

''' window.mainloop() tells Python to run the Tkinter event loop. This method 
listens for events, such as button clicks or keypresses, and blocks any code 
that comes after it from running until the window it’s called on is closed. '''
window.mainloop() 