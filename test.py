#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 16:15:14 2020

Code to test initialising a GUI 
Will contain explanations of some functions

@author: ppxee
"""

import tkinter as tk

def clicked():
    ''' This is a function that defines what happens when a button is clicked
    '''
    lbl.configure(text="Button was clicked !!"+ txt.get()) # changes text of the welcome label
    
window = tk.Tk() # define the main window of the GUI

window.geometry('350x200') # define the shape/size of the window

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



''' window.mainloop() tells Python to run the Tkinter event loop. This method 
listens for events, such as button clicks or keypresses, and blocks any code 
that comes after it from running until the window it’s called on is closed. '''
window.mainloop() 