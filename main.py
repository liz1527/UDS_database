#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 16:15:14 2020

Code to initialise a GUI that will become the UDS database programme

@author: ppxee
"""

### Import Modules and define constants ###
import tkinter as tk
from PIL import Image, ImageTk # for resizing images
from astropy.table import Table
import numpy as np
colnum = 4

### Import catalogues ###
DR11cat = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019.fits')

def ID_submit():
    ''' This is a function that defines what happens when an ID is submitted
    '''
    ID = enter_txt.get()
    checkID = np.isin(ID, DR11cat['ID'])
    if checkID == False:
        disp_text = 'Invalid ID entered. Please retry'
        ID_txt.configure(text=disp_text) # changes text of the welcome label
    else:
        disp_ID_text = 'Object ID: '+enter_txt.get()
        ID_txt.configure(text=disp_ID_text) 
        ### get and display basic obdata ###
        obdata = DR11cat[np.isin(DR11cat['ID'], int(ID))]
        RA_text = 'RA: '+str(obdata['ALPHA_J2000'][0])
        RA_txt.configure(text=RA_text) 
        dec_text = 'Dec: '+str(obdata['DELTA_J2000'][0])
        dec_txt.configure(text=dec_text) 
        z_p_text = 'z_p = '+str(obdata['z_p'][0])
        z_p_txt.configure(text=z_p_text) 
        z_spec = obdata['z_spec'][0] # check if it has a spec z
        if z_spec == -1:
            z_spec_text = 'z_spec = Unknown'
            z_spec_txt.configure(text=z_spec_text) 
        else:
            z_spec_text = 'z_spec = '+str(z_spec)
            z_spec_txt.configure(text=z_spec_text) 

    return
    
#%% Create main window ###
window = tk.Tk() # define the main window of the GUI

window.geometry('900x700') # define the shape/size of the window

window.title("The UDS Database") # Give it a title that will appear in the bar

window.config(bg='skyblue')

#%% create title bar with logo ###
### get and resize logo ###
org_logo_image =Image.open('UDS_logo.png')
org_logo_image = org_logo_image.resize((370,60), Image.ANTIALIAS)
logo_image = ImageTk.PhotoImage(org_logo_image)

### Add logo to page ###
logo = tk.Label(window, image=logo_image) # Add the logo
logo.grid(column=1, row=0, columnspan=2, padx=1, pady=3)

### Create and add title ###
lbl = tk.Label(window, text="UDS Object Database", 
               font=("Arial Bold", 20), width=38, height=2, bg='skyblue') # Add a title
lbl.grid(column=0, row=0, padx=1, pady=3)

#%% create main page info ###
### Add row for entering ID ###
sublbl = tk.Label(window, text="Please enter the DR11 ID of your desired object", 
               font=("Arial", 15), bg='skyblue') # Add search text
sublbl.grid(column=0, row=2, padx=1, pady=3)

enter_txt = tk.Entry(window,width=10) # add a box where you can enter text
enter_txt.grid(column=1, row=2, padx=1, pady=3)

enter_btn = tk.Button(window, text="Enter", command=ID_submit) 
enter_btn.grid(column=2, row=2, padx=1, pady=3)

### Add empty labels for info to fill ###
def make_text_info(column, row, text, font=("Arial", 13), anchor=tk.W):
    txt = tk.Label(window, text=text, anchor=anchor, width=35,
                   font=font, bg='skyblue') # Add ID confirmation label
    txt.grid(column=column, row=row, padx=1, pady=3)
    return txt

placeholder_txt = make_text_info(column=0, row=3, text='', font=("Arial", 17)) # to creat gap
ID_txt = make_text_info(column=0, row=4, text='', font=("Arial", 17),
                        anchor=tk.CENTER) # where the ID will be shown
RA_txt = make_text_info(column=0, row=5, text='',) # where RA will be shown
dec_txt = make_text_info(column=0, row=6, text='',)#  where dec will be shown
z_p_txt = make_text_info(column=0, row=7, text='') # where z will be shown
z_spec_txt = make_text_info(column=0, row=8, text='') # where z_spec will be shown


window.mainloop() 