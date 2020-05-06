#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 16:15:14 2020

Code to initialise a GUI that will become the UDS database programme

@author: ppxee
"""

#%% Import Modules and define constants ###
import tkinter as tk
from PIL import Image, ImageTk # for resizing images
from astropy.table import Table # for data handling
import numpy as np # for array handling
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg # to plot figures
from matplotlib.figure import Figure# to create figures
from matplotlib.patches import Rectangle # to create field plot
from os import path # to check if a file exists
from astropy import units as u
colnum = 4

### Import catalogues ###
DR11cat = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019_database.fits')

#%% define functions ### 
def ID_submit():
    ''' This is a function that defines what happens when an ID is submitted
    '''
    hide_details() # remove details currently entered (if any)
    ID = enter_txt.get()
    checkID = np.isin(ID, DR11cat['ID'])
    if checkID == False:
        disp_text = 'Invalid ID entered. Please retry'
        ID_txt.configure(text=disp_text) # changes text of the welcome label
        hide_buttons()
    else:
        obdata = DR11cat[np.isin(DR11cat['ID'], int(ID))]
        info_btn.grid(column=0, row=3, padx=1, pady=3)
        img_btn.grid(column=1, row=3, padx=1, pady=3)
        vary_btn.grid(column=2, row=3, padx=1, pady=3, columnspan=2)
        hide_images()
        hide_varys()
        show_basic_details(obdata)
        create_field_plot(obdata)
    return

def show_basic_details(obdata):
    ''' This function displays the basic data about an object from its DR11
    catalogue entry
    Inputs: 
        obdata = row of DR11cat that refers to the chosen ID
    '''
    ID = enter_txt.get()
    ### change button colours ###
    info_btn['fg'] = "red"
    img_btn['fg'] = "black"
    vary_btn['fg'] = "black"
    
    ### get and display basic obdata ###
    disp_ID_text = 'Object ID: '+str(ID)
    ID_txt.configure(text=disp_ID_text) 
    RA_text = 'RA: '+str(obdata['ALPHA_J2000'][0])
    RA_txt.configure(text=RA_text) 
    dec_text = 'Dec: '+str(obdata['DELTA_J2000'][0])
    dec_txt.configure(text=dec_text) 
    z_p_text = 'z_p = '+str(obdata['z_p'][0])
    z_p_txt.configure(text=z_p_text) 
    z_spec = obdata['z_spec'][0] # check if it has a spec z
    if z_spec == -1:
        z_spec_text = 'z_spec = Not Known'
    else:
        z_spec_text = 'z_spec = '+str(z_spec)
        ### find if ob has an available spectrum ###
        spec_data = Table.read('Spectra/UDS_spectral_measurements_DR2_v4.fits')
        spec_data = spec_data[spec_data['DR11-ID']==int(ID)]
        if len(spec_data) == 1: # check if found match
            create_spec_plot(spec_data)
        else:
            ax = specfig.add_subplot(111)
            ax.text(0.4,0.5, 'No Spectrum Available')
            speccanvas = FigureCanvasTkAgg(specfig, master=window)  # A tk.DrawingArea.
            speccanvas.draw()
            speccanvas.get_tk_widget().grid(column=0, row=20, columnspan=colnum, rowspan=1)
            
    z_spec_txt.configure(text=z_spec_text) 
    mag = obdata['KMAG_20'][0]
    mag_text = 'K-band mag: '+str(round(mag, 2))
    mag_txt.configure(text=mag_text) 
    mass = obdata['Mstar_z_p'][0]
    mass = format(mass,'5.2e')
    mass_text = 'Stellar mass: '+mass+' '+'M_solar'
    mass_txt.configure(text=mass_text) 
    xray_text = 'X-ray detection: '+str(obdata['X-ray'][0])
    xray_txt.configure(text=xray_text) 
        
    ### display the sample the object is in ###'
    if obdata['Stars-final'] == True and obdata['Best'] == True:
        samp_text = 'Best Star'
    elif obdata['Stars-final'] == True and obdata['Good'] == True:
        samp_text = 'Good Star'
    elif obdata['Stars-final'] == True:
        samp_text = 'Star'
    elif obdata['Best galaxies'] == True:
        samp_text = 'Best Galaxy'
    elif obdata['Good galaxies'] == True:
        samp_text = 'Good Galaxy'
    elif obdata['XTALK-DR11'] == True:
        samp_text = 'Cross Talk Contaminated'
    else:
        samp_text = 'Galaxy'
    samp_txt.configure(text='Designation: '+samp_text) 

def create_field_plot(obdata):
    ''' Function that creates a plot of the UDS field, and places a circle 
    where the object with the entered ID can be found.
    Inputs:
        obdata = row of the DR11 catalogue that relates to the ID entered.
    '''
    ### set up figure and axes ###
    fig = Figure(figsize=(4, 4), dpi=100)
    fig.patch.set_facecolor('skyblue')

    ax = fig.add_subplot(111)
#    fig.patch.set_facecolor('xkcd:aqua blue')
#    ax.set_facecolor('xkcd:aqua blue')
    
    # arrays set out as min RA, max RA, min Dec, max Dec
    UDS = np.array([33.98, 34.92, -5.66, -4.64])
    Chandra = np.array([34.07,34.72, -5.40, -4.93])
    
    ### Create and add rectangles ###
    udsrect = Rectangle([UDS[0],UDS[2]], 
                        width = UDS[1]-UDS[0], height = UDS[3]-UDS[2],
                        linewidth=2,edgecolor='b',facecolor='none')
    chanrect = Rectangle([Chandra[0],Chandra[2]],
                         width = Chandra[1]-Chandra[0], height = Chandra[3]-Chandra[2],
                         linewidth=2,edgecolor='r',facecolor='none', linestyle='--')
      
    ax.add_patch(udsrect)
    ax.add_patch(chanrect)
    
    ### Add point where object is ###
    ax.plot(obdata['ALPHA_J2000'], obdata['DELTA_J2000'], 'ko')
    
    ### adjust plot settings ###
    ax.set_xlim(xmin=33.8, xmax=35.2)
    ax.set_ylim(ymin=-5.70, ymax=-4.35)
    ax.axis('equal')
    ax.text(UDS[1]-0.01, UDS[2]+0.01, 'UDS', color='b')
    ax.text(Chandra[1]-0.01, Chandra[2]+0.01, 'Chandra', color='r')
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.invert_xaxis()# flip axes for RA
    fig.tight_layout()
    
    ### Add plot to GUI ###
    canvas = FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().grid(column=1, row=5, columnspan=3, rowspan=15)
    
def create_spec_plot(specdata):
    ''' Function that creates a plot of the spectrum of the sources, if one is
    available.
    Inputs: 
        specdata = a 1 row table which details at least the redshift and spec
                    ID of the object so that the spectrum can be imported and
                    plotted.
    '''
    ### set up figure and axes ###
    specfig = Figure(figsize=(9.2, 2), dpi=100)
    specfig.patch.set_facecolor('skyblue')
    ax = specfig.add_subplot(111)
    
    ### set up variables to plot spectra
    linedict = {'CIV': 1549.00,
                'SiII': 1526.70,
                'FeII': 1608.50,
                'CIII': 1909.00,
                'MgII': 2799.00,
                'OII': 3727.50,
                'OIII_doublet-1': 5006.80,
                'OIII_doublet-1/3': 4958.90,
                'MgI': 5175.40,
                'LyA': 1215.70,
                'Hbeta': 4861.30,
                'Hgamma': 4340.40,
                'Hdelta': 4101.70,
                'Balmer_break': 4000,
                'NV': 1240.00,
                'SiIV+OIV': 1397.00,
                'LyB': 1025.60,
                'OVI': 1033.01}
    cat_ID = specdata['cat-ID'][0]
    z = specdata['z-spec'][0]
    spec = Table.read('Spectra/UDS-spectra-Feb-19/UDS-'+str(cat_ID)+'-1Dspec.fits')
    
    ax.plot(spec['wavelength'], spec['flux'],'b')
    
    ylims = ax.get_ylim()
#    xlims = axes.get_xlim()
    
    for key in linedict:
        val = linedict[key]
        lineval = val*(z+1)
        if lineval < max(spec['wavelength']) and lineval > min(spec['wavelength']):
            ax.vlines(lineval, ymin=ylims[0], ymax=ylims[1], linestyles='dashed')
            if key == 'SiII' or key =='LyA' or key == 'LyB':
                ax.text(lineval-220, ylims[0]+0.5e-18, key)
            elif key == 'SiIV_1393+OIV1402':
                ax.text(lineval-400, ylims[0]+0.2e-18, key)
            elif key == 'Balmer_break':
                ax.text(lineval-300, ylims[1]+0.2e-18, key)
            elif key == 'OIII_doublet-1/3':
                ax.text(lineval-200, ylims[1]+0.2e-18, 'OIII_doublet')
            elif key == 'OIII_doublet-1':
                continue
            else:
                ax.text(lineval+10, ylims[0]+0.2e-18, key)
        
    ax.set_ylim(ylims)
    ax.set_ylabel('Flux')
    unit = u.AA
    ax.set_xlabel('Observed Wavelength, '+unit.to_string('latex'))
    
    ### Add plot to GUI ###
    speccanvas = FigureCanvasTkAgg(specfig, master=window)  # A tk.DrawingArea.
    speccanvas.draw()
    speccanvas.get_tk_widget().grid(column=0, row=20, columnspan=colnum, rowspan=1)
    
def hide_details():
    ''' This function hide object basic details '''
    ### hide basic obdata ###
    text = ''
    RA_txt.configure(text=text) 
    dec_txt.configure(text=text) 
    z_p_txt.configure(text=text) 
    z_spec_txt.configure(text=text) 
    mag_txt.configure(text=text) 
    mass_txt.configure(text=text) 
    xray_txt.configure(text=text) 
    samp_txt.configure(text=text) 
    
    ### hide plots ###
    fig = Figure(figsize=(4, 4), dpi=100)
    fig.patch.set_facecolor('skyblue')
    canvas= FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().grid(column=1, row=5, columnspan=3, rowspan=15)
    
    specfig = Figure(figsize=(9.2, 2), dpi=100)
    specfig.patch.set_facecolor('skyblue')
    speccanvas = FigureCanvasTkAgg(specfig, master=window)  # A tk.DrawingArea.
    speccanvas.draw()
    speccanvas.get_tk_widget().grid(column=0, row=20, columnspan=colnum, rowspan=1)
    
def show_images():
    ''' Code to show the images of an object when the Images button is pressed
    '''
    ### change button colours ###
    info_btn['fg'] = "black"
    img_btn['fg'] = "red"
    vary_btn['fg'] = "black"
    
    hide_details()
    hide_varys()
    ID = enter_txt.get()
    ### get and resize image ###
    ID = str(ID)
    filename = 'all_obs/'+ID+'_all_band_images.png'
    if path.exists(filename):
        ob_image =Image.open(filename)
        ob_image = ob_image.resize((400,400), Image.ANTIALIAS)
        tk_ob_image = ImageTk.PhotoImage(ob_image)
        ### Add img to page ###
        img = tk.Label(window, image=tk_ob_image)
        img.image = tk_ob_image
    else: 
        ### set to dummy text ###
        img = tk.Label(window, text='Image not yet available', font=('Arial', 20),
                       bg="skyblue")
        
    img.grid(column=0, row=7, columnspan=4, padx=1, pady=3)

def hide_images():
    ### get slaves ###
    widges = window.grid_slaves(row=7, column=0)
#    widges.grid_forget()
    for wid in widges:
        if wid.grid_info()['columnspan'] == 4: # identifies the image rather than text label
            wid.grid_forget()

def show_varys():
    ''' Function to show the variability data on the chosen object
    This uses the labels named outside the loop as need same positions - sorry 
    that this makes the variable names less clear
    '''
    ### hide any currently displayed data
    hide_details()
    hide_images()
    
    ### change button colours ###
    info_btn['fg'] = "black"
    img_btn['fg'] = "black"
    vary_btn['fg'] = "red"
    
    ### Get variability data and ID ###
    vary_cat = Table.read('UDS_catalogues/full_varystats_noneg.fits')
    ID = enter_txt.get()
    
    ### Reduce to relevant row/find if has data ###
    vary_cat = vary_cat[vary_cat['ID']==int(ID)]
    if len(vary_cat) == 1:
        ### Determine if object is variable and disp text ###
        if vary_cat['Chi_K'] > 30 or vary_cat['Chi_J'] > 32:
            vary_txt = 'Identified as variable'
        else:
            vary_txt = 'Not identified as variable'
        samp_txt.configure(text=vary_txt)
        
        ### Display stats ###
        kchi = str(round(vary_cat['Chi_K'][0],2))
        RA_txt.configure(text='Chi^2 in K = '+kchi) 
        ksig = str(round(vary_cat['sig_K'][0],2))+' +/- '+str(round(vary_cat['sig_K_err'][0],2))
        dec_txt.configure(text='sigma in K = '+ksig) 
        
        jchi = str(round(vary_cat['Chi_J'][0],2))
        z_p_txt.configure(text='Chi^2 in J = '+jchi) 
        jsig = str(round(vary_cat['sig_J'][0],2))+' +/- '+str(round(vary_cat['sig_J_err'][0],2))
        z_spec_txt.configure(text='sigma in J = '+jsig) 

        ### show light curves on fig ###
        fig = Figure(figsize=(4, 4), dpi=100)
        fig.patch.set_facecolor('skyblue')
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        
        ### set up time variable for plot ###
        t = np.linspace(1, 8, num=8)
        years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
        x = [1,3,4,5,6,7,8]
        
        ### set up flux variables ###
        kflux = vary_cat['Flux_K'][0,:]
        kfluxerr = vary_cat['Fluxerr_K'][0,:]
        
        jflux = vary_cat['Flux_J'][0,:]
        jfluxerr = vary_cat['Fluxerr_J'][0,:]
#        
#        ### Convert to magnitudes ###
#        kmag = 30 - 2.5*np.log10(kflux)
#        kmag += 1.9 # to get AB mag
#        kmagerr = 1.086/(kflux/kfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE
#
#        jmag = 30 - 2.5*np.log10(jflux)
#        jmag += 1.9 # to get AB mag
#        jmagerr = 1.086/(jflux/jfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

        ### plot light curves ###
        ax1.errorbar(x, kflux, yerr=kfluxerr,fmt='ro')
#        ax1.errorbar(x, kmag, yerr=kmagerr,fmt='ro')
        ax1.set_xticks(t)
        ax1.set_xticklabels(years)
        ax1.set_xlabel('Semester')
        ax1.set_ylabel('$K$-Band Flux')
#        ax1.set_ylabel('$K$-Band Magnitude')
#        ax1.invert_yaxis()
        
        ax2.errorbar(t, jflux, yerr=jfluxerr,fmt='bo')
#        ax2.errorbar(t, jmag, yerr=jmagerr,fmt='bo')
        ax2.set_xticks(t)
        ax2.set_xticklabels(years)
        ax2.set_xlabel('Semester')
        ax2.set_ylabel('$J$-Band Flux')
#        ax2.set_ylabel('$J$-Band Magnitude')
#        ax2.invert_yaxis()
        
        fig.tight_layout()
        
        ### Add to tkinter window ###
        canvas= FigureCanvasTkAgg(fig, master=window)  
        canvas.draw()
        canvas.get_tk_widget().grid(column=1, row=5, columnspan=3, rowspan=15)
        
    else:
        ### Show that object is not in catalogue ###
        vary_txt = 'Does not meet criteria for variability selection'
        vary = tk.Label(window, text=vary_txt, font=('Arial', 20),
                           bg="skyblue")
        vary.grid(column=0, row=7, columnspan=colnum, padx=1, pady=3)
        
def hide_varys():
    ### get slaves ###
    widges = window.grid_slaves(row=7) # get buttons details
#    widges.grid_forget()
    for wid in widges:
        if wid.grid_info()['columnspan'] == colnum: 
            wid.grid_forget()
    ### hide figure ###
    fig = Figure(figsize=(4, 4), dpi=100)
    fig.patch.set_facecolor('skyblue')
    canvas= FigureCanvasTkAgg(fig, master=window)  
    canvas.draw()
    canvas.get_tk_widget().grid(column=1, row=5, columnspan=3, rowspan=15)
    
def hide_buttons():
    ### get slaves ###
    widges = window.grid_slaves(row=3) # get buttons details
#    widges.grid_forget()
    for wid in widges:
        if wid.grid_info()['column'] == 0 or wid.grid_info()['column'] == 1 or wid.grid_info()['column'] == 2: 
            wid.grid_forget()

def return_basic_details():
    ''' This function hides the images so the basic details can be reshown,
    then returns to the ID_submit function so that the details plot'''
    hide_images()
    ID_submit()    
    
def make_text_info(column, row, text, font=("Arial", 13), anchor=tk.W):
    ''' This function creats a standard label for info in the basic details 
    Inputs:
        column = which column of the grid to create the label on
        row = which row of the grid to create the label on
        text = text string to be included on the label
        font = font details for the label, default is Arial in size 13
        anchor = alignment of the text within the label, default is tk.W which
                aligns text to the left of the box
    Outputs:
        txt = a Label object with these specifications
        '''
    txt = tk.Label(window, text=text, anchor=anchor, width=35,
                   font=font, bg='skyblue') # Add ID confirmation label
    txt.grid(column=column, row=row, padx=1, pady=3)
    return txt
#%% Create main window ###
window = tk.Tk() # define the main window of the GUI

window.geometry('950x800') # define the shape/size of the window

window.title("The UDS Database") # Give it a title that will appear in the bar

window.config(bg='skyblue')

#%% create title bar with logo ###
### get and resize logo ###
org_logo_image =Image.open('UDS_logo.png')
org_logo_image = org_logo_image.resize((370,60), Image.ANTIALIAS)
logo_image = ImageTk.PhotoImage(org_logo_image)

### Add logo to page ###
logo = tk.Label(window, image=logo_image) # Add the logo
logo.grid(column=2, row=0, columnspan=2, padx=1, pady=3)

### Create and add title ###
lbl = tk.Label(window, text="UDS Object Database", 
               font=("Arial Bold", 20), width=38, height=2, bg='skyblue') # Add a title
lbl.grid(column=0, row=0, padx=1, pady=3, columnspan=2)

#%% create main page info ###
### Add row for entering ID ###
sublbl = tk.Label(window, text="Please enter the DR11 ID of your desired object", 
               font=("Arial", 15), bg='skyblue') # Add search text
sublbl.grid(column=0, row=2, padx=1, pady=3, columnspan=2)

enter_txt = tk.Entry(window,width=10) # add a box where you can enter text
enter_txt.grid(column=2, row=2, padx=1, pady=3)

enter_btn = tk.Button(window, text="Enter", command=ID_submit) 
enter_btn.grid(column=3, row=2, padx=1, pady=3)

info_btn = tk.Button(window, text="Basic Info", command=return_basic_details,
                     bg="skyblue") 

img_btn = tk.Button(window, text="Images", command=show_images, bg="skyblue") 

vary_btn = tk.Button(window, text="Variability", command=show_varys, bg="skyblue") 
### Add empty labels for info to fill ###

placeholder_txt = make_text_info(column=0, row=4, text=' ') # to creat gap
ID_txt = make_text_info(column=0, row=5, text='', font=("Arial", 17),
                        anchor=tk.CENTER) # where the ID will be shown
samp_txt = make_text_info(column=0, row=6, text='', font=("Arial", 17),
                        anchor=tk.CENTER) # where the sample the object is in is shown
RA_txt = make_text_info(column=0, row=7, text='') # where RA will be shown
dec_txt = make_text_info(column=0, row=8, text='')#  where dec will be shown
z_p_txt = make_text_info(column=0, row=9, text='') # where z will be shown
z_spec_txt = make_text_info(column=0, row=10, text='') # where z_spec will be shown
mag_txt = make_text_info(column=0, row=11, text='') # where z_spec will be shown
mass_txt = make_text_info(column=0, row=12, text='') # where z_spec will be shown
xray_txt = make_text_info(column=0, row=13, text='') # where z_spec will be shown

### Add plot to GUI ###
fig = Figure(figsize=(4, 4), dpi=100)
fig.patch.set_facecolor('skyblue')
canvas = FigureCanvasTkAgg(fig, master=window)  # A tk.DrawingArea.
canvas.draw()
canvas.get_tk_widget().grid(column=1, row=5, columnspan=3, rowspan=15)

specfig = Figure(figsize=(9.2, 2), dpi=100)
specfig.patch.set_facecolor('skyblue')
speccanvas = FigureCanvasTkAgg(specfig, master=window)  # A tk.DrawingArea.
speccanvas.draw()
speccanvas.get_tk_widget().grid(column=0, row=20, columnspan=colnum, rowspan=1)

### Create label for image ###
img = tk.Label(window) # Add the logo

window.mainloop() 