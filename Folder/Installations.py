#!/usr/bin/env python
# coding: utf-8

# In[63]:

print('PREPARING THE INSTALLATIONS\n')
import tkinter as tk
from tkinter import * 
from tkinter import filedialog
import requests
import urllib.request
from urllib.request import urlopen
from tkinter import messagebox
import subprocess
import matplotlib as mpl
from matplotlib import cm
import os
import webbrowser
import tkinter
import tkinter.colorchooser


urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/icon_nevomics.ico',
                           'icon_nevomics.ico')

root = Tk()
root.attributes("-topmost", True)
root.title("NeVOmics")
#           ancho , alto
root.geometry("440x460")
root.iconbitmap(r'icon_nevomics.ico')


lab0 = Label(root, text="       ")# columna vacía
lab0.grid(column=0, row=0)

lab1 = Label(root, text="NeVOmics",
             font=("Arial", 20, "bold"),
             fg = 'green')
lab1.grid(column=1, columnspan=2, row=0, sticky= W)

lab11 = Label(root, text="Installations",
             font=("Arial", 15, "bold"),
             fg = 'grey')
lab11.grid(column=2, columnspan=3, row=0, sticky= W)


#### python

label = Label(root, text="If you don't have Python Modules installed yet",
              font=("Arial", 12, "bold"))
label.grid(column=1,columnspan = 6, row=1, sticky= S+W)

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/pymod10.png',
                           'pymod10.png')

image0 = tk.PhotoImage(file= 'pymod10.png')

label1 = Label(root, text="1. Click to install:",
             font=("Arial", 10))
label1.grid(column=1,columnspan = 6, row=2, sticky= S+W)

def inspymod():
    modulo = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/install_python_modules.py',
                                        'install_python_modules.py')
    os.system('start cmd /c python install_python_modules.py ^&^& del install_python_modules.py pymod10.png')
            
botonrun = tkinter.Button(root, image=image0,
                          command = inspymod, bd = 0,
                          activebackground= 'yellow',activeforeground= 'black', cursor="hand2")
botonrun.grid(column = 2, columnspan = 3, row = 2, rowspan = 2, sticky= W)


lab0 = Label(root, text="       ")# columna vacía
lab0.grid(column=0, row=4)

###################

label = Label(root, text="If you don't have R installed yet",
              font=("Arial", 12, "bold"))
label.grid(column=1,columnspan = 6, row=5, sticky= S+W)

label1 = Label(root, text="1. Click to download:",
             font=("Arial", 10))
label1.grid(column=1,columnspan = 6, row=6, sticky= S+W)

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Rlogo25.png',
                           'Rlogo25.png')

image1 = tk.PhotoImage(file= 'Rlogo25.png')

def callback(event):
    webbrowser.open_new(r"https://cran.r-project.org/bin/windows/base/old/3.5.3/R-3.5.3-win.exe")
link = Label(root, image=image1,
             font=("Arial", 10), fg="blue", cursor="hand2")
link.grid(column=2, row=6, sticky= W, columnspan=1)
link.bind("<Button-1>", callback)

label2 = Label(root, text="2. Search R-3.5.3-win.exe to install. ",
              font=("Arial", 10))
label2.grid(column=1,columnspan = 6, row=7, sticky= W)


lab0 = Label(root, text="       ")# columna vacía
lab0.grid(column=0, row=8)

label = Label(root, text="Once installed R",
              font=("Arial", 12, "bold"))
label.grid(column=1,columnspan = 6, row=9, sticky= S+W)

label3 = Label(root, text="If you don't have R library yet",
              font=("Arial", 10, "bold"))
label3.grid(column=1,columnspan = 6, row=10, sticky= S+W)

def parameters():
    urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Rlibrary_NeVOmics.py',
                               'Rlibrary_NeVOmics.py')
    os.system('start cmd /k python Rlibrary_NeVOmics.py ^&^& del Rlibrary_NeVOmics.py')

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/box1.png',
                           'box1.png')

    
image2 = tk.PhotoImage(file= 'box1.png')

label4 = Label(root, text="1. Click to download:",
              font=("Arial", 10))
label4.grid(column=1,row=11, sticky= S+W)

botonrun2 = tkinter.Button(root, image=image2,
                          command = parameters, bd = 0,
                          activebackground= 'yellow',activeforeground= 'black', cursor="hand2")
botonrun2.grid(column = 2, columnspan = 6, row = 11, rowspan = 2, sticky= W)

lab1 = Label(root, text="       ")# columna vacía
lab1.grid(column=0, row=14)

urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/NCBI.png',
                           'NCBI.png')

label33 = Label(root, text="If you don't have BLAST installed yet",
              font=("Arial", 12, "bold"))
label33.grid(column=1,columnspan = 6, row=15, sticky= S+W)


image3 = tk.PhotoImage(file= 'NCBI.png')

def callback1(event):
    webbrowser.open_new(r"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-win64.exe")
link2 = Label(root, image=image3,
             font=("Arial", 10), fg="blue", cursor="hand2")
link2.grid(column=2, row=16, sticky= W, columnspan=1)
link2.bind("<Button-1>", callback1)

label3 = Label(root, text="1. Click to download:",
             font=("Arial", 10))
label3.grid(column=1,columnspan = 13, row=16, sticky= S+W)

label4 = Label(root, text="2. Search ncbi-blast-2.8.1+-win64.exe to install. ",
              font=("Arial", 10))
label4.grid(column=1,columnspan = 6, row=17, sticky= W)



root.mainloop()

if os.path.exists('Rlogo25.png'): os.remove('Rlogo25.png')
if os.path.exists('icon_nevomics.ico'): os.remove('icon_nevomics.ico')
if os.path.exists('box1.png'): os.remove('box1.png')
if os.path.exists('NCBI.png'): os.remove('NCBI.png')
if os.path.exists('pymod10.png'): os.remove('pymod10.png')