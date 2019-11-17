#!/usr/bin/env python3

print('\nSTARTING NEVOMICS\n')
print('You can minimize this window.\n')
import tkinter
import tkinter as tk
from tkinter import * 
from tkinter import filedialog
import webbrowser
import pandas
import pandas as pd
from io import StringIO # cambio de "from pandas.compat import StringIO" a  "from io import StringIO" para compatibilidad con pandas 0.25.0
import requests
import urllib.request
from urllib.request import urlopen
import base64
import tkinter
from tkinter import ttk
import numpy as np
from tkinter import messagebox
import subprocess


from colormap import Colormap
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import *
import matplotlib as mpl
from matplotlib import cm
import os
import shutil
import numpy as np
import webbrowser
import tkinter
import tkinter.colorchooser
import matplotlib

valores = []
n = 0.005
while n < 0.205:
    valores.append(str(np.round(n * 100, 1)))
    n += 0.005
valores = ['0.05', '0.1', '0.2', '0.3', '0.4'] + valores

Pvalues_valores = []
n = 0.005
while n < 0.205:
    Pvalues_valores.append(str(np.round(n, 10)))
    n += 0.005
Pvalues_valores = ['0.0005', '0.001' , '0.002' , '0.003' , '0.004'] + Pvalues_valores

kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()
kegg_organism=pd.read_csv(StringIO(kegg_orgs),names=['T_number','Prefix','Organism','Group'],sep='\t')
kegg_organism = kegg_organism.sort_values(by ='Organism',ascending=True)


dict_org = {}
for index, row in kegg_organism.iterrows():
    dict_org[row.Organism] = [row.Prefix, row.T_number]

sequentials_colors = {'YlOrRd':cm.YlOrRd,'YlOrBr':cm.YlOrBr,'YlGnBu':cm.YlGnBu,
                      'YlGn':cm.YlGn,'Reds':cm.Reds,'PdPu':cm.RdPu,'Purples':cm.Purples,'PuRd':cm.PuRd,
                      'PuBuGn':cm.PuBuGn,'PuBu':cm.PuBu,'Oranges':cm.Oranges,'OrRd':cm.OrRd,
                      'Greys':cm.Greys,'Greens':cm.Greens,'GnBu':cm.GnBu,'BuPu':cm.BuPu,
                      'BuGn':cm.BuGn,'Blues':cm.Blues}

diverging_colors = {'GreBlaRed':Colormap().cmap_linear('green', 'black', 'red'),
                    'RedBlaGre':Colormap().cmap_linear('red', 'black', 'green'),
                    'RedBlaBlu':Colormap().cmap_linear('red', 'black', 'blue'),
                    'BluBlaRed':Colormap().cmap_linear('blue', 'black', 'red'),
                    'RedYlBlu':Colormap().cmap_linear('red', 'yellow', 'blue'),
                    'BluYlRed':Colormap().cmap_linear('blue', 'yellow', 'red'),
                      'seismic':cm.seismic,'coolwarm':cm.coolwarm,'bwr':cm.bwr,
                      'Spectral':cm.Spectral,'RdYlGn':cm.RdYlGn,'RdYlBu':cm.RdYlBu,
                      'RdGy':cm.RdGy,'RdBu':cm.RdBu,
                      'PuOr':cm.PuOr,'PiYG':cm.PiYG,'PRGn':cm.PRGn,'BrBG':cm.BrBG}

uniform_sequential = {'viridis':cm.viridis,'viridis_rev':cm.viridis.reversed(),
                      'plasma':cm.plasma,'plasma_rev':cm.plasma.reversed(),
                      'inferno':cm.inferno,'inferno_rev':cm.inferno.reversed(),
                      'magma':cm.magma,'magma_rev':cm.magma.reversed(),
                      'cividis':cm.cividis,'cividis_rev':cm.cividis.reversed()}

qualitative_colors = {'Pastel1':cm.Pastel1,'Pastel2':cm.Pastel2,'Paired':cm.Paired,
                      'Accent':cm.Accent,'Dark2':cm.Dark2,'Set1':cm.Set1,'Set2':cm.Set2,'Set3':cm.Set3,
                      'tab10':cm.tab10,'tab20':cm.tab20,'tab20b':cm.tab20b,'tab20c':cm.tab20c}

tab20 = [matplotlib.colors.to_hex(i) for i in cm.tab20(np.arange(20)/20.)]
tab20b = [matplotlib.colors.to_hex(i) for i in cm.tab20b(np.arange(20)/20.)]
tab20c = [matplotlib.colors.to_hex(i) for i in cm.tab20c(np.arange(20)/20.)]
Set3 = [matplotlib.colors.to_hex(i) for i in cm.Set3(np.arange(12)/12.)]
Set2 = [matplotlib.colors.to_hex(i) for i in cm.Set2(np.arange(8)/8.)]
Set1 = [matplotlib.colors.to_hex(i) for i in cm.Set1(np.arange(9)/9.)]
Pastel2 = [matplotlib.colors.to_hex(i) for i in cm.Pastel2(np.arange(8)/8.)]
Pastel1 = [matplotlib.colors.to_hex(i) for i in cm.Pastel1(np.arange(9)/9.)]
Dark2 = [matplotlib.colors.to_hex(i) for i in cm.Dark2(np.arange(8)/8.)]
Paired = [matplotlib.colors.to_hex(i) for i in cm.Paired(np.arange(12)/12.)]
Accent = [matplotlib.colors.to_hex(i) for i in cm.Accent(np.arange(8)/8.)]
Spectral = [matplotlib.colors.to_hex(i) for i in cm.Spectral(np.arange(11)/11.)]

edge_colors = {'Colormap1':tab20 + tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral,
               'Colormap2':Spectral + tab20 + tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent,
               'Colormap3':Accent + Spectral + tab20 + tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired,
               'Colormap4':Paired + Accent + Spectral + tab20 + tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2,
               'Colormap5':Dark2 + Paired + Accent + Spectral + tab20 + tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1,
               'Colormap6':Pastel1 + Dark2 + Paired + Accent + Spectral + tab20 + tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2,
               'Colormap7':Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral + tab20 + tab20b + tab20c + Set3 + Set2 + Set1,
               'Colormap8':Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral + tab20 + tab20b + tab20c + Set3 + Set2,
               'Colormap9':Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral + tab20 + tab20b + tab20c + Set3,
               'Colormap10':Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral + tab20 + tab20b + tab20c,
               'Colormap11':tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral + tab20 + tab20b,
               'Colormap12':tab20b + tab20c + Set3 + Set2 + Set1 + Pastel2 + Pastel1 + Dark2 + Paired + Accent + Spectral + tab20}


R_COLORS = {'GreBlaRed':['green', 'black', 'red'],
            'RedBlaGre':['red', 'black', 'green'],
            'RedBlaBlu':['red', 'black', 'blue'],
            'BluBlaRed':['blue', 'black', 'red'],
            'RedYlBlu':['red', 'yellow', 'blue'],
            'BluYlRed':['blue', 'yellow', 'red'],
            'BrBG':[matplotlib.colors.to_hex(i) for i in cm.BrBG(np.arange(11)/11.)],
            'PiYG':[matplotlib.colors.to_hex(i) for i in cm.PiYG(np.arange(11)/11.)],
            'PRGn':[matplotlib.colors.to_hex(i) for i in cm.PRGn(np.arange(11)/11.)],
            'PuOr':[matplotlib.colors.to_hex(i) for i in cm.PuOr(np.arange(11)/11.)],
            'RdBu':[matplotlib.colors.to_hex(i) for i in cm.RdBu(np.arange(11)/11.)],
            'RdGy':[matplotlib.colors.to_hex(i) for i in cm.RdGy(np.arange(11)/11.)],
            'RdYlBu':[matplotlib.colors.to_hex(i) for i in cm.RdYlBu(np.arange(11)/11.)],
            'RdYlGn':[matplotlib.colors.to_hex(i) for i in cm.RdYlGn(np.arange(11)/11.)],
            'Spectral':[matplotlib.colors.to_hex(i) for i in cm.Spectral(np.arange(11)/11.)],
            'Blues':[matplotlib.colors.to_hex(i) for i in cm.Blues(np.arange(9)/9.)],
            'BuGn':[matplotlib.colors.to_hex(i) for i in cm.BuGn(np.arange(9)/9.)],
            'BuPu':[matplotlib.colors.to_hex(i) for i in cm.BuPu(np.arange(9)/9.)],
            'GnBu':[matplotlib.colors.to_hex(i) for i in cm.GnBu(np.arange(9)/9.)],
            'Greens':[matplotlib.colors.to_hex(i) for i in cm.Greens(np.arange(9)/9.)],
            'Greys':[matplotlib.colors.to_hex(i) for i in cm.Greys(np.arange(9)/9.)],
            'Oranges':[matplotlib.colors.to_hex(i) for i in cm.Oranges(np.arange(9)/9.)],
            'OrRd':[matplotlib.colors.to_hex(i) for i in cm.OrRd(np.arange(9)/9.)],
            'PuBu':[matplotlib.colors.to_hex(i) for i in cm.PuBu(np.arange(9)/9.)],
            'PuBuGn':[matplotlib.colors.to_hex(i) for i in cm.PuBuGn(np.arange(9)/9.)],
            'PuRd':[matplotlib.colors.to_hex(i) for i in cm.PuRd(np.arange(9)/9.)],
            'Purples':[matplotlib.colors.to_hex(i) for i in cm.Purples(np.arange(9)/9.)],
            'RdPu':[matplotlib.colors.to_hex(i) for i in cm.RdPu(np.arange(9)/9.)],
            'Reds':[matplotlib.colors.to_hex(i) for i in cm.Reds(np.arange(9)/9.)],
            'YlGn':[matplotlib.colors.to_hex(i) for i in cm.YlGn(np.arange(9)/9.)],
            'YlGnBu':[matplotlib.colors.to_hex(i) for i in cm.YlGnBu(np.arange(9)/9.)],
            'YlOrBr':[matplotlib.colors.to_hex(i) for i in cm.YlOrBr(np.arange(9)/9.)],
            'YlOrRd':[matplotlib.colors.to_hex(i) for i in cm.YlOrRd(np.arange(9)/9.)]}


# para que no haya errores deben existir estas imagenes dentro del folder NeVOmics_img
imagenes = ['Accent.png', 'BluBlaRed.png', 'Blues.png','BluYlRed.png', 'BrBG.png', 'BuGn.png',
 'BuPu.png', 'bwr.png', 'cividis.png', 'cividis_rev.png', 'Colormap1.png', 'Colormap10.png',
 'Colormap11.png', 'Colormap12.png', 'Colormap2.png', 'Colormap3.png', 'Colormap4.png',
 'Colormap5.png', 'Colormap6.png', 'Colormap7.png', 'Colormap8.png', 'Colormap9.png',
 'colormap_figures.txt', 'coolwarm.png', 'Dark2.png', 'diverging_figures.txt', 'GnBu.png',
 'GreBlaRed.png', 'Greens.png', 'Greys.png', 'inferno.png', 'inferno_rev.png','magma.png',
 'magma_rev.png', 'Oranges.png', 'OrRd.png', 'Paired.png', 'Pastel1.png', 'Pastel2.png',
 'PdPu.png', 'PiYG.png', 'plasma.png', 'plasma_rev.png', 'PRGn.png', 'PuBu.png',
 'PuBuGn.png', 'PuOr.png', 'PuRd.png', 'Purples.png', 'qualitative_figures.txt', 'RdBu.png',
 'RdYlBu.png', 'RdGy.png', 'RdYlGn.png', 'RedBlaBlu.png', 'RedBlaGre.png', 'Reds.png',
 'RedYlBlu.png', 'seismic.png', 'sequentials_figures.txt', 'Set1.png', 'Set2.png', 'Set3.png',
 'Spectral.png', 'tab10.png', 'tab20.png', 'tab20b.png', 'tab20c.png', 'uniform_figures.txt',
 'viridis.png', 'viridis_rev.png', 'YlGn.png', 'YlGnBu.png', 'YlOrBr.png', 'YlOrRd.png']




def create_colormaps():
    x1 = 0
    x2 = .085 # ancho
    y2 = .0085 # alto
    fig = plt.figure(figsize=(5, 3))
    y1 = 0
    diverging_figures = []
    for i in diverging_colors:
        
        ax = fig.add_axes([x1, y1, x2, y2])
    
        cb = mpl.colorbar.ColorbarBase(ax, cmap=diverging_colors[i],orientation='horizontal')
        cb.outline.set_linewidth(0)
        ax.axis('off')
        ax.text(1.01, 0.5, i, size = 1, ha = 'left', va = 'center', color = 'black')
        
        y1+=y2 +0.001
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        diverging_figures.append('NeVOmics_img/'+i+'.png')
        plt.savefig('NeVOmics_img/'+i+'.png', dpi = 600, bbox_inches=extent, pad_inches = 0)
    plt.close()
    params = open('NeVOmics_img/diverging_figures.txt','w')
    params.write('\n'.join(diverging_figures))
    params.close()
    #.................................................................................................
    fig = plt.figure(figsize=(5, 3))
    y1 = 0
    sequentials_figures = []
    for i in sequentials_colors:
        
        ax = fig.add_axes([x1, y1, x2, y2])
        
        cb = mpl.colorbar.ColorbarBase(ax, cmap=sequentials_colors[i],orientation='horizontal')
        cb.outline.set_linewidth(0)
        ax.axis('off')
        ax.text(1.01, 0.5,i,size=1,ha='left',va='center',color= 'black')
        
        y1+=y2 +0.001
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        sequentials_figures.append('NeVOmics_img/'+i+'.png')
        plt.savefig('NeVOmics_img/'+i+'.png', dpi = 600, bbox_inches=extent, pad_inches = 0)
    plt.close()
    params = open('NeVOmics_img/sequentials_figures.txt','w')
    params.write('\n'.join(sequentials_figures))
    params.close()
    #....................................................................
    fig = plt.figure(figsize=(5, 3))
    y1 = 0
    uniform_figures = []
    for i in uniform_sequential:
        
        ax = fig.add_axes([x1, y1, x2, y2])
        
        cb = mpl.colorbar.ColorbarBase(ax, cmap=uniform_sequential[i],orientation='horizontal')
        cb.outline.set_linewidth(0)
        ax.axis('off')
        ax.text(1.01, 0.5, i, size = 1, ha = 'left', va = 'center', color = 'black')
        
        y1+=y2 +0.001
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        uniform_figures.append('NeVOmics_img/'+i+'.png')
        plt.savefig('NeVOmics_img/'+i+'.png', dpi = 600, bbox_inches=extent, pad_inches = 0)
    plt.close()
    params = open('NeVOmics_img/uniform_figures.txt','w')
    params.write('\n'.join(uniform_figures))
    params.close()
    #-------------------------------------------------------------------------------------------------
    fig = plt.figure(figsize=(5, 3))
    y1 = 0
    qualitative_figures = []
    for i in qualitative_colors:
        
        ax = fig.add_axes([x1, y1, x2, y2])
        
        cb = mpl.colorbar.ColorbarBase(ax, cmap=qualitative_colors[i],orientation='horizontal')
        cb.outline.set_linewidth(0)
        ax.axis('off')
        ax.text(1.01, 0.5, i, size = 1, ha = 'left', va = 'center', color = 'black')
        
        y1+=y2 +0.001
        extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        qualitative_figures.append('NeVOmics_img/'+i+'.png')
        plt.savefig('NeVOmics_img/'+i+'.png', dpi = 600, bbox_inches=extent, pad_inches = 0)
    plt.close()
    params = open('NeVOmics_img/qualitative_figures.txt','w')
    params.write('\n'.join(qualitative_figures))
    params.close()
    #-------------------------------------------------------------------------------------------------
    y1 = 0
    colormap_figures = []
    for k in edge_colors:
        fig = plt.figure(figsize=(5, 3))
        aa = fig.add_axes([x1, y1, x2, y2])
        n = 0
        for j in edge_colors[k][0:30]:
            plt.plot([n,n],[0,0.5], color=j,  linewidth=1.1)
            n += 1
            aa.margins(x=0.018)
            aa.axis('off')
            extent = aa.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        colormap_figures.append('NeVOmics_img/'+k+'.png')
        plt.savefig('NeVOmics_img/'+k+'.png', dpi = 600, bbox_inches=extent, pad_inches = 0)    
        plt.close()
    params = open('NeVOmics_img/colormap_figures.txt','w')
    params.write('\n'.join(colormap_figures))
    params.close()
    #-------------------------------------------------------------------------------------------------

if os.path.exists("NeVOmics_img") == True:
    faltantes = []
    for k in imagenes:
        if os.path.exists('NeVOmics_img/'+k) == True:
            continue
        else:
            print('The following file is missing:',k)
            faltantes.append(k)
    if len(faltantes) >= 1:
        create_colormaps()
else:
    ## Create a folder
    os.makedirs('NeVOmics_img')
    # descargar las imagenes de la interfase y depositarlas en esta carpeta
    create_colormaps()

#*urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/icon_nevomics.ico',
#*                           'icon_nevomics.ico')
cir = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/circo.png',
                           'NeVOmics_img/circo.png')
netw2 = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Network2.png',
                           'NeVOmics_img/Network2.png')
net0 = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/net.png',
                           'NeVOmics_img/net.png')


params = open('NeVOmics_params.txt','w')
params.close()

# default values
with open("NeVOmics_params.txt", "a") as f:
    f.write('edgecolor=Colormap8\n')
    f.write('networkcolor=GreBlaRed\n')
    f.write('usertext=LogFC\n')
    f.write('uniquecolor=#0000FF\n')
    f.close()
#-----------------------------------------   
    


opciones = [" Gene Ontology Enrichment                         ",
            " KEGG Pathways Enrichment                         ",
            " KEGG Blast Pathways Enrichment                   "]

print('NeVOmics History:')

root = Tk()
#root.attributes("-topmost", True)
root.title("NeVOmics")
#           ancho , alto
#root.geometry("1020x630")
#*root.iconbitmap(r'icon_nevomics.ico')


lab0 = Label(root, text="       ")# columna vacía
lab0.grid(column=0, row=0)

lab1 = Label(root, text="NeVOmics",
             font=("Courier New", 35, "bold"), underline=10,
             fg = 'green')
lab1.grid(column=1, columnspan=6, row=0, sticky= W+S+N)


uno = Label(root, text="1. Functional Analysis", font=("Courier New", 15, "bold"))
uno.grid(column=1, columnspan=4, row=2, sticky= W+N)

colors = {0:'black', 1:'black', 2:'black'} # bg="#000000", fg="white",
analysis = IntVar()#--------------------------------------------------------------
pos = {0:3,1:8,2:12}
for i, option in enumerate(opciones):
    radio = Radiobutton(root, text=option, font=("Courier New", 10, "bold"), bg = colors[i],
                        activebackground = 'yellow', activeforeground = 'black', selectcolor='black',
                        fg = 'white', cursor="hand2", variable=analysis, value=i)
    radio.grid(column=1, columnspan=3, row=pos[i], sticky= W+N)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_adj = LabelFrame(root, text="P-value adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_adj.grid(column=2, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_plots = LabelFrame(root, text = "Plots") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_plots.grid(column=3, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_aspect = LabelFrame(root, text = "Aspects") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect.grid(column=1, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#group_pval = LabelFrame(root, text = "Not adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_pval.grid(column=2, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=============================================================================
#oror = Label(root,text = " OR ", font = ("Courier New", 11, "bold"))
#oror.grid(column=3, row=5, sticky= W+E+S)

#yyy = Label(root, text=" & ", font=("Courier New", 12, "bold"))
#yyy.grid(column=5, row=5, sticky= W+E+S)
#=============================================================================

uno = Label(group_aspect, text = 'Biological Process: ', font=("Courier New", 10, "bold"))
uno.grid(column = 1, row = 4, sticky = E)

#dos = Label(group_pval, text= 'P-value: ', font=("Courier New", 10))
#dos.grid(column = 1, row = 4, sticky= S+W)

#bppval = StringVar()#--------------------------------------------------------------
#tres = ttk.Combobox(group_pval, textvariable = bppval, font=("Courier New", 10),
#                     values = Pvalues_valores, width=5)
#tres.grid(column=2, row=4, sticky= S+W)
#tres.current(14)

#vacia = Label(group_pval, text = ' ', font = ("Courier New", 10))
#vacia.grid(column=3, row=4, sticky= S+W)


cuatro = Label(group_adj, text= 'FDR: ', font=("Courier New", 10))
cuatro.grid(column=1, row=4, sticky= S+W)

bpz = StringVar()#--------------------------------------------------------------
cinco = ttk.Combobox(group_adj, textvariable = bpz, font=("Courier New", 10),
                     values = valores, width=4)
cinco.grid(column=2, row=4, sticky= S+W)
cinco.current(6)

seis = Label(group_adj, text=" % ", font=("Courier New", 10))
seis.grid(column=3, row=4, sticky= S+W)

bpgr = IntVar()#--------------------------------------------------------------
siete = Checkbutton(group_plots, text='', font=("Courier New", 10),activebackground= 'yellow',
                    cursor="hand2", variable=bpgr)
siete.grid(column = 1, row = 4, sticky = W)

######

uno1 = Label(group_aspect, text = 'Molecular Function: ', font=("Courier New", 10, "bold"))
uno1.grid(column = 1, row = 5, sticky = E)

#dos1 = Label(group_pval, text= 'P-value: ', font=("Courier New", 10))
#dos1.grid(column = 1, row = 5, sticky= S+W)

#mfpval = StringVar()#--------------------------------------------------------------
#tres1 = ttk.Combobox(group_pval, textvariable = mfpval, font=("Courier New", 10),
#                     values = Pvalues_valores, width=5)
#tres1.grid(column=2, row=5, sticky= S+W)
#tres1.current(14)

#vacia1 = Label(group_pval, text = ' ', font = ("Courier New", 10))
#vacia1.grid(column=3, row=5, sticky= S+W)


cuatro1 = Label(group_adj, text= 'FDR: ', font=("Courier New", 10))
cuatro1.grid(column=1, row=5, sticky= S+W)

mfz = StringVar()#--------------------------------------------------------------
cinco1 = ttk.Combobox(group_adj, textvariable = mfz, font=("Courier New", 10),
                     values = valores, width=4)
cinco1.grid(column=2, row=5, sticky= S+W)
cinco1.current(6)

seis1 = Label(group_adj, text=" % ", font=("Courier New", 10))
seis1.grid(column=3, row=5, sticky= S+W)

mfgr = IntVar()#--------------------------------------------------------------
siete1 = Checkbutton(group_plots, text='', font=("Courier New", 10),activebackground= 'yellow',
                    cursor="hand2", variable=mfgr)
siete1.grid(column = 1, row = 5, sticky = W)

######

uno2 = Label(group_aspect, text = 'Cellular Component: ', font=("Courier New", 10, "bold"))
uno2.grid(column = 1, row = 6, sticky = E)

#dos2 = Label(group_pval, text= 'P-value: ', font=("Courier New", 10))
#dos2.grid(column = 1, row = 6, sticky= S+W)

#ccpval = StringVar()#--------------------------------------------------------------
#tres2 = ttk.Combobox(group_pval, textvariable = ccpval, font=("Courier New", 10),
#                     values = Pvalues_valores, width=5)
#tres2.grid(column=2, row=6, sticky= S+W)
#tres2.current(14)

#vacia2 = Label(group_pval, text = ' ', font = ("Courier New", 10))
#vacia2.grid(column=3, row=6, sticky= S+W)


cuatro2 = Label(group_adj, text= 'FDR: ', font=("Courier New", 10))
cuatro2.grid(column=1, row=6, sticky= S+W)

ccz = StringVar()#--------------------------------------------------------------
cinco2 = ttk.Combobox(group_adj, textvariable = ccz, font=("Courier New", 10),
                     values = valores, width=4)
cinco2.grid(column=2, row=6, sticky= S+W)
cinco2.current(6)

seis2 = Label(group_adj, text=" % ", font=("Courier New", 10))
seis2.grid(column=3, row=6, sticky= S+W)

ccgr = IntVar()#--------------------------------------------------------------
siete2 = Checkbutton(group_plots, text='', font=("Courier New", 10),activebackground= 'yellow',
                    cursor="hand2", variable=ccgr)
siete2.grid(column = 1, row = 6, sticky = W)

######

#fila_vacia = Label(root, text=" ")# file vacía
#fila_vacia.grid(column=1, row=7)



#=====================================================================

anotacion = LabelFrame(root, text="Default GO Annotation  &  Complete GO Annotation") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
anotacion.grid(column=1, columnspan = 5, row=7, sticky = W)

fila_ann1 = Label(anotacion, text=" ")# file vacía
fila_ann1.grid(column=1, row=7)
fila_ann2 = Label(anotacion, text="           ")# file vacía
fila_ann2.grid(column=3, row=7)
fila_ann3 = Label(anotacion, text="             ")# file vacía
fila_ann3.grid(column=5, row=7)

check_annotation = Checkbutton(anotacion, text= "    UniProtKB    ",
                               state=DISABLED, bg = '#90ee90', fg = 'black')
check_annotation.grid(column=2, row=7)



# aviso al usuario de que la descarga de la anotación completa tardará en descargar
def advertencia():
    if aceptar.get() == 1:
        messagebox.showinfo('Status',
                            'The full annotation will be downloaded, this may '\
                            'take some time, approximately 1-3 h.\n\n'\
                            'Once the Complete Annotation is downloaded, it will no'\
                            ' longer be downloaded for later analysis, if you later '\
                            'wish to update the Complete Annotation (for example every month),'\
                            ' simply delete '\
                            'the file labeled with: Complete_Annotation_.*')
    if aceptar.get() == 0:
        pass

aceptar = IntVar()
check = Checkbutton(anotacion, text="    QuickGO    ", bg = '#90ee90', fg = 'black',
                    variable = aceptar,
                   command =  advertencia)
check.grid(column=4, row=7)

def aclaracion(event):
    webbrowser.open_new(r"https://www.uniprot.org/help/complete_go_annotation")
acla = Label(anotacion, text="?",
             font=("Courier New", 12, "bold"), fg="blue", cursor="hand2")
acla.grid(column=6, row=7)
acla.bind("<Button-1>", aclaracion)




#comentario = Button(anotacion, text="?", font=("Courier New", 10, "bold"))
#comentario.grid(column=6, row=7)

#======================================================================

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_adj1 = LabelFrame(root, text="P-value adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_adj1.grid(column=2, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_plots1 = LabelFrame(root, text = "Plots") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_plots1.grid(column=3, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_aspect1 = LabelFrame(root, text = "Classification") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect1.grid(column=1, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#group_pval1 = LabelFrame(root, text = "Not adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_pval1.grid(column=2, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=============================================================================



#=============================================================================

uno3 = Label(group_aspect1, text = '      KEGG Pathways: ', font=("Courier New", 10, "bold"))
uno3.grid(column = 1, row = 9, sticky = E)

#dos3 = Label(group_pval1, text= 'P-value: ', font=("Courier New", 10))
#dos3.grid(column = 1, row = 9, sticky= S+W)

#keggpval = StringVar()#--------------------------------------------------------------
#tres3 = ttk.Combobox(group_pval1, textvariable = keggpval, font=("Courier New", 10),
#                     values = Pvalues_valores, width=5)
#tres3.grid(column=2, row=9, sticky= S+W)
#tres3.current(14)

#vacia3 = Label(group_pval1, text = ' ', font = ("Courier New", 10))
#vacia3.grid(column=3, row=9, sticky= S+W)


cuatro3 = Label(group_adj1, text= 'FDR: ', font=("Courier New", 10))
cuatro3.grid(column=1, row=9, sticky= S+W)

keggz = StringVar()#--------------------------------------------------------------
cinco3 = ttk.Combobox(group_adj1, textvariable = keggz, font=("Courier New", 10),
                     values = valores, width=4)
cinco3.grid(column=2, row=9, sticky= S+W)
cinco3.current(6)

seis3 = Label(group_adj1, text=" % ", font=("Courier New", 10))
seis3.grid(column=3, row=9, sticky= S+W)

kegggr = IntVar()#--------------------------------------------------------------
siete3 = Checkbutton(group_plots1, text='', font=("Courier New", 10),activebackground= 'yellow',
                    cursor="hand2", variable=kegggr)
siete3.grid(column = 1, row = 9)

###

group_aspect11 = LabelFrame(root, text = "Reference") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect11.grid(column=1, row=10)

group_org = LabelFrame(root, text = "Organism-specific") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_org.grid(column=1, columnspan=3, row=10)


#ocho = Label(group_aspect11, text= '    KEGG Organisms: ', font=("Courier New", 10, "bold"))
#ocho.grid(column=1, row=10, sticky= E)    

org_kegg = StringVar()#--------------------------------------------------------------
nueve = ttk.Combobox(group_org, textvariable = org_kegg, font=("Courier New", 10),
                     cursor="hand2", values = kegg_organism.Organism.tolist(), height=20, width=50)
nueve.grid(column=1, row=10, sticky= W)
nueve.current(1)


######

fila_vacia1 = Label(root, text=" ")# file vacía
fila_vacia1.grid(column=1, row=11)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_adj2 = LabelFrame(root, text="P-value adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_adj2.grid(column=2, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_plots2 = LabelFrame(root, text = "Plots") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_plots2.grid(column=3, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_aspect2 = LabelFrame(root, text = "Classification") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect2.grid(column=1, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#group_pval2 = LabelFrame(root, text = "Not adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_pval2.grid(column=2, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=============================================================================
#oror2 = Label(root,text = " OR ", font = ("Courier New", 11, "bold"))
#oror2.grid(column=3, row=13, sticky= W+E+S)

#yyy2 = Label(root, text=" & ", font=("Courier New", 12, "bold"))
#yyy2.grid(column=5, row=13, sticky= W+E+S)
#=============================================================================

uno4 = Label(group_aspect2, text = '      KEGG Pathways: ', font=("Courier New", 10, "bold"))
uno4.grid(column = 1, row = 13, sticky = E)

#dos4 = Label(group_pval2, text= 'P-value: ', font=("Courier New", 10))
#dos4.grid(column = 1, row = 13, sticky= S+W)

#keggblpval = StringVar()#--------------------------------------------------------------
#tres4 = ttk.Combobox(group_pval2, textvariable = keggblpval, font=("Courier New", 10),
#                     values = Pvalues_valores, width=5)
#tres4.grid(column=2, row=13, sticky= S+W)
#tres4.current(14)

#vacia4 = Label(group_pval2, text = ' ', font = ("Courier New", 10))
#vacia4.grid(column=3, row=13, sticky= S+W)


cuatro4 = Label(group_adj2, text= 'FDR: ', font=("Courier New", 10))
cuatro4.grid(column=1, row=13, sticky= S+W)

keggblz = StringVar()#--------------------------------------------------------------
cinco4 = ttk.Combobox(group_adj2, textvariable = keggblz, font=("Courier New", 10),
                     values = valores, width=4)
cinco4.grid(column=2, row=13, sticky= S+W)
cinco4.current(6)

seis4 = Label(group_adj2, text=" % ", font=("Courier New", 10))
seis4.grid(column=3, row=13, sticky= S+W)

keggblgr = IntVar()#--------------------------------------------------------------
siete4 = Checkbutton(group_plots2, text='', font=("Courier New", 10), activebackground= 'yellow',
                    cursor="hand2", variable=keggblgr)
siete4.grid(column = 1, row = 13)

###

group_aspect22 = LabelFrame(root, text = "Reference") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect22.grid(column=1, row=14)

group_org1 = LabelFrame(root, text = "Organism-specific") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_org1.grid(column=1, columnspan=3, row=14)


#ocho1 = Label(group_aspect22, text= '    KEGG Organisms: ', font=("Courier New", 10, "bold"))
#ocho1.grid(column=1, row=14, sticky= E)    

org_keggbl = StringVar()#--------------------------------------------------------------
nueve1 = ttk.Combobox(group_org1, textvariable = org_keggbl, font=("Courier New", 10),
                     cursor="hand2", values = kegg_organism.Organism.tolist(), height=20, width=50)
nueve1.grid(column=1, row=14, sticky= W)
nueve1.current(1)

###

group_aspect22 = LabelFrame(root, text = "Tool") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect22.grid(column=1, row=15)

group_method = LabelFrame(root, text = "Method") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_method.grid(column=2, row=15)



diez = Label(group_aspect22, text= '             Local Blast: ', font=("Courier New", 10, "bold"))
diez.grid(column=1, row=15, sticky= E)

keggblmethod = IntVar()#--------------------------------------------------------------
tipos = ['Blastp', 'Blastx']
#mets = {0:2,1:3}
mets = {0:1,1:2}
#span = {0:1,1:1}
for i, tips in enumerate(tipos):
    once = Radiobutton(group_method, text=tips, font=("Courier New", 10), cursor="hand2",
                       activebackground = 'black', activeforeground = 'lightgray',
                       variable=keggblmethod, value=i)#.pack(anchor='sw')
    once.grid(column=mets[i], #columnspan = span[i],
                row=15, sticky= W)


######
columna_vacia = Label(root, text="      ")# file vacía
columna_vacia.grid(column=7, row=2)

######

#group_label = LabelFrame(root, text = "Method") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_label.grid(column=8, row=1, sticky= W)
######

#label = Label(root,
#             text="2. Threshold level of significance", font=("Courier New", 15, "bold"))
#label.grid(column=8, row=0, sticky= W+S)

#umbralelegido = IntVar()#--------------------------------------------------------------
#umbrales = ['P-value', 'FDR']
#mets = {0:8,1:9}
#for i, tips in enumerate(umbrales):
#    radio10 = Radiobutton(group_label, text = tips, font=("Courier New", 10), cursor="hand2",
#                         variable = umbralelegido, value=i)
#    radio10.grid(column=mets[i], row=1, sticky= W)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fila_vacia1 = Label(root, text="    ")# file vacía
fila_vacia1.grid(column=4, row=2)

############


def inputfile():
    seleccionado = ''.join(re.findall('[A-Z].*[a-z]', opciones[analysis.get()]))
    if seleccionado == 'KEGG Blast Pathways Enrichment':
        file_path = ()
        while file_path == ():
            file_path = filedialog.askopenfilename()
            if file_path == ():
                messagebox.showinfo('Status',
                                    'It must be a file separated by tabs and with extension .fasta')
            else:
                if file_path.split('.')[-1] == 'fasta':
                    new = "filelocation="+file_path+"\n"
                    op = open("NeVOmics_params.txt", "r")
                    op = op.read()
                    if re.findall('filelocation.*', op):
                        newfile = re.sub('filelocation.*\n', new, op)
                        params = open('NeVOmics_params.txt','w')
                        params.write(newfile)
                        params.close()
                    else:
                        with open("NeVOmics_params.txt", "a") as f:
                            f.write(new)
                            f.close()
                else:
                    messagebox.showinfo('Status',
                                        'It must be a file in fasta format and with extension .fasta')
                    file_path = ()
    else:
        file_path = ()
        while file_path == ():
            file_path = filedialog.askopenfilename()

            if file_path == ():
                messagebox.showinfo('Status',
                                    'It must be a file separated by tabs and with extension .tsv')
                #file_path = ()
            else:
                if file_path.split('.')[-1] == 'tsv':
                    new = "filelocation="+file_path+"\n"
                    op = open("NeVOmics_params.txt", "r")
                    op = op.read()
                    if re.findall('filelocation.*', op):
                        newfile = re.sub('filelocation.*\n', new, op)
                        params = open('NeVOmics_params.txt','w')
                        params.write(newfile)
                        params.close()
                    else:
                        with open("NeVOmics_params.txt", "a") as f:
                            f.write(new)
                            f.close()
                else:
                    messagebox.showinfo('Status',
                                       'It must be a file separated by tabs and with extension .tsv')
                    file_path = ()
                

subir_archivo = Label(root, text="2. Select File ", font=("Courier New", 15, "bold"))
subir_archivo.grid(column=5, row=0, sticky= W+S)
subir_archivo1 = Button(root, text="           Upload           ", bg="#000000", fg="white",
                activebackground= 'yellow',activeforeground= 'black',
                font=("Courier New", 10, "bold"), cursor="hand2", command=inputfile)
subir_archivo1.grid(column=5, row=1, sticky= W+N)





#..........................................................
"""
visualizaciones
"""

subir_archivo = Label(root, text="3. Visualizations", font=("Courier New", 15, "bold"))
subir_archivo.grid(column=5, row=2, sticky= W+S)

#%%%%%%%%%%%%%%%

plot_space = LabelFrame(root) # 
plot_space.grid(column=5, columnspan = 4, row=3, rowspan=8, sticky= W)



textonet1 = Label(plot_space, text="Networks:", font=("Courier New", 13, "bold"))
textonet1.grid(column=2, row=3, sticky= W)

def informacion(event):
    messagebox.showinfo('Information', 'Color of nodes for genes/proteins with associated numerical values.')
textonet2 = Label(plot_space, text="?",
             font=("Courier New", 12, "bold"), fg="blue", cursor="hand2")
textonet2.grid(column=3, row=4, sticky= W)
textonet2.bind("<Button-1>", informacion)

def informacion2(event):
    messagebox.showinfo('Information', 'Color of nodes for terms and edges.\n')
textonet3 = Label(plot_space, text="?",
             font=("Courier New", 12, "bold"), fg="blue", cursor="hand2")
textonet3.grid(column=3, row=7, sticky= W)
textonet3.bind("<Button-1>", informacion2)

def informacion3(event):
    messagebox.showinfo('Information', 'Color of nodes for genes/proteins without associated numerical values.')
textonet4 = Label(plot_space, text="?",
             font=("Courier New", 12, "bold"), fg="blue", cursor="hand2")
textonet4.grid(column=3, row=10, sticky= W)
textonet4.bind("<Button-1>", informacion3)

def informacion4(event):
    messagebox.showinfo('Information', 'Colormap title (default: LogFC).')
textonet5 = Label(plot_space, text="?",
             font=("Courier New", 12, "bold"), fg="blue", cursor="hand2")
textonet5.grid(column=3, row=11, sticky= W)
textonet5.bind("<Button-1>", informacion4)



image1 = tkinter.PhotoImage(file= 'NeVOmics_img/net.png')
netplot = IntVar()#--------------------------------------------------------------
siete4 = Checkbutton(plot_space, image=image1, activebackground= 'yellow',
                    cursor="hand2", variable=netplot)
siete4.grid(column = 3, columnspan=3, row = 3, sticky= W)


imagee = tkinter.PhotoImage(file= 'NeVOmics_img/GreBlaRed.png')
lbl = Label(plot_space, image=imagee)
lbl.grid(column=2, columnspan=4, row=5, rowspan=2, sticky= W)


def color_palettes():
    
    newwin = Toplevel()
    
    lab1 = Label(newwin, text="Python colors",
             font=("Courier New", 30, "bold"),
             fg = 'grey')
    lab1.grid(column=0, columnspan=4, row=0, sticky= W+S+N)


    def callback(event):
        webbrowser.open_new(r"https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html")
    link = Label(newwin, text="https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html\n",
             font=("Courier New", 8), fg="blue", cursor="hand2")
    link.grid(column=0, row=1, sticky= W)
    link.bind("<Button-1>", callback)
    
    
    diverging = LabelFrame(newwin, text="Diverging colormaps") 
    diverging.grid(column=0, columnspan = 2, row=2, rowspan=18, sticky= W)

    sequentials = LabelFrame(newwin, text = "Sequential colormaps")
    sequentials.grid(column=2, columnspan = 2, row=2, rowspan=18, sticky= W)

    uniforme = LabelFrame(newwin, text = "Uniform Sequential colormaps")
    uniforme.grid(column=4, columnspan = 2, row=2, rowspan=10, sticky= W) 
    
    div = open('NeVOmics_img/diverging_figures.txt', 'r')
    dive = div.read()
    div.close()
    diverging_figures = dive.split('\n')
    seq = open('NeVOmics_img/sequentials_figures.txt', 'r')
    sequ = seq.read()
    seq.close()
    sequentials_figures = sequ.split('\n')
    uni = open('NeVOmics_img/uniform_figures.txt', 'r')
    unif = uni.read()
    uni.close()
    uniform_figures = unif.split('\n')
    

    # diccionario para cambiar de root, y organizar dos conjuntos de valores
    clases = dict(zip(diverging_figures, np.repeat(diverging, len(diverging_figures))))
    clase2_seq = dict(zip(sequentials_figures, np.repeat(sequentials, len(sequentials_figures))))
    clase2_uni = dict(zip(uniform_figures, np.repeat(uniforme, len(uniform_figures))))
    clases.update(clase2_seq)
    clases.update(clase2_uni)
    

    columna = 0
    fila = 2
    color = IntVar()
    #color.set(0)
    l2 = []
    for i, j in enumerate(clases.keys()):
        image1 = tkinter.PhotoImage(file= j)
    
        l1 = tkinter.Radiobutton(clases[j], image=image1, cursor="hand2",
                   activebackground= 'white',activeforeground= 'white',
                         variable = color, value = i)
        
        l1.image = image1
        l1.grid(column = 1, row = fila)
        l2.append(l1)

        seis = Label(clases[j], text=j.split('/')[1].split('.')[0], font=("Courier New", 10, "bold"))
        seis.grid(column=0, row=fila, sticky= E)
        
    
        fila +=1
   
    fila_vacia = Label(newwin, text="   \n   ")# file vacía
    fila_vacia.grid(column=0, row=21)
    
    
    #-------------------------------------------------------------------------------
    def ShowChoice():
        colname = list(clases.keys())[color.get()].split('/')[1].split('.')[0]
        new = 'networkcolor='+colname+'\n' 
        op = open("NeVOmics_params.txt", "r")
        op = op.read()
        if re.findall('networkcolor.*', op):
            newfile = re.sub('networkcolor.*\n', new, op)
            params = open('NeVOmics_params.txt','w')
            params.write(newfile)
            params.close()
        else:
            with open("NeVOmics_params.txt", "a") as f:
                f.write(new)
                f.close()
        newwin.destroy()

        #-------------------------------------------------------------------------------
        #parametros = open('NeVOmics_params.txt', 'r')
        #parametros = parametros.read()
        #colormap_definido = re.search('networkcolor.*', parametros).group().split('=')[1]
        #for i in os.listdir('NeVOmics_img'):
        #    if re.search(colormap_definido, i):
        #        palette = i
        img = tkinter.PhotoImage(file= 'NeVOmics_img/'+colname+'.png')
        lbl.configure(image=img)
        lbl.image = img
        
        b1.configure(text=colname)
    
    boton = Button(newwin, text="    OK    ", cursor="hand2",
                activebackground= 'black',activeforeground= 'black',
                bg="green", fg="white",font=("Courier New", 12, "bold"),
                    command = ShowChoice) # newwin.destroy
    boton.grid(column=5, row=13, sticky= W)
    
    newwin.mainloop()




b1 = Button(plot_space, text='GreBlaRed', bg="#000000", fg="white", font=("Courier New", 10),
            command=color_palettes,
           activebackground = 'yellow', activeforeground = 'black')
b1.grid(column = 2, row = 4, sticky= W)


#fila_vaciax = Label(plot_space, text="  ")# file vacía
#fila_vaciax.grid(column=7, row=3)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# colormap para edges

imagee2 = tkinter.PhotoImage(file= 'NeVOmics_img/Colormap8.png')
lbl2 = Label(plot_space, image=imagee2)
lbl2.grid(column=2, columnspan=5, row=8, rowspan=2, sticky= W)


def edgeCOLORMAP():
    
    newwin = Toplevel()
    
    lab1 = Label(newwin, text="Python colors",
             font=("Courier New", 30, "bold"),
             fg = 'grey')
    lab1.grid(column=0, columnspan=4, row=0, sticky= W+S+N)


    def callback(event):
        webbrowser.open_new(r"https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html")
    link = Label(newwin, text="https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html\n",
             font=("Courier New", 8), fg="blue", cursor="hand2")
    link.grid(column=0, row=1, sticky= W)
    link.bind("<Button-1>", callback)
    
    
    qualitative = LabelFrame(newwin, text="Qualitative colormaps") 
    qualitative.grid(column=0, columnspan = 2, row=2, rowspan=18, sticky= W)
 
    qua = open('NeVOmics_img/colormap_figures.txt', 'r')
    qual = qua.read()
    qua.close()
    qualitative_figures = qual.split('\n')

    # diccionario para cambiar de root, y organizar dos conjuntos de valores
    clases = dict(zip(qualitative_figures, np.repeat(qualitative, len(qualitative_figures))))

    columna = 0
    fila = 2
    color = IntVar()
    #color.set(0)
    l2 = []
    for i, j in enumerate(clases.keys()):
        image1 = tkinter.PhotoImage(file= j)
    
        l1 = tkinter.Radiobutton(clases[j], image=image1, cursor="hand2",
                   activebackground= 'white',activeforeground= 'white',
                         variable = color, value = i)
        
        l1.image = image1
        l1.grid(column = 1, row = fila)
        l2.append(l1)

        seis = Label(clases[j], text=j.split('/')[1].split('.')[0], font=("Courier New", 10, "bold"))
        seis.grid(column=0, row=fila, sticky= E)
        
    
        fila +=1
    
    fila_vacia = Label(newwin, text="   \n   ")# file vacía
    fila_vacia.grid(column=0, row=21)
    #-------------------------------------------------------------------------------
    def ShowChoice():
        colorname = list(clases.keys())[color.get()].split('/')[1].split('.')[0]
        new = 'edgecolor='+colorname+'\n' 
        op = open("NeVOmics_params.txt", "r")
        op = op.read()
        if re.findall('edgecolor.*', op):
            newfile = re.sub('edgecolor.*\n', new, op)
            params = open('NeVOmics_params.txt','w')
            params.write(newfile)
            params.close()
        else:
            with open("NeVOmics_params.txt", "a") as f:
                f.write(new)
                f.close()
                    
        newwin.destroy()
            
        img2 = tkinter.PhotoImage(file= 'NeVOmics_img/'+colorname+'.png')
        lbl2.configure(image=img2)
        lbl2.image = img2
        
        colmap.configure(text=colorname)
    #-------------------------------------------------------------------------------
    
    boton = Button(newwin, text="    OK    ", cursor="hand2",
                activebackground= 'black',activeforeground= 'black',
                bg="green", fg="white",font=("Courier New", 12, "bold"),
                    command = ShowChoice) # newwin.destroy
    boton.grid(column=5, row=10, sticky= W)
    newwin.mainloop()



colmap = Button(plot_space, text='Colormap8', bg="#000000", fg="white", font=("Courier New", 10),
                command= edgeCOLORMAP, activebackground = 'yellow', activeforeground = 'black')
colmap.grid(column = 2, row = 7, sticky= W)

####################33

# presentar colores individuales #######################
# modifica el archivo de parámetros, se actualiza

def color_button():
    col = tkinter.colorchooser.askcolor(parent=root)
    sinback.configure(bg=col[1])
    sinback.configure(text=col[1])
    if col[0] == None:
        new = 'uniquecolor=#0000FF\n' # default = salmon
        op = open("NeVOmics_params.txt", "r")
        op = op.read()
        if re.findall('uniquecolor.*', op):
            newfile = re.sub('uniquecolor.*\n', new, op)
            params = open('NeVOmics_params.txt','w')
            params.write(newfile)
            params.close()
        else:
            with open("NeVOmics_params.txt", "a") as f:
                f.write(new)
                f.close()
    else:
        new = 'uniquecolor='+col[1]+"\n"
        op = open("NeVOmics_params.txt", "r")
        op = op.read()
        if re.findall('uniquecolor.*', op):
            newfile = re.sub('uniquecolor.*\n', new, op)
            params = open('NeVOmics_params.txt','w')
            params.write(newfile)
            params.close()
        else:
            with open("NeVOmics_params.txt", "a") as f:
                f.write(new)
                f.close()

                
                
sinback = tkinter.Button(plot_space, text='#0000ff', font=("Courier New", 10), bg= '#0000ff',
command=color_button)
sinback.grid(column=2, row=10, sticky= W)



def USERTEXT():
    newwin = Toplevel()
    
    e = Entry(newwin, bd =3) # para introducir titulo de la barra
    e.grid(column = 0, columnspan=2, row = 0, sticky= W)
    e.focus_set()
    
    def ShowChoice():
        
        new = 'usertext='+e.get()+'\n' 
        op = open("NeVOmics_params.txt", "r")
        op = op.read()
        if re.findall('usertext.*', op):
            newfile = re.sub('usertext.*\n', new, op)
            params = open('NeVOmics_params.txt','w')
            params.write(newfile)
            params.close()
        else:
            with open("NeVOmics_params.txt", "a") as f:
                f.write(new)
                f.close()
        if e.get() == '':
            new = 'usertext=LogFC\n' 
            op = open("NeVOmics_params.txt", "r")
            op = op.read()
            newfile = re.sub('usertext.*\n', 'usertext=LogFC\n', op)
            params = open('NeVOmics_params.txt','w')
            params.write(newfile)
            params.close()
                    
        newwin.destroy()
    #-------------------------------------------------------------------------------
    
        
    boton = Button(newwin, text="OK", cursor="hand2",
                activebackground= 'black',activeforeground= 'black',
                bg="green", fg="white",font=("Courier New", 10, "bold"),
                    command = ShowChoice) # newwin.destroy
    boton.grid(column=0, row=1, sticky= E)
    newwin.mainloop()
    
    

changetitle = Button(plot_space, text=" Title ", bg="#000000", fg="white", font=("Courier New", 10),
                     activebackground = 'yellow', activeforeground = 'black', command=USERTEXT)
changetitle.grid(column = 2, row = 11, sticky= W)


#######   CIRCOS

plot_space2 = LabelFrame(root) # 
plot_space2.grid(column=5, columnspan = 4, row=10, rowspan=5, sticky= W)

def informacion6(event):
    messagebox.showinfo('Information', 'To use this application you need to install R\n(Check the Installations.txt file).')
textonet6 = Label(plot_space2, text="! Warning",
             font=("Courier New", 12, "bold"), fg="red", cursor="hand2")
textonet6.grid(column=0, row=10, sticky= W)
textonet6.bind("<Button-1>", informacion6)


def informacion7(event):
    messagebox.showinfo('Information', 'The same colors of Networks will be used in Chord Plot.\nYou can modify the parameters above to configure this plot.')
textonet7 = Label(plot_space2, text=" ?",
             font=("Courier New", 12, "bold"), fg="blue", cursor="hand2")
textonet7.grid(column=1, row=10, sticky= W)
textonet7.bind("<Button-1>", informacion7)

textonet11 = Label(plot_space2, text="Chords: ", font=("Courier New", 13, "bold"))
textonet11.grid(column=0, row=11, sticky= W)

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,


#textonet22 = Label(plot_space2, text="Locate the R.exe program:", font=("Courier New", 8))
#textonet22.grid(column=5,columnspan=2, row=14, sticky= W)

#textonet33 = Label(plot_space2, text="Locate the Rlibrary_NeVOmics\ndirectory:", font=("Courier New", 8))
#textonet33.grid(column=5,columnspan=2, row=15, sticky= W)

image2 = tkinter.PhotoImage(file= 'NeVOmics_img/circo.png')
circoplot = IntVar()#--------------------------------------------------------------
siete5 = Checkbutton(plot_space2, image=image2, activebackground = 'yellow',
                    cursor="hand2", variable=circoplot)
siete5.grid(column = 1, row = 11, sticky= W)



### botones para localizar r

#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

#########


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

############




############
group_label5 = LabelFrame(root, text = "Identifier") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_label5.grid(column=10, row=1, rowspan=2, sticky= W)
######

label = Label(root,
             text="4. Node label", font=("Courier New", 15, "bold"))
label.grid(column=10, row=0, sticky= W)

node_lab = IntVar()#--------------------------------------------------------------
etiquetas = ['Gene Name', 'UniProt ID']
mets = {0:10,1:11}

for i, tips in enumerate(etiquetas):
    radio7 = Radiobutton(group_label5, text = tips, font=("Courier New", 10), cursor="hand2",
                         activebackground = 'black', activeforeground = 'lightgray',
                         variable = node_lab, value=i)
    radio7.grid(column=mets[i], row=1, sticky= W)

########

def parameters():
    
    boxplots = (1 in [bpgr.get(), mfgr.get(), ccgr.get(), kegggr.get(), keggblgr.get()])
    tipoplots = (1 in [netplot.get(), circoplot.get()])
    
    if boxplots == True and tipoplots == False:
        #print('seleccionar al menos uno')
        messagebox.showinfo('Status',
                            'You have activated at least one Plots box, therefore it is required to choose at least one type of graphic: Networks or Chords.')
    elif boxplots == True and tipoplots == True:
        #print('selecionaron box plot, y también seleccionarion un tipo')
        # esta opción se para correr el programa con al menos un tipo de gráfico ***********
        if circoplot.get() == 1: # solo si eligen la aplicacion "Chords"
            locrexe = os.path.exists("NeVOmics_locRexe.txt")
            locrlib = os.path.exists("NeVOmics_locRlib.txt")
            if True == True:
                #print('listo para ejecutar, eligieron chords, si no eligieron networks el script no los hace, pero si lo eligieron se harán los dos tipos')
                #-------------------------------->
                file = open('NeVOmics_params.txt', 'r')
                if bool(re.search('filelocation.*\n', file.read())) == True:
                    #print('listo para ejecutar') # --------------->
                    file = open('NeVOmics_params.txt', 'r')
                    file = file.read()
                    # para heredar los parametros previos
                    uno = "".join(re.findall('filelocation.*\n', file))
                    dos = "".join(re.findall('edgecolor.*\n', file))
                    tres = "".join(re.findall('networkcolor.*\n', file))
                    cuatro = "".join(re.findall('usertext.*\n', file))
                    cinco = "".join(re.findall('uniquecolor.*\n', file))
                    #######                           
                    ## Control of Uniprot and GOA directories
                    ls=[]
                    for i in os.listdir("./"):
                        if re.search('job_NeVOmics_[0-9]{1,3}',str(i)):
                            ls.append(i)
                    if ls == []:
                        new_folder='job_NeVOmics_1'
                        xoxo='1'
                        os.makedirs('job_NeVOmics_1',exist_ok=True)
                        print(new_folder)

                    else:
                        n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
                        xoxo=str(n+1)
                        n=str(n)
                        old_folder=''.join(re.findall('job_NeVOmics_'+n,str(ls)))
                        new_folder=re.sub(n,xoxo,old_folder)
                        os.makedirs(new_folder,exist_ok=True)
                        print(new_folder)

                    def abre_folder():
                        os.system('xdg-open '+new_folder)

                    label = Button(root, text= new_folder, font=("Courier New", 10, "bold"),cursor="hand2",
                    command = abre_folder, activebackground= 'yellow')
                    label.grid(column = 10, row = 5, rowspan=3, sticky= W)
        
                    pars = ["#====================\n"\
                            "analysis="+' '.join(re.findall('\w+', opciones[analysis.get()]))+"\n"\
                            "anotacion_goa="+str(aceptar.get())+"\n"\
                            "#=====\n"\
                            "#**GO ENRICHMENT**\n"\
                            "bpfdr="+str(bpz.get())+"\n"\
                            "mffdr="+str(mfz.get())+"\n"\
                            "ccfdr="+str(ccz.get())+"\n"\
                            "bpplots="+str(bpgr.get())+"\n"\
                            "mfplots="+str(mfgr.get())+"\n"\
                            "ccplots="+str(ccgr.get())+"\n"\
                            "#====================\n"\
                            "#**KEGG ENRICHMENT**\n"\
                            "keggfdr="+str(keggz.get())+"\n"\
                            "keggplots="+str(kegggr.get())+"\n"\
                            "keggorganism="+org_kegg.get()+"\n"\
                            "keggprefix="+dict_org[org_kegg.get()][0]+"\n"\
                            "keggTnumber="+dict_org[org_kegg.get()][1]+"\n"\
                            "#====================\n"\
                            "#**KEGG BLAST ENRICHMENT**\n"\
                            "keggblastfdr="+str(keggblz.get())+"\n"\
                            "keggblastplots="+str(keggblgr.get())+"\n"\
                            "keggmethodblast="+tipos[keggblmethod.get()]+"\n"\
                            "keggblastorganism="+org_keggbl.get()+"\n"\
                            "keggblastprefix="+dict_org[org_keggbl.get()][0]+"\n"\
                            "keggblastTnumber="+dict_org[org_keggbl.get()][1]+"\n"\
                            "#====================\n"\
                            "labelnode="+etiquetas[node_lab.get()]+"\n"\
                            "networksplots="+str(netplot.get())+"\n"\
                            "circosplots="+str(circoplot.get())+"\n"]
                        
                    params = open('NeVOmics_params.txt','w')
                    params.close()
                    with open("NeVOmics_params.txt", "a") as f:
                        f.write(uno+dos+tres+cuatro+cinco+''.join(pars))
                        f.close()
                    
                    shutil.copyfile('NeVOmics_params.txt', new_folder+'/NeVOmics_params.txt')
                    analisis = ' '.join(re.findall('\w+', opciones[analysis.get()]))
            
                    if analisis == 'Gene Ontology Enrichment':
                        print('Run: Gene Ontology Enrichment')
                        print(uno)
                        go_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_GO_Ubuntu.py',
                                                                 new_folder+'/short_GO_Ubuntu.py')
                        #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                        comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                        run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_GO_Ubuntu.py"', shell = True)
                        run.wait()
                        #os.system("start cmd /k cd "+comando+ " ^&^& python short_GO.py")
                    if analisis == 'KEGG Pathways Enrichment':
                        print('Run: KEGG Pathways Enrichment')
                        print(uno)
                        kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Ubuntu.py',
                                                                 new_folder+'/short_KEGG_Ubuntu.py')
                        #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                        comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                        run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_KEGG_Ubuntu.py"', shell = True)
                        run.wait()
                        #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG.py")
                        # >>>>>>>> fin
                    if analisis == 'KEGG Blast Pathways Enrichment':
                        print('Run: KEGG Blast Pathways Enrichment')
                        print(uno)
                        #print('!!! At the moment this analysis is in maintenance !!!')
                        kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Ubuntu_Blast.py',
                                                                 new_folder+'/short_KEGG_Ubuntu_Blast.py')
                        #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                        comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                        run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_KEGG_Ubuntu_Blast.py"', shell = True)
                        run.wait()
                    
                    
                else:
                    messagebox.showinfo('Status',
                                        'No file were selected.')
            else:
                messagebox.showinfo('Status',
                        'The location of R.exe programm or Rlibrary_NeVOmics directory is not yet defined e.g.\nC:/Users/home/Documents/R-3.5.3/bin/R.exe\n'+\
                        'C:/Users/home/Downloads/Rlibrary_NeVOmics.')
        elif netplot.get() == 1:
            #print('aqui la opcion si solo decidieron hacer hacer networks, sin chords, aqui no usa R')
            #--------------------------------->
            file = open('NeVOmics_params.txt', 'r')
            if bool(re.search('filelocation.*\n', file.read())) == True:
                #print('listo para ejecutar')# --------------->
                file = open('NeVOmics_params.txt', 'r')
                file = file.read()
                # para heredar los parametros previos
                uno = "".join(re.findall('filelocation.*\n', file))
                dos = "".join(re.findall('edgecolor.*\n', file))
                tres = "".join(re.findall('networkcolor.*\n', file))
                cuatro = "".join(re.findall('usertext.*\n', file))
                cinco = "".join(re.findall('uniquecolor.*\n', file))
                #######                           
                ## Control of Uniprot and GOA directories
                ls=[]
                for i in os.listdir("./"):
                    if re.search('job_NeVOmics_[0-9]{1,3}',str(i)):
                        ls.append(i)
                if ls == []:
                    new_folder='job_NeVOmics_1'
                    xoxo='1'
                    os.makedirs('job_NeVOmics_1',exist_ok=True)
                    print(new_folder)
                else:
                    n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
                    xoxo=str(n+1)
                    n=str(n)
                    old_folder=''.join(re.findall('job_NeVOmics_'+n,str(ls)))
                    new_folder=re.sub(n,xoxo,old_folder)
                    os.makedirs(new_folder,exist_ok=True)
                    print(new_folder)
                
                def abre_folder():
                    os.system('xdg-open '+new_folder)

                label = Button(root, text= new_folder, font=("Courier New", 10, "bold"),cursor="hand2",
                command = abre_folder, activebackground= 'yellow')
                label.grid(column = 10, row = 5, rowspan=3, sticky= W)
        
                pars = ["#====================\n"\
                        "analysis="+' '.join(re.findall('\w+', opciones[analysis.get()]))+"\n"\
                        "anotacion_goa="+str(aceptar.get())+"\n"\
                        "#=====\n"\
                        "#**GO ENRICHMENT**\n"\
                        "bpfdr="+str(bpz.get())+"\n"\
                        "mffdr="+str(mfz.get())+"\n"\
                        "ccfdr="+str(ccz.get())+"\n"\
                        "bpplots="+str(bpgr.get())+"\n"\
                        "mfplots="+str(mfgr.get())+"\n"\
                        "ccplots="+str(ccgr.get())+"\n"\
                        "#====================\n"\
                        "#**KEGG ENRICHMENT**\n"\
                        "keggfdr="+str(keggz.get())+"\n"\
                        "keggplots="+str(kegggr.get())+"\n"\
                        "keggorganism="+org_kegg.get()+"\n"\
                        "keggprefix="+dict_org[org_kegg.get()][0]+"\n"\
                        "keggTnumber="+dict_org[org_kegg.get()][1]+"\n"\
                        "#====================\n"\
                        "#**KEGG BLAST ENRICHMENT**\n"\
                        "keggblastfdr="+str(keggblz.get())+"\n"\
                        "keggblastplots="+str(keggblgr.get())+"\n"\
                        "keggmethodblast="+tipos[keggblmethod.get()]+"\n"\
                        "keggblastorganism="+org_keggbl.get()+"\n"\
                        "keggblastprefix="+dict_org[org_keggbl.get()][0]+"\n"\
                        "keggblastTnumber="+dict_org[org_keggbl.get()][1]+"\n"\
                        "#====================\n"\
                        "labelnode="+etiquetas[node_lab.get()]+"\n"\
                        "networksplots="+str(netplot.get())+"\n"\
                        "circosplots="+str(circoplot.get())+"\n"]
                        
                params = open('NeVOmics_params.txt','w')
                params.close()
                with open("NeVOmics_params.txt", "a") as f:
                    f.write(uno+dos+tres+cuatro+cinco+''.join(pars))
                    f.close()
                
                shutil.copyfile('NeVOmics_params.txt', new_folder+'/NeVOmics_params.txt')
                analisis = ' '.join(re.findall('\w+', opciones[analysis.get()]))
            
                if analisis == 'Gene Ontology Enrichment':
                    print('Run: Gene Ontology Enrichment')
                    print(uno)
                    go_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_GO_Ubuntu.py',
                                                             new_folder+'/short_GO_Ubuntu.py')
                    #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                    comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                    run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_GO_Ubuntu.py"', shell = True)
                    run.wait()
                    #os.system("start cmd /c cd "+comando+ " ^&^& python short_GO.py")
                if analisis == 'KEGG Pathways Enrichment':
                    print('Run: KEGG Pathways Enrichment')
                    print(uno)
                    kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Ubuntu.py',
                                                             new_folder+'/short_KEGG_Ubuntu.py')
                    #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                    comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                    run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_KEGG_Ubuntu.py"', shell = True)
                    run.wait()
                    #os.system("start cmd /c cd "+comando+ " ^&^& python short_KEGG.py")
                    # >>>>>>>> fin
                if analisis == 'KEGG Blast Pathways Enrichment':
                    print('Run: KEGG Blast Pathways Enrichment')
                    print(uno)
                    #print('!!! At the moment this analysis is in maintenance !!!')
                    kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Ubuntu_Blast.py',
                                                                new_folder+'/short_KEGG_Ubuntu_Blast.py')
                    #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                    comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                    run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_KEGG_Ubuntu_Blast.py"', shell = True)
                    run.wait()
                
                
                
            else:
                messagebox.showinfo('Status',
                                    'No file were selected.')
    elif boxplots == False and tipoplots == False:
        #print('no seleccionaron nada, no se harán gráficos')
        # esta opción es para correr el programa sin generar ningún tipo de gráficos ************
        #------------------------------>
        file = open('NeVOmics_params.txt', 'r')
        if bool(re.search('filelocation.*\n', file.read())) == True:
            print('listo para ejecutar') # --------------->
            file = open('NeVOmics_params.txt', 'r')
            file = file.read()
            # para heredar los parametros previos
            uno = "".join(re.findall('filelocation.*\n', file))
            dos = "".join(re.findall('edgecolor.*\n', file))
            tres = "".join(re.findall('networkcolor.*\n', file))
            cuatro = "".join(re.findall('usertext.*\n', file))
            cinco = "".join(re.findall('uniquecolor.*\n', file))
            #######                           
            ## Control of Uniprot and GOA directories
            ls=[]
            for i in os.listdir("./"):
                if re.search('job_NeVOmics_[0-9]{1,3}',str(i)):
                    ls.append(i)
            if ls == []:
                new_folder='job_NeVOmics_1'
                xoxo='1'
                os.makedirs('job_NeVOmics_1',exist_ok=True)
                print(new_folder)
            else:
                n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
                xoxo=str(n+1)
                n=str(n)
                old_folder=''.join(re.findall('job_NeVOmics_'+n,str(ls)))
                new_folder=re.sub(n,xoxo,old_folder)
                os.makedirs(new_folder,exist_ok=True)
                print(new_folder)
            
            def abre_folder():
                os.system('xdg-open '+new_folder)

            label = Button(root, text= new_folder, font=("Courier New", 10, "bold"),cursor="hand2",
            command = abre_folder, activebackground= 'yellow')
            label.grid(column = 10, row = 5, rowspan=3, sticky= W)
        
            pars = ["#====================\n"\
                    "analysis="+' '.join(re.findall('\w+', opciones[analysis.get()]))+"\n"\
                    "anotacion_goa="+str(aceptar.get())+"\n"\
                    "#=====\n"\
                    "#**GO ENRICHMENT**\n"\
                    "bpfdr="+str(bpz.get())+"\n"\
                    "mffdr="+str(mfz.get())+"\n"\
                    "ccfdr="+str(ccz.get())+"\n"\
                    "bpplots="+str(bpgr.get())+"\n"\
                    "mfplots="+str(mfgr.get())+"\n"\
                    "ccplots="+str(ccgr.get())+"\n"\
                    "#====================\n"\
                    "#**KEGG ENRICHMENT**\n"\
                    "keggfdr="+str(keggz.get())+"\n"\
                    "keggplots="+str(kegggr.get())+"\n"\
                    "keggorganism="+org_kegg.get()+"\n"\
                    "keggprefix="+dict_org[org_kegg.get()][0]+"\n"\
                    "keggTnumber="+dict_org[org_kegg.get()][1]+"\n"\
                    "#====================\n"\
                    "#**KEGG BLAST ENRICHMENT**\n"\
                    "keggblastfdr="+str(keggblz.get())+"\n"\
                    "keggblastplots="+str(keggblgr.get())+"\n"\
                    "keggmethodblast="+tipos[keggblmethod.get()]+"\n"\
                    "keggblastorganism="+org_keggbl.get()+"\n"\
                    "keggblastprefix="+dict_org[org_keggbl.get()][0]+"\n"\
                    "keggblastTnumber="+dict_org[org_keggbl.get()][1]+"\n"\
                    "#====================\n"\
                    "labelnode="+etiquetas[node_lab.get()]+"\n"\
                    "networksplots="+str(netplot.get())+"\n"\
                    "circosplots="+str(circoplot.get())+"\n"]
                        
            params = open('NeVOmics_params.txt','w')
            params.close()
            with open("NeVOmics_params.txt", "a") as f:
                f.write(uno+dos+tres+cuatro+cinco+''.join(pars))
                f.close()
        
            shutil.copyfile('NeVOmics_params.txt', new_folder+'/NeVOmics_params.txt')
            analisis = ' '.join(re.findall('\w+', opciones[analysis.get()]))
            
            if analisis == 'Gene Ontology Enrichment':
                print('Run: Gene Ontology Enrichment')
                print(uno)
                go_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_GO_Ubuntu.py',
                                                         new_folder+'/short_GO_Ubuntu.py')
                #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_GO_Ubuntu.py"', shell = True)
                run.wait()
                #os.system("start cmd /c cd "+comando+ " ^&^& python short_GO.py")
            if analisis == 'KEGG Pathways Enrichment':
                print('Run: KEGG Pathways Enrichment')
                print(uno)
                kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Ubuntu.py',
                                                         new_folder+'/short_KEGG_Ubuntu.py')
                #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_KEGG_Ubuntu.py"', shell = True)
                run.wait()
                #os.system("start cmd /c cd "+comando+ " ^&^& python short_KEGG.py")
                # >>>>>>>> fin
            if analisis == 'KEGG Blast Pathways Enrichment':
                print('Run: KEGG Blast Pathways Enrichment')
                print(uno)
                #print('!!! At the moment this analysis is in maintenance !!!')
                kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Ubuntu_Blast.py',
                                                            new_folder+'/short_KEGG_Ubuntu_Blast.py')
                #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                run = subprocess.Popen('cd '+comando+'/ && gnome-terminal -q -e "python3 short_KEGG_Ubuntu_Blast.py"', shell = True)
                run.wait()
            
            
            
            
        else:
            messagebox.showinfo('Status',
                                'No file were selected.')
    elif boxplots == False and tipoplots == True:
        #print('seleccionaron un tipo pero no seleccionaron un box plots')
        messagebox.showinfo('Status',
                            'You have activated at least one type of graphic: Networks or Chords, therefore it is required to choose at least one Plots box.')


#-----------------------------------------

botonrun = tkinter.Button(root, text = ' RUN ', bg="green", fg="white",
                font=("Courier New", 20, "bold"),
                          command = parameters,
                activebackground= 'yellow',activeforeground= 'black', cursor="hand2")
botonrun.grid(column = 10, row = 4, rowspan=2, sticky= W+E)



fila_vaciax1 = Label(plot_space, text="  ")# file vacía
fila_vaciax1.grid(column=10, row=7)

####>>>>>
ejemplo = LabelFrame(root) #, text = "Example Network created by NeVOmics" 
ejemplo.grid(column=10, columnspan = 1, row=8, rowspan=8, sticky= W+E)

imagenet = tkinter.PhotoImage(file= 'NeVOmics_img/Network2.png')
    
labnet = Label(ejemplo, image=imagenet)
        
labnet.image = imagenet
labnet.grid(column = 10, row = 9, rowspan=8, sticky= W+E)




def callback(event):
    webbrowser.open_new(r"https://doi.org/10.3390/genes9120569")
link = Label(root, text="If you use NeVOmics, please cite.  ",
             font=("Courier New", 10), fg="blue", cursor="hand2")
link.grid(column=10, row=16, sticky= W, columnspan=3)
link.bind("<Button-1>", callback)

final = Label(root, text = '  ', font = ("Courier New", 10))
final.grid(column=10, row=17)

root.mainloop()
