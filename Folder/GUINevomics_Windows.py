import subprocess
import re
import os

"""
os.makedirs('NeVOmics_PyMod',exist_ok=True)

modulos = {'numpy':'numpy==1.18.5',
           'pandas':'pandas==1.0.4', # 0.24.2
           'matplotlib':'matplotlib==3.0.3',
           'scipy':'scipy==1.4.1',
           'openpyxl':'openpyxl==3.0.3',
           'colormap':'colormap==1.0.2',
           'easydev':'easydev==0.9.38',
           'networkx':'networkx==2.2',
           'bioservices':'bioservices==1.6.0',
           'xlsxwriter':'xlsxwriter==1.2.7',
           'requests':'requests==2.23.0'} # 2.22.0


listmodules = subprocess.Popen('python -mpip list --path NeVOmics_PyMod/', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
lista = listmodules.communicate()[0].decode()
for modname in list(modulos.keys()):
    if modulos[modname] in [re.sub(' * ', '==', i.rstrip().lower()) for i in lista.split('\n')]:
        pass
    else:
        print('>>>>>>>>>>>>>>', modname, '|  Version: ', modulos[modname])
        subprocess.call('python -mpip install --force-reinstall --no-color '+modulos[modname]+' -t NeVOmics_PyMod/', shell = True)


import sys
sys.path.append("NeVOmics_PyMod/")
"""
modulos = {'numpy':'numpy==1.18.5',
           'pandas':'pandas==1.0.4', # 0.24.2
           'matplotlib':'matplotlib==3.0.3',
           'scipy':'scipy==1.4.1',
           'openpyxl':'openpyxl==3.0.3',
           'colormap':'colormap==1.0.2',
           'easydev':'easydev==0.9.38',
           'networkx':'networkx==2.2',
           'bioservices':'bioservices==1.6.0',
           'xlsxwriter':'xlsxwriter==1.2.7',
           'requests':'requests==2.23.0'} # 2.22.0

listmodules = subprocess.Popen('python -mpip list', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
lista = listmodules.communicate()[0].decode()
for modname in list(modulos.keys()):
    if modulos[modname] in [re.sub(' * ', '==', i.rstrip().lower()) for i in lista.split('\n')]:
        pass
    else:
        subprocess.call('python -mpip install --force-reinstall --no-color '+modulos[modname], shell = True)



print('\n______________________________________________\n')
print('STARTING NEVOMICS\n')
import requests
from matplotlib import cm
from colormap import Colormap
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import urllib.request
import webbrowser
from tkinter import ttk
from tkinter import *
import tkinter
from tkinter import messagebox
from tkinter import filedialog
import tkinter.colorchooser
import shutil


valores = []
n = 0.005
while n < 0.205:
    valores.append(str(round(n * 100, 1)))
    n += 0.005
valores = ['0.05', '0.1', '0.2', '0.3', '0.4'] + valores
valores.reverse()

Pvalues_valores = []
n = 0.005
while n < 0.205:
    Pvalues_valores.append(str(round(n, 10)))
    n += 0.005
Pvalues_valores = ['0.0005', '0.001' , '0.002' , '0.003' , '0.004'] + Pvalues_valores

# about KEGG





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

if os.path.exists('NeVOmics_img/icon_nevomics.ico') == True:
    pass
else:
    urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/icon_nevomics.ico',
                               'NeVOmics_img/icon_nevomics.ico')
if os.path.exists('NeVOmics_img/circo.png') == True:
    pass
else:
    urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/circo.png',
                               'NeVOmics_img/circo.png')
if os.path.exists('NeVOmics_img/Network2.png') == True:
    pass
else:
    urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Network2.png',
                               'NeVOmics_img/Network2.png')
if os.path.exists('NeVOmics_img/net.png') == True:
    pass
else:
    urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/net.png',
                               'NeVOmics_img/net.png')
if os.path.exists('NeVOmics_img/tkHyperlinkManager.py') == True:
    pass
else:
    urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/tkHyperlinkManager.py',
                               'NeVOmics_img/tkHyperlinkManager.py')

import sys
sys.path.append("NeVOmics_img/")
from tkHyperlinkManager import HyperlinkManager
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

params = open('NeVOmics_params.txt','w')
params.close()

# default values
with open("NeVOmics_params.txt", "a") as f:
    f.write('edgecolor=Colormap8\n')
    f.write('networkcolor=GreBlaRed\n')
    f.write('usertext=LogFC\n')
    f.write('uniquecolor=#0000FF\n')
    f.close()


if os.path.isfile('NeVOmics_img/KEGG_Organisms.txt'):
    han = open('NeVOmics_img/KEGG_Organisms.txt', 'r')
    dict_org = {}
    for line in han:
        line = line.rstrip()
        if re.search('^#', line):
            pass
        else:
            separados = line.split('\t')
            dict_org[' '+separados[2]] = [separados[1], separados[0]]
    han.close()
    orgs = list(dict_org.keys())
    orgs = sorted(orgs, key=str.upper)

else:
    orgs = []



opciones = ["1.-  Gene Ontology Enrichment",
            "2.-  KEGG Pathways Enrichment",
            "3.-  KEGG Blast Pathways Enrichment"]

print('NeVOmics History\n')

root = Tk()
#root.tk.call('tk', 'scaling', 1.5)
root.title("NeVOmics")
root.iconbitmap(r'NeVOmics_img/icon_nevomics.ico')
root.geometry("760x630")
#root.resizable(0, 0) # fija las dimensiones de la ventana
root.configure(background='white')

columna0 = 0
columna1 = 1
columna2 = 2
columna3 = 3
columna4 = 4
columna5 = 5
columna6 = 6

fila0 = 0
fila1 = 1
fila2 = 2
fila3 = 3
fila5 = 5
fila7 = 7


actualizaciones = LabelFrame(root, text=" Updates, Downloads and Installations", borderwidth=2, font=("Arial", 8, "bold")) # 
actualizaciones.grid(column=1, columnspan = 11, row=0, sticky = W+E+S+N)
actualizaciones.configure(background='white')


#$$$$$$$$$$$$$$$$$$$$$$$$

def updateGO():
    url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    go_file = 'NeVOmics_img/go-basic.obo'
    with open(go_file, 'wb') as f:
        response = requests.get(url, stream=True)
        total_length = response.headers.get('content-length')
        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=8192):
                dl += len(data)
                f.write(data)
                done = int(40 * dl / total_length)
                sys.stdout.write("\rDownloading the ontology [%s%s] %s MB" % ('>' * done, ' ' * (40-done), round(dl/1000000,2)), )    
                sys.stdout.flush()
        f.close()
    print('\nFinished')
    gobasic = open('NeVOmics_img/go-basic.obo', 'r')
    for line in gobasic:
        if re.search('data-version: .*', line):
            go_version = re.search('data-version: .*', line).group()
            break
            gobasic.close()

    from datetime import datetime
    xx = datetime.now()
    hora = '('+'{}'.format(xx).split('.')[0].split(' ')[1]+')'

    infogo0 = Label(root, text=' Gene Ontology: '+go_version, font=("Arial", 7, "bold"), bg = 'white', compound = LEFT)
    infogo0.grid(column=1, row=20, sticky= W, columnspan=4)


upgo = Button(actualizaciones, text=' Gene Ontology (GO) v1.2', bg='#ffb75f', fg="blue", font=("Arial", 7, "bold"), borderwidth=0,
            command = updateGO, 
            activeforeground = 'orange', cursor="hand2")
upgo.grid(column = 0, row = 0, sticky= W+N)

###$$$$$$$$$$$$$$$$$$$$$$$$$

c00 = Label(actualizaciones, text=" ", bg = 'white').grid(column=1, row=0)


def updateKEGGorgs():
    print('Updating KEGG Organisms...')
    kegg_organismos = requests.get("http://rest.kegg.jp/list/organism")
    kegginfo = re.findall('Release .*', requests.get('http://rest.kegg.jp/info/hsa').content.decode())[0]
    kegg_orgs = kegg_organismos.text.rstrip()
    kegg_orgs = re.sub('"', '', kegg_orgs)
    kegg_orgs = re.sub("'", '', kegg_orgs)
    with open('NeVOmics_img/KEGG_Organisms.txt', 'w') as fq:
        for h in [i for i in kegg_orgs.split('\n')]:
            fq.write(h+'\n')
        fq.close()
    from datetime import datetime
    xx = datetime.now()
    hora = '('+'{}'.format(xx).split('.')[0].split(' ')[1]+')'
    with open('NeVOmics_img/KEGG_Organisms.txt', 'a') as fq:
        fq.write('#'+kegginfo)
        fq.close()
    han = open('NeVOmics_img/KEGG_Organisms.txt', 'r')
    dict_org = {}
    for line in han:
        line = line.rstrip()
        if re.search('^#', line):
            pass
        else:
            separados = line.split('\t')
            dict_org[' '+separados[2]] = [separados[1], separados[0]]
    han.close()
    orgs = list(dict_org.keys())
    orgs = sorted(orgs, key=str.upper)
    print('Finished')

    ###
    def callback1(event):
        webbrowser.open_new(r"https://www.kegg.jp/dbget-bin/www_bget?gn:"+dict_org[org_kegg.get()][1])

    nueve = ttk.Combobox(group_org, font=("Arial", 8),
                     cursor="hand2", textvariable = org_kegg, values = orgs, height=20, width=42)
    nueve.grid(column=1, row=1, sticky= W)
    nueve.current(0)

    vacia22 = Label(group_org, textvariable = org_kegg, bg = 'white', fg = 'blue', font=("Arial", 8),
    cursor="hand2")# file vacía
    vacia22.grid(column=1, columnspan=4, row=2, sticky= W)
    vacia22.bind("<Button-1>", callback1)

    nueve1 = ttk.Combobox(group_org1, textvariable = org_keggbl, font=("Arial", 8),
                     cursor="hand2", values = orgs, height=20, width=42)
    nueve1.grid(column=1, row=1, sticky= W)
    nueve1.current(0)

    infokegg1 = Label(root, text=' KEGG: '+kegginfo, font=("Arial", 7, "bold"), bg = 'white', compound = LEFT)
    infokegg1.grid(column=1, row=21, sticky= W, columnspan=4)





uporgkegg = Button(actualizaciones, text=' KEGG Organisms ', bg='#ffb75f', fg="blue", font=("Arial", 7, "bold"), borderwidth=0,
            command = updateKEGGorgs, 
            activeforeground = 'orange', cursor="hand2")
uporgkegg.grid(column = 2, row = 0, sticky= N)

c000 = Label(actualizaciones, text=" ", bg = 'white').grid(column=3, row=0)

def callback_R(event):
    webbrowser.open_new(r"https://cran.r-project.org/bin/windows/base/old/3.5.3/R-3.5.3-win.exe")
downR = Button(actualizaciones, text = '  Download R v3.5.3  ', bg='#ffb75f', font=("Arial", 7, "bold"), fg="blue", borderwidth=0, cursor="hand2",
activeforeground = 'orange')
downR.grid(column = 4, row = 0, sticky= N)
downR.bind("<Button-1>", callback_R)

c001 = Label(actualizaciones, text=" ", bg = 'white').grid(column=5, row=0)

def Rlibrary(event):
    #webbrowser.open_new(r"https://drive.google.com/file/d/1u9PkZ5UWcKeOaJGjhS-VYm0WIz0vd6cP/view?usp=sharing")
    #webbrowser.open_new(r"https://drive.google.com/uc?id=1u9PkZ5UWcKeOaJGjhS-VYm0WIz0vd6cP&export=download")
    webbrowser.open_new(r"https://drive.google.com/uc?export=download&confirm=wuyZ&id=1u9PkZ5UWcKeOaJGjhS-VYm0WIz0vd6cP")

downR = Button(actualizaciones, text = '  Download R Library v3.5.3  ', bg='#ffb75f', font=("Arial", 7, "bold"), fg="blue", borderwidth=0, cursor="hand2",
activeforeground = 'orange')
downR.grid(column = 6, row = 0, sticky= N)
downR.bind("<Button-1>", Rlibrary)


c001 = Label(actualizaciones, text=" ", bg = 'white').grid(column=7, row=0)


def callback_Bl(event):
    webbrowser.open_new(r"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-win64.exe")
downR = Button(actualizaciones, text = '  Download Blast v2.8.1  ', bg='#ffb75f', font=("Arial", 7, "bold"), fg="blue", borderwidth=0, cursor="hand2",
activeforeground = 'orange')
downR.grid(column = 8, row = 0, sticky= N)
downR.bind("<Button-1>", callback_Bl)


#---------------------------------
c001 = Label(actualizaciones, text="                        ", bg = 'white').grid(column=9, row=0)

def ayuda():
    newwin0 = Tk()
    newwin0.title("NeVOmics")
    newwin0.iconbitmap(r'NeVOmics_img/icon_nevomics.ico')
    newwin0.geometry("515x600")
    newwin0.configure(background='white')
    newwin0.resizable(0, 0)

    c00 = Label(newwin0, text="    ", bg = 'white')
    c00.grid(column=0, row=0)

    text2 = Text(newwin0, height=40, width=80, bg = 'white', bd = 0, font=('Arial', 8))
    scroll = Scrollbar(newwin0, command=text2.yview)
    text2.configure(yscrollcommand=scroll.set)

    text2.tag_configure('big', font=('Arial', 9, 'bold'))
    # titulos

    # contenidos
    text2.tag_configure('color', foreground='red', font=('Arial', 8, 'bold'))
    text2.tag_configure('color2', foreground='darkorange', font=('Arial', 8, 'bold'))

    text2.tag_configure('text', font=('Arial', 8))
    text2.tag_configure('text2', font=('Arial', 8, 'bold'))

    hyperlink = HyperlinkManager(text2)

    #------------------------------------------------------------
    text2.insert(END,'\n◼◼ Gene Ontology (GO) ◼◼\n\n', 'big')
    text2.insert(END, "If you want to update the GO at any time, click on button:\n", 'text')
    text2.insert(END, 'Gene Ontology (GO) v1.2\n\n', 'color2')


    text2.insert(END, "The GO is downloaded from:\n", 'text')
    def click1():
        webbrowser.open_new(r"http://geneontology.org/docs/download-ontology/")
    text2.insert(END, "http://geneontology.org/docs/download-ontology/", hyperlink.add(click1))

    text2.insert(END, "\n\nFor more information about GO visit:\n", 'text')

    def click2():
        webbrowser.open_new(r"http://geneontology.org/")
    text2.insert(END, "http://geneontology.org/\n", hyperlink.add(click2))

    #------------------------------------------------------------
    text2.insert(END,'\n\n◼◼ KEGG (Kyoto Encyclopedia of Genes and Genomes) Organisms ◼◼\n\n', 'big')

    text2.insert(END, "If you want to update the KEGG organisms list click on button:\n", 'text')
    text2.insert(END, 'KEGG Organisms\n\n', 'color2')


    text2.insert(END, 'This list belongs to organisms with complete genomes \ndeposited in the '\
                'KEGG database, and is downloaded from:\n', 'text')
    def click3():
        webbrowser.open_new(r"https://www.kegg.jp/kegg/catalog/org_list.html")
    text2.insert(END, "https://www.kegg.jp/kegg/catalog/org_list.html", hyperlink.add(click3))

    text2.insert(END, "\n\nFor more information about KEGG visit:\n", 'text')
    def click4():
        webbrowser.open_new(r"https://www.kegg.jp/kegg/")
    text2.insert(END, "https://www.kegg.jp/kegg/\n", hyperlink.add(click4))

    #------------------------------------------------------------

    text2.insert(END,'\n\n◼◼ Download and Install R v3.5.3 ◼◼\n\n', 'big')

    text2.insert(END, "If you want to install R v3.5.3 click on button:\n", 'text')
    text2.insert(END, 'Download R v3.5.3\n\n', 'color2')

    text2.insert(END, 'How to install R v3.5.3 on Windows 10?\n', 'text2')

    text2.insert(END, '  1.- Search Download R-3.5.3-win.exe to install\n', 'text')
    text2.insert(END, '  2.- Double-click on R-3.5.3-win.exe\n', 'text')
    text2.insert(END, '  3.- OK (Select the language [English])\n', 'text')
    text2.insert(END, '  4.- Next (Information about R)\n', 'text')
    text2.insert(END, '  5.- Next (Select destination location, [default C:\Program Files\R\R-3.5.3])\n', 'text')
    text2.insert(END, '  6.- Next (Select components [default options])\n', 'text')
    text2.insert(END, '  7.- Next (Startup options [default options])\n', 'text')
    text2.insert(END, '  8.- Next (Select additional tasks [default options])\n', 'text')
    text2.insert(END, '  9.- Finish (Setup has finished installing R v3.5.3 for Windows)\n\n', 'text')

    text2.insert(END, 'NeVOmics is compatible only with R v3.5.3\n', 'color')


    text2.insert(END, "\nFor more information about R for Windows visit:\n", 'text')
    def click5():
        webbrowser.open_new(r"https://cran.r-project.org/bin/windows/base/")
    text2.insert(END, "https://cran.r-project.org/bin/windows/base/\n", hyperlink.add(click5))

    #------------------------------------------------------------

    text2.insert(END,'\n\n◼◼ Download R Library ◼◼\n\n', 'big')

    text2.insert(END, "If you want to download R Library v3.5.3 click on button:\n", 'text')
    text2.insert(END, 'Download R Library v3.5.3\n\n', 'color2')

    text2.insert(END, 'This library was built with R packages compatible with NeVOmics\n', 'text2')

    text2.insert(END, '\nThis R Library is compatible only with R v3.5.3\n', 'color')

    #------------------------------------------------------------

    text2.insert(END,'\n\n◼◼ Download and Install Blast v2.8.1 ◼◼\n\n', 'big')

    text2.insert(END, "If you want to install Blast v2.8.1 click on button:\n", 'text')
    text2.insert(END, 'Download Blast v2.8.1\n\n', 'color2')

    text2.insert(END, 'How to install Blast v2.8.1 on Windows 10?\n', 'text2')
    text2.insert(END, '  1.- Search Download ncbi-blast-2.8.1+-win64.exe to install\n', 'text')
    text2.insert(END, '  2.- Double-click on ncbi-blast-2.8.1+-win64.exe\n', 'text')
    text2.insert(END, '  3.- I gree (License agreement)\n', 'text')
    text2.insert(END, '  4.- Install (Choose install location [default C:\Program Files\\NCBI\\blast-2.8.1+]) \n', 'text')
    text2.insert(END, '  5.- Close (Installation complete)\n\n', 'text')

    text2.insert(END, 'The Blast parameters used by NeVOmics are:\n', 'text2')
    text2.insert(END, '-evalue = 1E-6, -outfmt = "6", -max_target_seqs = 50, -max_hsps = 50,\n', 'text')
    text2.insert(END, '-num_threads = all cores found on the computer.\n\n', 'text')
    text2.insert(END, 'The best hits (close to 1) are obtained using the following custom approach:\n', 'text')
    text2.insert(END, 'Qscore = mean of ((nident / qlen), (nident / length), ((length-gaps) / qlen))\n', 'text')
    text2.insert(END, 'where: nident means Number of identical matches,\n', 'text')
    text2.insert(END, 'qlen means Query sequence length, length means Alignment length,\n', 'text')
    text2.insert(END, 'and gaps means Total number of gaps.\n', 'text')

    text2.insert(END, "\nFor more information about BLAST+ executables\n", 'text')
    def click6():
        webbrowser.open_new(r"https://www.ncbi.nlm.nih.gov/books/NBK279671/")
    text2.insert(END, "https://www.ncbi.nlm.nih.gov/books/NBK279671/\n", hyperlink.add(click6))

    text2.grid(column=1, row=0)
    scroll.grid(column=1, row=0, sticky = E+S+N)
    text2.configure(state=DISABLED)
    newwin0.mainloop()


duda = Button(actualizaciones, text="  Help  ", borderwidth=0, activeforeground = 'blue',
             font=("Arial", 7, "bold"), fg="white", cursor="hand2", bg = 'darkblue', command = ayuda)
duda.grid(column=10, row=0, sticky = N)
#duda.bind("<Button-1>", aclaraciones)











#________________________________________________________________
def ORG_GO():
    if "KEGG" in opciones[analysis.get()]:
        if os.path.isfile('NeVOmics_img/KEGG_Organisms.txt'):
            pass
        else:
            messagebox.showinfo('Status', 'You must update the list of deposited organisms from the KEGG database')
    if "Ontology" in opciones[analysis.get()]:
        if os.path.isfile('NeVOmics_img/go-basic.obo'):
            pass
        else:
            messagebox.showinfo('Status', 'You must update or download the Genetic Ontology')



c0 = Label(root, text="  ", bg = 'white')
c0.grid(column=columna0, row=fila1)

uno = Label(root, text="1. Functional Analysis", font=("Arial", 10, "bold"), bg = 'white')
uno.grid(column=columna1, columnspan=4,  # este "columnspan" debe ser igual que "columnspan" del radiobutton de abajo
         row=1, sticky= W+N)

colors = {0:'tan', 1:'tan', 2:'tan'}
analysis = IntVar()#
pos = {0:2,1:8,2:12}
for i, option in enumerate(opciones):
    radio = Radiobutton(root, text=option, font=("Arial", 8, "bold"), bg = colors[i],
                        indicatoron=0, selectcolor='black', borderwidth=0, command = ORG_GO,
                        fg = 'white', cursor="hand2", variable=analysis, value=i)
    radio.grid(column=columna1, columnspan=3, row=pos[i], sticky= W+E)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_aspect = LabelFrame(root, text = "Aspects", font=("Arial", 8))
group_aspect.grid(column=columna1, row=3, rowspan=3, sticky= W+E+N+S)
group_aspect.configure(background='white')

group_adj = LabelFrame(root, text="P-value adjusted", font=("Arial", 8))
group_adj.grid(column=columna2, row=3, rowspan=3, sticky= W+E+N+S) 
group_adj.configure(background='white')

group_plots = LabelFrame(root, text = "Plots", font=("Arial", 8))
group_plots.grid(column=columna3, row=3, rowspan=3, sticky= W+E+N+S)
group_plots.configure(background='white')

#--------------------------------------------------------------
uno = Label(group_aspect, text = 'Biological Process: ', font=("Arial", 8, "bold"), bg = 'white')
uno.grid(column = 1, row = 1, sticky = E)


cuatro = Label(group_adj, text= 'FDR (%): ', font=("Arial", 8), bg = 'white')
cuatro.grid(column= 1, row=1, sticky= S+W)

bpz = StringVar()
cinco = ttk.Combobox(group_adj, textvariable = bpz, font=("Arial", 8),
                     values = valores, width=4)
cinco.grid(column=2, row=1, sticky= S+W)
cinco.current(38)

seis = Label(group_adj, text=" ", font=("Arial", 8), bg = 'white')
seis.grid(column=3, row=1, sticky= S+W)

bpgr = IntVar()
siete = Checkbutton(group_plots, text='', font=("Arial", 5), bg = 'white',
                    cursor="hand2", variable=bpgr)
siete.grid(column = 1, row = 1, sticky = N)

###--------------------------------------------------------------

uno1 = Label(group_aspect, text = 'Molecular Function: ', font=("Arial", 8, "bold"), bg = 'white')
uno1.grid(column = 1, row = 2, sticky = E)


cuatro1 = Label(group_adj, text= 'FDR (%): ', font=("Arial", 8), bg = 'white')
cuatro1.grid(column=1, row=2, sticky= S+W)

mfz = StringVar()
cinco1 = ttk.Combobox(group_adj, textvariable = mfz, font=("Arial", 8),
                     values = valores, width=4)
cinco1.grid(column=2, row=2, sticky= S+W)
cinco1.current(38)

seis1 = Label(group_adj, text=" ", font=("Arial", 8), bg = 'white')
seis1.grid(column=3, row=2, sticky= S+W)

mfgr = IntVar()
siete1 = Checkbutton(group_plots, text='', font=("Arial", 5),borderwidth=2, bg = 'white',
                    cursor="hand2", variable=mfgr)
siete1.grid(column = 1, row = 2, sticky = N)

###--------------------------------------------------------------

uno2 = Label(group_aspect, text = 'Cellular Component: ', font=("Arial", 8, "bold"), bg = 'white')
uno2.grid(column = 1, row = 3, sticky = E)

cuatro2 = Label(group_adj, text= 'FDR (%): ', font=("Arial", 8), bg = 'white')
cuatro2.grid(column=1, row=3, sticky= S+W)

ccz = StringVar()
cinco2 = ttk.Combobox(group_adj, textvariable = ccz, font=("Arial", 8),
                     values = valores, width=4)
cinco2.grid(column=2, row=3, sticky= S+W)
cinco2.current(38)

seis2 = Label(group_adj, text=" ", font=("Arial", 8), bg = 'white')
seis2.grid(column=3, row=3, sticky= S+W)

ccgr = IntVar()
siete2 = Checkbutton(group_plots, text='', font=("Arial", 5),borderwidth=2, bg = 'white', 
                     cursor="hand2", variable=ccgr)
siete2.grid(column = 1, row = 3, sticky = N)

#=====================================================================

anotacion = LabelFrame(root, text="Default GO Annotation & Complete GO Annotation", font=("Arial", 8)) # 
anotacion.grid(column=columna1, columnspan = 3, row=6, sticky= W+E+N+S)
anotacion.configure(background='white')

check_annotation = Checkbutton(anotacion, text= "    UniProtKB    ", font=("Arial", 7, "bold"), borderwidth=0, 
                               state=DISABLED, bg = '#90ee90', fg = 'black')
check_annotation.grid(column=1, row=1)


vacia = Label(anotacion, text="  ", bg = 'white')# file vacía
vacia.grid(column=2, row=1)

# aviso al usuario de que la descarga de la anotacion completa demora
def advertencia():
    if aceptar.get() == 1:
        messagebox.showinfo('Status',
                            'The full annotation will be downloaded, this may '\
                            'take some time, from 10 min to 2 h.\n\n'\
                            'Once the **Complete Annotation** is downloaded, it will no'\
                            ' longer be downloaded for later analysis, if you later '\
                            'wish to update the **Complete Annotation** (for example every month),'\
                            ' simply delete '\
                            'the file labeled with: **Complete_Annotation_...**')
    if aceptar.get() == 0:
        pass

aceptar = IntVar()
check = Checkbutton(anotacion, text="    QuickGO    ", bg = '#90ee90', fg = 'black',
                    variable = aceptar, cursor="hand2", borderwidth=0, font=("Arial", 7, "bold"),
                   command = advertencia)
check.grid(column= 3, row= 1)

vacia2 = Label(anotacion, text="   ", bg = 'white', fg = 'white')# file vacía
vacia2.grid(column=4, row=1)

def aclaracion(event):
    webbrowser.open_new(r"https://www.uniprot.org/help/complete_go_annotation")
acla = Label(anotacion, text="?",
             font=("Arial", 9, "bold"), fg="blue", cursor="hand2", bg = 'white')
acla.grid(column=5, row=1, sticky = W)
acla.bind("<Button-1>", aclaracion)

##########################################################################

vacia1 = LabelFrame(root, text="", bd = 0) # 
vacia1.grid(column=1, columnspan = 3, row=7, sticky = W+E+S+N)
vacia1.configure(background='white')

vacia = Label(vacia1, text="    ", bg = 'white')# file vacía
vacia.grid(column=0, row=0)

#############################################################################
fila8 = 9
group_kegg = LabelFrame(root, text = "Classification", font=("Arial", 8))
group_kegg.grid(column=columna1, row=fila8, sticky= W+E+N+S)
group_kegg.configure(background='white')

group_adj1 = LabelFrame(root, text="P-value adjusted", font=("Arial", 8))
group_adj1.grid(column=columna2, row=fila8, sticky= W+E+N+S)
group_adj1.configure(background='white')

group_plots1 = LabelFrame(root, text = "Plots", font=("Arial", 8))
group_plots1.grid(column=columna3, row=fila8, sticky= W+E+N+S)
group_plots1.configure(background='white')


uno33 = Label(group_kegg, text = '        KEGG Pathways: ', font=("Arial", 8, "bold"), bg = 'white')
uno33.grid(column = 1, row = 1, sticky = E)

cuatro22 = Label(group_adj1, text= 'FDR (%): ', font=("Arial", 8), bg = 'white')
cuatro22.grid(column=1, row=1, sticky= S+W)

keggz = StringVar()
cinco22 = ttk.Combobox(group_adj1, textvariable = keggz, font=("Arial", 8),
                     values = valores, width=4)
cinco22.grid(column=2, row=1, sticky= S+W)
cinco22.current(38)

seis22 = Label(group_adj1, text=" ", font=("Arial", 8), bg = 'white')
seis22.grid(column=3, row=1, sticky= S+W)

kegggr = IntVar()
siete22 = Checkbutton(group_plots1, text='', font=("Arial", 5),borderwidth=4, bg = 'white', 
                     cursor="hand2", variable=kegggr)
siete22.grid(column = 1, row = 1, sticky = N)


fila9 = 10
group_org = LabelFrame(root, text = "Specific organism and KEGG Information*", borderwidth=2, bg = 'white', font=("Arial", 8))
group_org.grid(column=columna1, columnspan=3, row=fila9, rowspan=2, sticky= W+E+N+S)
group_org.configure(background='white')




org_kegg = StringVar()
if orgs == []:
    nueve = ttk.Combobox(group_org, font=("Arial", 8),
                     cursor="hand2", values = orgs, height=20, width=42)
    nueve.grid(column=1, row=1, sticky= W)
else:
    def callback1(event):
        webbrowser.open_new(r"https://www.kegg.jp/dbget-bin/www_bget?gn:"+dict_org[org_kegg.get()][1])

    vacia22 = Label(group_org, textvariable = org_kegg, bg = 'white', font=("Arial", 8),
        fg = 'blue', cursor="hand2")# file vacía
    vacia22.grid(column=1, row=2, sticky= W)
    vacia22.bind("<Button-1>", callback1)





    nueve = ttk.Combobox(group_org, textvariable = org_kegg, font=("Arial", 8),
                     cursor="hand2", values = orgs, height=20, width=42)
    nueve.grid(column=1, row=1, sticky= W)
    
    nueve.current(0)


vacia22 = Label(group_org, text= '', bg = 'white', fg = 'blue', font=("Arial", 8))# file vacía
vacia22.grid(column=1, row=2, sticky= W)

#############################################################################
fila13 = 13
group_kegg2 = LabelFrame(root, text = "Classification", font=("Arial", 8))
group_kegg2.grid(column=columna1, row=fila13, sticky= W+E+N+S)
group_kegg2.configure(background='white')

group_adj2 = LabelFrame(root, text="P-value adjusted", font=("Arial", 8))
group_adj2.grid(column=columna2, row=fila13, sticky= W+E+N+S)
group_adj2.configure(background='white')

group_plots2 = LabelFrame(root, text = "Plots", font=("Arial", 8))
group_plots2.grid(column=columna3, row=fila13, sticky= W+E+N+S)
group_plots2.configure(background='white')

uno333 = Label(group_kegg2, text = '        KEGG Pathways: ', font=("Arial", 8, "bold"), bg = 'white')
uno333.grid(column = 1, row = 1, sticky = E)

cuatro222 = Label(group_adj2, text= 'FDR (%): ', font=("Arial", 8), bg = 'white')
cuatro222.grid(column=1, row=1, sticky= S+W)

keggblz = StringVar()
cinco222 = ttk.Combobox(group_adj2, textvariable = keggblz, font=("Arial", 8),
                     values = valores, width=4)
cinco222.grid(column=2, row=1, sticky= S+W)
cinco222.current(38)

seis222 = Label(group_adj2, text=" ", font=("Arial", 8), bg = 'white')
seis222.grid(column=3, row=1, sticky= S+W)

keggblgr = IntVar()
siete222 = Checkbutton(group_plots2, text='', font=("Arial", 5),borderwidth=2, bg = 'white',
                     cursor="hand2", variable=keggblgr)
siete222.grid(column = 1, row = 1, sticky = N)


fila14 = 14
group_org1 = LabelFrame(root, text = "Specific organism", bg = 'white', borderwidth=2, font=("Arial", 8))
group_org1.grid(column=columna1, columnspan=3, row=fila14, sticky= W+E+N)
group_org1.configure(background='white')


org_keggbl = StringVar()
if orgs == []:
    nueve1 = ttk.Combobox(group_org1, font=("Arial", 8),
                     cursor="hand2", values = orgs, height=20, width=42)
    nueve1.grid(column=1, row=1, sticky= W)
else:
    nueve1 = ttk.Combobox(group_org1, textvariable = org_keggbl, font=("Arial", 8),
                     cursor="hand2", values = orgs, height=20, width=42)
    nueve1.grid(column=1, row=1, sticky= W)
    nueve1.current(0)



fila15 = 15
tool = LabelFrame(root, text = "Tool", font=("Arial", 8))
tool.grid(column=columna1, row=fila15, sticky= W+E+N+S)
tool.configure(background='white')

method = LabelFrame(root, text = "Method", font=("Arial", 8))
method.grid(column=columna2, columnspan = 2, row=fila15, sticky= W+E+N+S)
method.configure(background='white')

diez = Label(tool, text= '                Local Blast: ', font=("Arial", 8, "bold"), bg = 'white')
diez.grid(column=1, row=1, sticky= E)

keggblmethod = IntVar()
tipos = ['Blastp', 'Blastx']
mets = {0:1,1:2}
for i, tips in enumerate(tipos):
    once = Radiobutton(method, text='      '+tips+'     ', font=("Arial", 8, "bold"), cursor="hand2", borderwidth=0,
                       indicatoron=0, selectcolor='salmon', bg = '#ffeae6',  fg = 'black',
                       variable=keggblmethod, value=i)
    once.grid(column=mets[i], row=1)

fila16 = 16
review1 = LabelFrame(root, text = "Proteome", font=("Arial", 8))
review1.grid(column=columna1, row=fila16, rowspan=2, sticky= W+E+N+S)
review1.configure(background='white')

review2 = LabelFrame(root, text = "Entries", font=("Arial", 8))
review2.grid(column=columna2, columnspan =2, row=fila16, rowspan=2, sticky= W+E+N+S)
review2.configure(background='white')

review3 = Label(review1, text= '                  UniProtKB: \n', font=("Arial", 8, "bold"), bg = 'white',)
review3.grid(column=1, row=1, sticky= E)

REVISADO = IntVar()
reviewed = [' Reviewed\n  (Swiss-Prot)', ' Unreviewed     \n(TrEMBL)']
metsREV = {0:1,1:2}
for i, tipss in enumerate(reviewed):
    review4 = Radiobutton(review2, text=' '+tipss+' ', font=("Arial", 7, "bold"), cursor="hand2", borderwidth=0,
                       indicatoron=0, selectcolor='salmon', bg = '#ffeae6',  fg = 'black',
                       variable=REVISADO, value=i)
    review4.grid(column=metsREV[i], row=1, sticky= W)


###############################################################
archivo=StringVar()


def inputfile():
    seleccionado = ''.join(re.findall('[A-Z].*[a-z]', opciones[analysis.get()]))
    if seleccionado == 'KEGG Blast Pathways Enrichment':
        file_path = ''
        while file_path == '':
            file_path = filedialog.askopenfilename()
            if (file_path.split('.')[-1] == 'fasta') or (file_path.split('.')[-1] == 'fa'):
                archivo.set(file_path.split('/')[-1])# imprime el nombre del archivo seleccionado
                root.update()
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
                messagebox.showinfo('Status', 'For '+seleccionado+' analysis,'\
                                    'it must be a file in fasta format and with extension .fasta or .fa\nExample:\n'\
                                    '>Seq1\nHGIKPVISTQLLLNGSLAEEEIIIRSKNITDNTKTII\n'\
                                    '>Seq2\nWFGITNWLWYIRIFIMIVGGLIGLRIIFAVLSIVNRV')
                file_path = ''
    else:
        file_path = ''
        while file_path == '':
            file_path = filedialog.askopenfilename()
            if (file_path.split('.')[-1] == 'tsv') or (file_path.split('.')[-1] == 'txt'):
                archivo.set(file_path.split('/')[-1])# imprime el nombre del archivo seleccionado
                root.update()
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
                messagebox.showinfo('Status', 'For '+seleccionado+' analysis,'\
                                    'it must be a file separated by tabs (\\t) and with extension .txt or .tsv')
                file_path = ''
    # imprime el nombre del archivo seleccionado
    labebl = Label(ffile, textvariable= archivo, font=("Arial", 8),fg="red", bg = 'white')
    labebl.grid(column = 0, columnspan = 3, row = 3, sticky= W)

columna5 = 5

subir_archivo = Label(root, text="2. Select File", font=("Arial", 10, "bold"), bg = 'white')
subir_archivo.grid(column=columna5, row=1, sticky= W)

ffile = LabelFrame(root, text = "")
ffile.grid(column=columna5, columnspan = 3, row=2, rowspan=2, sticky= W+E+N+S)
ffile.configure(background='white')


subir_archivo1 = Radiobutton(ffile, text="           Upload :          ", bg = 'white',
                indicatoron=0, selectcolor='#3fbfff', fg="black", borderwidth = 0,
                font=("Arial", 8, "bold"), cursor="hand2", command=inputfile)
subir_archivo1.grid(column=0, row=2, sticky= W+N)

labebl = Label(ffile, text= '', font=("Arial", 8), fg="red", bg = 'white')
labebl.grid(column = 0, columnspan = 3, row = 3, sticky= W)

#..........................................................
#..........................................................
#..........................................................


#######   NETWORK
#######   NETWORK
#######   NETWORK
#######   NETWORK
#######   NETWORK
#######   NETWORK

subir_archivo = Label(root, text="3. Visualizations", font=("Arial", 10, "bold"), bg = 'white')
subir_archivo.grid(column=5, row=4, sticky= W)

plot_space = LabelFrame(root) # 
plot_space.grid(column=5, columnspan = 3, row=5, rowspan=7, sticky= W+E+N+S)
plot_space.configure(background='white')



image1 = tkinter.PhotoImage(file= 'NeVOmics_img/net.png')
netplot = IntVar()#--------------------------------------------------------------
siete4 = Checkbutton(plot_space, text="Networks :", image=image1, bg = 'white',
                    cursor="hand2", variable=netplot, compound='right', font=("Arial", 9, "bold"))
siete4.grid(column = 1, row = 1, sticky= W)




def color_palettes():
    
    newwin = Toplevel()
    newwin.title("NeVOmics")
    newwin.iconbitmap(r'NeVOmics_img/icon_nevomics.ico')
    newwin.configure(background='white')
    
    lab1 = Label(newwin, text="Choosing Colormaps", bg = 'white',
             font=("Arial", 20, "bold"),
             fg = 'grey')
    lab1.grid(column=0, columnspan=4, row=0, sticky= W+S+N)


    def callback(event):
        webbrowser.open_new(r"https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html")
    link = Label(newwin, text="https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html\n",
             font=("Arial", 8), fg="blue", cursor="hand2", bg = 'white')
    link.grid(column=0, row=1, sticky= W)
    link.bind("<Button-1>", callback)
    
    
    diverging = LabelFrame(newwin, text="Diverging colormaps") 
    diverging.grid(column=0, columnspan = 2, row=2, rowspan=18, sticky= W)
    diverging.configure(background='white')

    sequentials = LabelFrame(newwin, text = "Sequential colormaps")
    sequentials.grid(column=2, columnspan = 2, row=2, rowspan=18, sticky= W)
    sequentials.configure(background='white')

    uniforme = LabelFrame(newwin, text = "Uniform Sequential colormaps")
    uniforme.grid(column=4, columnspan = 2, row=2, rowspan=10, sticky= W)
    uniforme.configure(background='white')
    
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
    
        l1 = tkinter.Radiobutton(clases[j], image=image1, cursor="hand2", bg = 'white',
                   activebackground= 'white',activeforeground= 'white',
                         variable = color, value = i)
        
        l1.image = image1
        l1.grid(column = 1, row = fila)
        l2.append(l1)

        seis = Label(clases[j], text=j.split('/')[1].split('.')[0], font=("Arial", 8, "bold"), bg = 'white')
        seis.grid(column=0, row=fila, sticky= E)
        
    
        fila +=1
   
    fila_vacia = Label(newwin, text="   \n   ", bg = 'white')# file vacía
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
                bg="green", fg="white",font=("Arial", 10, "bold"),
                    command = ShowChoice) # newwin.destroy
    boton.grid(column=5, row=13, sticky= W)
    
    newwin.mainloop()


textonet2 = Label(plot_space, text="Node colors if there are values:", font=("Arial", 8), bg = 'white')
textonet2.grid(column=1, row=2, sticky= W)

b1 = Button(plot_space, text='GreBlaRed', bg='black', fg="white", font=("Arial", 7, "bold"), borderwidth=0,
            command=color_palettes, activeforeground = 'black', cursor="hand2")
b1.grid(column = 2, row = 2, sticky= E)

imagee = tkinter.PhotoImage(file= 'NeVOmics_img/GreBlaRed.png')
lbl = Label(plot_space, image=imagee,  bg = 'white')
lbl.grid(column=1, columnspan=2, row=3, sticky= W)

# colormap para edges




def edgeCOLORMAP():
    
    newwin = Toplevel()
    newwin.iconbitmap(r'NeVOmics_img/icon_nevomics.ico')
    newwin.configure(background='white')
    
    lab1 = Label(newwin, text="Choosing Colormaps", bg = 'white',
             font=("Arial", 20, "bold"),
             fg = 'grey')
    lab1.grid(column=0, columnspan=4, row=0, sticky= W+S+N)


    def callback(event):
        webbrowser.open_new(r"https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html")
    link = Label(newwin, text="https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html\n",
             font=("Arial", 8), fg="blue", cursor="hand2", bg = 'white')
    link.grid(column=0, row=1, sticky= W)
    link.bind("<Button-1>", callback)
    
    
    qualitative = LabelFrame(newwin, text="Qualitative colormaps", bg = 'white') 
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
    
        l1 = tkinter.Radiobutton(clases[j], image=image1, cursor="hand2", bg = 'white',
                   activebackground= 'white',activeforeground= 'white',
                         variable = color, value = i)
        
        l1.image = image1
        l1.grid(column = 1, row = fila)
        l2.append(l1)

        seis = Label(clases[j], text=j.split('/')[1].split('.')[0], font=("Arial", 8, "bold"), bg = 'white')
        seis.grid(column=0, row=fila, sticky= E)
        
    
        fila +=1
    
    fila_vacia = Label(newwin, text="   \n   ", bg = 'white')# file vacía
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
                bg="green", fg="white",font=("Arial", 12, "bold"),
                    command = ShowChoice) # newwin.destroy
    boton.grid(column=5, row=10, sticky= W)
    newwin.mainloop()

textonet3 = Label(plot_space, text="Terms and Edges color:                    ",
                  font=("Arial", 8), bg = 'white')
textonet3.grid(column=1, row=4, sticky= W)

colmap = Button(plot_space, text='Colormap8', bg="#000000", fg="white", font=("Arial", 7, "bold"), borderwidth=0,
                command= edgeCOLORMAP, cursor="hand2", activeforeground = 'black')
colmap.grid(column = 2, row = 4, sticky= E)

imagee2 = tkinter.PhotoImage(file= 'NeVOmics_img/Colormap8.png')
lbl2 = Label(plot_space, image=imagee2, bg = 'white')
lbl2.grid(column=1, columnspan=2, row=5, sticky= W)    


### color si no hay valores

def color_button():
    col = tkinter.colorchooser.askcolor(parent=root)
    sinback.configure(bg=col[1])
    sinback.configure(text='   '+col[1]+'   ')
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

textonet4 = Label(plot_space, text="Node colors if there are no values:", font=("Arial", 8), bg = 'white')
textonet4.grid(column=1,  row=6, sticky= W)                
                
sinback = Button(plot_space, text='   '+'#0000ff'+'   ', font=("Arial", 7, "bold"), borderwidth=0,
                         bg= '#0000ff', fg = 'white', command=color_button, cursor="hand2")
sinback.grid(column=2, row=6, sticky= E)


# texto

def USERTEXT():
    newwin = Toplevel()
    newwin.iconbitmap(r'NeVOmics_img/icon_nevomics.ico')
    newwin.configure(background='white')
    newwin.geometry("250x80")
    
    e = Entry(newwin, bd =3, width = 30) # para introducir titulo de la barra
    e.grid(column = 0, columnspan=4, row = 0, sticky= W)
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
        changetitle.configure(text='    '+e.get()+'    ', bg="white", font=("Arial", 7))    
        newwin.destroy()
    #-------------------------------------------------------------------------------
    
        
    boton = Button(newwin, text="OK", cursor="hand2",
                activebackground= 'black',activeforeground= 'black',
                bg="green", fg="white",font=("Arial", 10, "bold"),
                    command = ShowChoice) # newwin.destroy
    boton.grid(column=0, row=1, sticky= E)
    newwin.mainloop()
    
textonet5 = Label(plot_space, text="Colormap title (default: LogFC):   ", font=("Arial", 8), bg = 'white')
textonet5.grid(column=1, row=7, sticky= W)    

changetitle = Button(plot_space, text="     Title    ", bg="white", fg="black", font=("Arial", 7, "bold"),
                     activebackground = 'yellow', cursor="hand2", activeforeground = 'black', command=USERTEXT)
changetitle.grid(column = 2, row = 7, sticky= E)


#######   CIRCOS
#######   CIRCOS
#######   CIRCOS
#######   CIRCOS
#######   CIRCOS

plot_space2 = LabelFrame(root) # 
plot_space2.grid(column=5, columnspan = 3, row=12, rowspan=6, sticky= W+E+S+N)
plot_space2.configure(background='white')

nota = Label(plot_space2, text="[ To use this application you need to install R ]",
             font=("Arial", 8), fg="red", bg = 'white')
nota.grid(column=1, row=1, columnspan = 3, sticky= W)

def callback(event):
    webbrowser.open_new(r"https://github.com/bioinfproject/bioinfo")
link = Label(plot_space2, text="Check the documentation here.",
             font=("Arial", 7), fg="blue", cursor="hand2", bg = 'white')
link.grid(column=1, row=2, sticky= W)
link.bind("<Button-1>", callback)

nota2 = Label(plot_space2, text="[ The same colors of the networks will be used,\nyou can modify the parameters above ]",
              font=("Arial", 8), bg = 'white')
nota2.grid(column=1, row=3, columnspan = 3, sticky= W)

image2 = tkinter.PhotoImage(file= 'NeVOmics_img/circo.png')
circoplot = IntVar()#--------------------------------------------------------------
siete5 = Checkbutton(plot_space2, text="Circos :", image=image2,font=("Arial", 9, "bold"), bg = 'white',
                    cursor="hand2", variable=circoplot, compound='right')
siete5.grid(column = 1, row = 4, sticky= W)


### botones para localizar r
def r_exe_loc():
    if os.path.exists('NeVOmics_locRexe.txt'):
        rexe = open('NeVOmics_locRexe.txt', 'r')
        rexe = rexe.read()
        messagebox.showinfo('Status',
                            'The location of R.exe program is already defined in file:\nNeVOmics_locRexe.txt with location in:\n\n'+rexe+'\n\n\
        You can now run the analysis')
        pass
    else:
        R_exe = ''
        while R_exe == '':
            R_exe = filedialog.askopenfilename()
            if R_exe.split('/')[-1] == 'R.exe':
                f= open('NeVOmics_locRexe.txt','w')
                f.write(R_exe)
                f.close()
                #lbl = Label(root, text= R_exe)
                #lbl.grid(column=8, columnspan=5, row=10, sticky= N)
            else:
                messagebox.showinfo('Status',
                                    'It must be R.exe program\ne.g. C:/Users/home/Documents/R-3.5.3/bin/R.exe')
               
                R_exe = ''

textonet22 = Label(plot_space2, text="Locate the R.exe program:", font=("Arial", 8), bg = 'white')
textonet22.grid(column=1, row=5, sticky= W)

boton22 = Button(plot_space2, text=" R.exe       ", bg="#000000", fg = 'white', borderwidth=0,
                font=("Arial", 7, "bold"), activeforeground= 'black', 
                cursor="hand2", command = r_exe_loc)
boton22.grid(column=3, row=5, sticky= E)

# localizar la libreria

def r_lib_loc():
    if os.path.exists('NeVOmics_locRlib.txt'):
        rlib = open('NeVOmics_locRlib.txt', 'r')
        rlib = rlib.read()
        messagebox.showinfo('Status',
                            'The location of Rlibrary_NeVOmics directory is already defined in file:\nNeVOmics_locRlib.txt with location in:\n\n'+rlib+'\n\n\
 You can now run the analysis')
        pass
    else:
        R_lib = ''
        while R_lib == '':
            R_lib = filedialog.askdirectory()
            if R_lib.split('/')[-1] == 'Rlibrary_NeVOmics':
                ff= open('NeVOmics_locRlib.txt','w')
                ff.write(R_lib)
                ff.close()
                #lbl = Label(root, text= R_exe)
                #lbl.grid(column=8, columnspan=5, row=10, sticky= N)
            else:
                messagebox.showinfo('Status',
                                    'It must be Rlibrary_NeVOmics directory\ne.g. C:/Users/home/Downloads/Rlibrary_NeVOmics')
               
                R_lib = ''

textonet33 = Label(plot_space2, text="Locate the Rlibrary_NeVOmics directory:", font=("Arial", 8), bg = 'white')
textonet33.grid(column=1,columnspan=2, row=6, sticky= W)
                
boton33 = Button(plot_space2, text=" R Library ", bg="#000000", fg="white", borderwidth=0,
                font=("Arial", 7, "bold"),
                activebackground= 'yellow',activeforeground= 'black',
                cursor="hand2", command = r_lib_loc)
boton33.grid(column=3, row=6, sticky= E)


## tipo de identificadores


label = Label(root, text="4. Node label", font=("Arial", 10, "bold"), bg = 'white')
label.grid(column=9, row=1, sticky= W)



group_label5 = LabelFrame(root, text = "Identifier", font=("Arial", 8))
group_label5.grid(column=9, columnspan = 2, row=2, rowspan=2, sticky= W+E+N+S)
group_label5.configure(background='white')


labevacio = Label(group_label5, text=" ", font=("Arial", 10, "bold"), bg = 'white')
labevacio.grid(column=0, row=1, sticky= W)
labevacio = Label(group_label5, text=" ", font=("Arial", 10, "bold"), bg = 'white')
labevacio.grid(column=3, row=1, sticky= W)

node_lab = IntVar()
etiquetas = ['Gene Name', 'UniProt ID']
mets = {0:1,1:4}

for i, tips in enumerate(etiquetas):
    radio7 = Radiobutton(group_label5, text = '  '+tips+'  ', font=("Arial", 8, "bold"), cursor="hand2", borderwidth=0,
                         indicatoron=0, selectcolor='salmon', bg = '#ffeae6', fg = 'black',
                         variable = node_lab, value=i)
    radio7.grid(column=mets[i], row=1, sticky= W+E)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
varFolder=StringVar()
running = StringVar()
running.set('RUN')
root.update()

def parameters():

    running.set('RUNING...')
    botonrun.config(background= 'red')
    root.update()
    
    boxplots = (1 in [bpgr.get(), mfgr.get(), ccgr.get(), kegggr.get(), keggblgr.get()])
    tipoplots = (1 in [netplot.get(), circoplot.get()])
    
    if boxplots == True and tipoplots == False:
        #print('seleccionar al menos uno')
        messagebox.showinfo('Status',
                            'You have activated at least one Plots box, therefore it is required to choose at '\
                            'least one type of graphic: Networks or Chords.')
    elif boxplots == True and tipoplots == True:
        #print('selecionaron box plot, y también seleccionarion un tipo')
        # esta opción se para correr el programa con al menos un tipo de gráfico ***********
        if circoplot.get() == 1: # solo si eligen la aplicacion "Chords"
            locrexe = os.path.exists("NeVOmics_locRexe.txt")
            locrlib = os.path.exists("NeVOmics_locRlib.txt")
            if locrexe == True and locrlib == True:
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
                        #print(new_folder)

                    else:
                        n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
                        xoxo=str(n+1)
                        n=str(n)
                        old_folder=''.join(re.findall('job_NeVOmics_'+n,str(ls)))
                        new_folder=re.sub(n,xoxo,old_folder)
                        os.makedirs(new_folder,exist_ok=True)
                        #print(new_folder)
                
                    def open_folder():
                        os.system('start '+new_folder)
                        #subprocess.call('start '+new_folder, shell = True)

                    labelfolder.configure(text= 'Open:   '+new_folder+'  ', font=("Arial", 8, "bold"), command = open_folder,bg='wheat', borderwidth=0)
        
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
                            "keggorganism="+re.sub('^ ', '', org_kegg.get())+"\n"\
                            "#====================\n"\
                            "#**KEGG BLAST ENRICHMENT**\n"\
                            "keggblastfdr="+str(keggblz.get())+"\n"\
                            "keggblastplots="+str(keggblgr.get())+"\n"\
                            "keggmethodblast="+tipos[keggblmethod.get()]+"\n"\
                            "reviewed="+reviewed[REVISADO.get()]+"\n"\
                            "keggblastorganism="+re.sub('^ ', '', org_keggbl.get())+"\n"\
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
            
                    if  'Gene Ontology Enrichment' in analisis:
                        if os.path.isfile('NeVOmics_img/go-basic.obo'):
                            if ('fasta' or 'fa') in uno.split('/')[-1]:
                                messagebox.showinfo('Status', 'For Gene Ontology Enrichment analysis,'\
                                            'it must be a file separated by tabs (\\t) and with extension .txt or .tsv')
                            else:
                                print('________________________________________\n________________________________________\nRun: Gene Ontology Enrichment')
                                print('Job number:', new_folder.split('_')[-1])
                                print('New folder:', new_folder)
                                print(uno)
                                go_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_GO.py',
                                                                        new_folder+'/short_GO.py')
                                #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                                comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                                #os.system("start cmd /k cd "+comando+ " ^&^& python short_GO.py")

                                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                subprocess.call('cd '+comando+' & python short_GO.py', shell = True)
                                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                print('\n........................................\n')
                        else:
                            messagebox.showinfo('Status', 'You must update or download the Genetic Ontology')
                        
                    if  'KEGG Pathways Enrichment' in analisis:
                        print('________________________________________\n________________________________________\nRun: KEGG Pathways Enrichment')
                        print('Job number:', new_folder.split('_')[-1])
                        print('New folder:', new_folder)
                        print(uno)
                        kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG.py',
                                                                 new_folder+'/short_KEGG.py')
                        #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                        comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                        #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG.py")

                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        subprocess.call('cd '+comando+' & python short_KEGG.py', shell = True)
                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        print('\n........................................\n')
                        # >>>>>>>> fin
                    if 'KEGG Blast Pathways Enrichment' in analisis:
                        if ('tsv' or 'txt') in uno:
                            messagebox.showinfo('Status', 'For KEGG Blast Pathways Enrichment analysis,'\
                                    'it must be a file in fasta format and with extension .fasta or .fa\nExample:\n'\
                                    '>Seq1\nHGIKPVISTQLLLNGSLAEEEIIIRSKNITDNTKTII\n'\
                                    '>Seq2\nWFGITNWLWYIRIFIMIVGGLIGLRIIFAVLSIVNRV')
                        else:
                            print('________________________________________\n________________________________________\nRun: KEGG Blast Pathways Enrichment')
                            print('Job number:', new_folder.split('_')[-1])
                            print('New folder:', new_folder)
                            print(uno)
                            kegg_blast_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Blast.py',
                                                                    new_folder+'/short_KEGG_Blast.py')
                            #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                            comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                            #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG_Blast.py")

                            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            subprocess.call('cd '+comando+' & python short_KEGG_Blast.py', shell = True)
                            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            print('\n........................................\n')
                    
                else:
                    messagebox.showinfo('Status',
                                        'No file were selected.')
            else:
                messagebox.showinfo('Status',
                        'The location of R.exe programm or Rlibrary_NeVOmics directory is not yet defined '\
                                    'e.g.\nC:/Users/home/Documents/R-3.5.3/bin/R.exe\n'+\
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
                    #print(new_folder)
                else:
                    n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
                    xoxo=str(n+1)
                    n=str(n)
                    old_folder=''.join(re.findall('job_NeVOmics_'+n,str(ls)))
                    new_folder=re.sub(n,xoxo,old_folder)
                    os.makedirs(new_folder,exist_ok=True)
                    #print(new_folder)
                
                def open_folder():
                    os.system('start '+new_folder)
                    #subprocess.call('start '+new_folder, shell = True)

                labelfolder.configure(text= 'Open:   '+new_folder+'  ', font=("Arial", 8, "bold"), command = open_folder,bg='wheat', borderwidth=0)
        
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
                        "keggorganism="+re.sub('^ ', '', org_kegg.get())+"\n"\
                        "#====================\n"\
                        "#**KEGG BLAST ENRICHMENT**\n"\
                        "keggblastfdr="+str(keggblz.get())+"\n"\
                        "keggblastplots="+str(keggblgr.get())+"\n"\
                        "keggmethodblast="+tipos[keggblmethod.get()]+"\n"\
                        "reviewed="+reviewed[REVISADO.get()]+"\n"\
                        "keggblastorganism="+re.sub('^ ', '', org_keggbl.get())+"\n"\
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
            
                if  'Gene Ontology Enrichment' in analisis:
                    if os.path.isfile('NeVOmics_img/go-basic.obo'):
                        if ('fasta' or 'fa') in uno.split('/')[-1]:
                            messagebox.showinfo('Status', 'For Gene Ontology Enrichment analysis,'\
                                        'it must be a file separated by tabs (\\t) and with extension .txt or .tsv')
                        else:
                            print('________________________________________\n________________________________________\nRun: Gene Ontology Enrichment')
                            print('Job number:', new_folder.split('_')[-1])
                            print('New folder:', new_folder)
                            print(uno)
                            go_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_GO.py',
                                                                    new_folder+'/short_GO.py')
                            #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                            comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                            #os.system("start cmd /k cd "+comando+ " ^&^& python short_GO.py")

                            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            subprocess.call('cd '+comando+' & python short_GO.py', shell = True)
                            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            print('\n........................................\n')
                    else:
                        messagebox.showinfo('Status', 'You must update or download the Genetic Ontology')
                    

                if  'KEGG Pathways Enrichment' in analisis:
                    print('________________________________________\n________________________________________\nRun: KEGG Pathways Enrichment')
                    print('Job number:', new_folder.split('_')[-1])
                    print('New folder:', new_folder)
                    print(uno)
                    kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG.py',
                                                             new_folder+'/short_KEGG.py')
                    #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                    comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                    #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG.py")

                    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    subprocess.call('cd '+comando+' & python short_KEGG.py', shell = True)
                    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    print('\n........................................\n')
                    # >>>>>>>> fin
                if  'KEGG Blast Pathways Enrichment' in analisis:
                    if ('tsv' or 'txt') in uno:
                        messagebox.showinfo('Status', 'For KEGG Blast Pathways Enrichment analysis,'\
                                    'it must be a file in fasta format and with extension .fasta or .fa\nExample:\n'\
                                    '>Seq1\nHGIKPVISTQLLLNGSLAEEEIIIRSKNITDNTKTII\n'\
                                    '>Seq2\nWFGITNWLWYIRIFIMIVGGLIGLRIIFAVLSIVNRV')
                    else:
                        print('________________________________________\n________________________________________\nRun: KEGG Blast Pathways Enrichment')
                        print('Job number:', new_folder.split('_')[-1])
                        print('New folder:', new_folder)
                        print(uno)
                        kegg_blast_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Blast.py',
                                                                    new_folder+'/short_KEGG_Blast.py')
                        #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                        comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                        #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG_Blast.py")
                    
                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        subprocess.call('cd '+comando+' & python short_KEGG_Blast.py', shell = True)
                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        print('\n........................................\n')
                
            else:
                messagebox.showinfo('Status',
                                    'No file were selected.')
    elif boxplots == False and tipoplots == False:
        #print('no seleccionaron nada, no se harán gráficos')
        # esta opción es para correr el programa sin generar ningún tipo de gráficos ************
        #------------------------------>
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
            else:
                n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
                xoxo=str(n+1)
                n=str(n)
                old_folder=''.join(re.findall('job_NeVOmics_'+n,str(ls)))
                new_folder=re.sub(n,xoxo,old_folder)
                os.makedirs(new_folder,exist_ok=True)
            
            def open_folder():
                os.system('start '+new_folder)
                #subprocess.call('start '+new_folder, shell = True)

            labelfolder.configure(text= 'Open:   '+new_folder+'  ', font=("Arial", 8, "bold"), command = open_folder, 
            bg='wheat', borderwidth=0)
        
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
                    "keggorganism="+re.sub('^ ', '', org_kegg.get())+"\n"\
                    "#====================\n"\
                    "#**KEGG BLAST ENRICHMENT**\n"\
                    "keggblastfdr="+str(keggblz.get())+"\n"\
                    "keggblastplots="+str(keggblgr.get())+"\n"\
                    "keggmethodblast="+tipos[keggblmethod.get()]+"\n"\
                    "reviewed="+reviewed[REVISADO.get()]+"\n"\
                    "keggblastorganism="+re.sub('^ ', '', org_keggbl.get())+"\n"\
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
            
            if 'Gene Ontology Enrichment' in analisis:
                if os.path.isfile('NeVOmics_img/go-basic.obo'):
                    if ('fasta' or 'fa') in uno.split('/')[-1]:
                        messagebox.showinfo('Status', 'For Gene Ontology Enrichment analysis,'\
                                    'it must be a file separated by tabs (\\t) and with extension .txt or .tsv')
                    else:
                        print('________________________________________\n________________________________________\nRun: Gene Ontology Enrichment')
                        print('Job number:', new_folder.split('_')[-1])
                        print('New folder:', new_folder)
                        print(uno)

                        go_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_GO.py',
                                                                new_folder+'/short_GO.py')
                        #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                        comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                        #os.system("start cmd /k cd "+comando+ " ^&^& python short_GO.py")
                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        subprocess.call('cd '+comando+' & python short_GO.py', shell = True)
                        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        print('\n........................................\n')
                else:
                    messagebox.showinfo('Status', 'You must update or download the Genetic Ontology')


            if 'KEGG Pathways Enrichment' in analisis:
                print('________________________________________\n________________________________________\nRun: KEGG Pathways Enrichment')
                print('Job number:', new_folder.split('_')[-1])
                print('New folder:', new_folder)
                print(uno)

                kegg_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG.py',
                                                         new_folder+'/short_KEGG.py')
                #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG.py")
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                subprocess.call('cd '+comando+' & python short_KEGG.py', shell = True)
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                print('\n........................................\n')

                # >>>>>>>> fin
            if 'KEGG Blast Pathways Enrichment' in analisis:
                if ('tsv' or 'txt') in uno:
                    messagebox.showinfo('Status', 'For KEGG Blast Pathways Enrichment analysis,'\
                                    'it must be a file in fasta format and with extension .fasta or .fa\nExample:\n'\
                                    '>Seq1\nHGIKPVISTQLLLNGSLAEEEIIIRSKNITDNTKTII\n'\
                                    '>Seq2\nWFGITNWLWYIRIFIMIVGGLIGLRIIFAVLSIVNRV')
                else:
                    print('________________________________________\n________________________________________\nRun: KEGG Blast Pathways Enrichment')
                    print('Job number:', new_folder.split('_')[-1])
                    print('New folder:', new_folder)
                    print(uno)
                    

                    kegg_blast_script = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/short_KEGG_Blast.py',
                                                                new_folder+'/short_KEGG_Blast.py')
                    #print(re.sub('\\\\', '/', os.path.abspath(new_folder)))
                    comando = re.sub('\\\\', '/', os.path.abspath(new_folder))
                    #os.system("start cmd /k cd "+comando+ " ^&^& python short_KEGG_Blast.py")
                    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    subprocess.call('cd '+comando+' & python short_KEGG_Blast.py', shell = True)
                    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                    print('\n........................................\n')
            
            
        else:
            messagebox.showinfo('Status',
                                'No file were selected.')
    elif boxplots == False and tipoplots == True:
        #print('seleccionaron un tipo pero no seleccionaron un box plots')
        messagebox.showinfo('Status',
                            'You have activated at least one type of graphic: Networks or Chords, '\
                            'therefore it is required to choose at least one Plots box.')

    running.set('RUN')
    botonrun.config(background= 'green')
    root.update()

    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

botonrun = Button(root, textvariable= running, bg="green", fg="white", borderwidth=0,
                activebackground = 'red',
                font=("Arial", 15, "bold"), command = parameters, cursor="hand2")
botonrun.grid(column = 9, columnspan=2, row = 5, sticky= W+E+N+S)
####>>>>>

labelfolder = Button(root, text= '  New working folder  ', font=("Arial", 8, "bold"), cursor = "hand2", bg="oldlace",  borderwidth=0)
labelfolder.grid(column = 9, columnspan=2, row = 6, rowspan=2, sticky= W+E)



infokegg0 = Label(root, text=" Databases information:", font=("Arial", 7, "bold"), bg = 'white')
infokegg0.grid(column=1, row=19, sticky= W)




##############################
if os.path.isfile('NeVOmics_img/go-basic.obo'):
    from datetime import datetime
    xx = datetime.now()
    hora0 = '('+'{}'.format(xx).split('.')[0].split(' ')[1]+')'
    gobasic = open('NeVOmics_img/go-basic.obo', 'r')
    for line in gobasic:
        if re.search('data-version: .*', line):
            go_version = re.search('data-version: .*', line).group()
            
            break
            gobasic.close()
else:
    go_version = ''



infogo0 = Label(root, text=" Gene Ontology: "+go_version, font=("Arial", 7, "bold"), bg = 'white')
infogo0.grid(column=1, row=20, sticky= W+S, columnspan=4)




if os.path.isfile('NeVOmics_img/KEGG_Organisms.txt'):
    hann = reversed(open('NeVOmics_img/KEGG_Organisms.txt').readlines())
    for line in hann:
        line.rstrip()
        if re.search('^#', line):
            Kegginfoversion = re.sub('#', '', line)
            break
            hann.close()
else:
    Kegginfoversion = ''


infokegg1 = Label(root, text= ' KEGG: '+Kegginfoversion, font=("Arial", 7, "bold"), bg = 'white')
infokegg1.grid(column=1, row=21, sticky= W, columnspan=4)

def callback(event):
    webbrowser.open_new(r"https://doi.org/10.3390/genes9120569")
link = Label(root, text="  If you use NeVOmics, please cite."+''.join([' ']*40),compound = LEFT,
             font=("Arial", 7), fg="blue", cursor="hand2")
link.grid(column=1, row=22, sticky= W+N+S+E, columnspan=3)
link.bind("<Button-1>", callback)

####>>>>>
ejemplo = LabelFrame(root) #, text = "Example Network created by NeVOmics" 
ejemplo.grid(column=9, columnspan = 2, row=8, rowspan=11, sticky= W+E+S+N)
ejemplo.configure(background='white')

imagenet = tkinter.PhotoImage(file= 'NeVOmics_img/Network2.png')
labnet = Label(ejemplo, image=imagenet, bg="white")
        
labnet.image = imagenet
labnet.grid(column = 1, row =1, sticky= W+S)


root.mainloop()
