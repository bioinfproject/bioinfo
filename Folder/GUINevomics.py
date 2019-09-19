#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tkinter as tk
from tkinter import * 
from tkinter import filedialog
import webbrowser
import pandas as pd
from pandas.compat import StringIO
import requests
import urllib.request
from urllib.request import urlopen
import base64
from tkinter import ttk
import numpy as np
import os
from tkinter import messagebox
import subprocess
import tkinter


# In[ ]:





# In[2]:


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


# In[3]:


kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()
kegg_organism=pd.read_csv(StringIO(kegg_orgs),names=['T_number','Prefix','Organism','Group'],sep='\t')
kegg_organism = kegg_organism.sort_values(by ='Organism',ascending=True)


dict_org = {}
for index, row in kegg_organism.iterrows():
    dict_org[row.Organism] = [row.Prefix, row.T_number]


# In[21]:


params = open('NeVOmics_params.txt','w')
params.close()

opciones = [" Gene Ontology (GO) Enrichment                                  ",
            " KEGG Pathways Enrichment                                         ",
            " KEGG Blast Pathways Enrichment                                "]

root = Tk()
root.title("NeVOmics")
#           ancho , alto
root.geometry("750x570")
root.iconbitmap(r'icon_nevomics.ico')




lab0 = Label(root, text="       ")# columna vacía
lab0.grid(column=0, row=0)

lab1 = Label(root, text="NeVOmics",
             font=("Arial", 35, "bold"), underline=10,
             fg = 'green')
lab1.grid(column=1, columnspan=8, row=0, sticky= W+S+N)


uno = Label(root, text="1. Functional Analysis", font=("Arial", 15, "bold"))
uno.grid(column=1, columnspan=2, row=2, sticky= W+N)

colors = {0:'lightgray', 1:'lightgray', 2:'lightgray'}
analysis = IntVar()#--------------------------------------------------------------
pos = {0:3,1:8,2:12}
for i, option in enumerate(opciones):
    radio = Radiobutton(root, text=option, font=("Arial", 10, "bold"), bg = colors[i],
                        fg = 'black', cursor="hand2", variable=analysis, value=i)
    radio.grid(column=1, columnspan=3, row=pos[i], sticky= W+N)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_adj = LabelFrame(root, text="P-value adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_adj.grid(column=2, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_plots = LabelFrame(root, text = "Create") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_plots.grid(column=3, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_aspect = LabelFrame(root, text = "Aspects") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect.grid(column=1, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#group_pval = LabelFrame(root, text = "Not adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_pval.grid(column=2, row=4, rowspan=3)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=============================================================================
#oror = Label(root,text = " OR ", font = ("Arial", 11, "bold"))
#oror.grid(column=3, row=5, sticky= W+E+S)

#yyy = Label(root, text=" & ", font=("Arial", 12, "bold"))
#yyy.grid(column=5, row=5, sticky= W+E+S)
#=============================================================================

uno = Label(group_aspect, text = 'Biological Process: ', font=("Arial", 10, "bold"))
uno.grid(column = 1, row = 4, sticky = E)

#dos = Label(group_pval, text= 'P-value: ', font=("Arial", 10))
#dos.grid(column = 1, row = 4, sticky= S+W)

#bppval = StringVar()#--------------------------------------------------------------
#tres = ttk.Combobox(group_pval, textvariable = bppval, font=("Arial", 10),
#                     values = Pvalues_valores, width=5)
#tres.grid(column=2, row=4, sticky= S+W)
#tres.current(14)

#vacia = Label(group_pval, text = ' ', font = ("Arial", 10))
#vacia.grid(column=3, row=4, sticky= S+W)


cuatro = Label(group_adj, text= 'FDR: ', font=("Arial", 10))
cuatro.grid(column=1, row=4, sticky= S+W)

bpz = StringVar()#--------------------------------------------------------------
cinco = ttk.Combobox(group_adj, textvariable = bpz, font=("Arial", 10),
                     values = valores, width=4)
cinco.grid(column=2, row=4, sticky= S+W)
cinco.current(6)

seis = Label(group_adj, text=" % ", font=("Arial", 10))
seis.grid(column=3, row=4, sticky= S+W)

bpgr = IntVar()#--------------------------------------------------------------
siete = Checkbutton(group_plots, text=' Plots', font=("Arial", 10),
                    cursor="hand2", variable=bpgr)
siete.grid(column = 1, row = 4, sticky = W)

######

uno1 = Label(group_aspect, text = 'Molecular Function: ', font=("Arial", 10, "bold"))
uno1.grid(column = 1, row = 5, sticky = E)

#dos1 = Label(group_pval, text= 'P-value: ', font=("Arial", 10))
#dos1.grid(column = 1, row = 5, sticky= S+W)

#mfpval = StringVar()#--------------------------------------------------------------
#tres1 = ttk.Combobox(group_pval, textvariable = mfpval, font=("Arial", 10),
#                     values = Pvalues_valores, width=5)
#tres1.grid(column=2, row=5, sticky= S+W)
#tres1.current(14)

#vacia1 = Label(group_pval, text = ' ', font = ("Arial", 10))
#vacia1.grid(column=3, row=5, sticky= S+W)


cuatro1 = Label(group_adj, text= 'FDR: ', font=("Arial", 10))
cuatro1.grid(column=1, row=5, sticky= S+W)

mfz = StringVar()#--------------------------------------------------------------
cinco1 = ttk.Combobox(group_adj, textvariable = mfz, font=("Arial", 10),
                     values = valores, width=4)
cinco1.grid(column=2, row=5, sticky= S+W)
cinco1.current(6)

seis1 = Label(group_adj, text=" % ", font=("Arial", 10))
seis1.grid(column=3, row=5, sticky= S+W)

mfgr = IntVar()#--------------------------------------------------------------
siete1 = Checkbutton(group_plots, text=' Plots', font=("Arial", 10),
                    cursor="hand2", variable=mfgr)
siete1.grid(column = 1, row = 5, sticky = W)

######

uno2 = Label(group_aspect, text = 'Cellular Component: ', font=("Arial", 10, "bold"))
uno2.grid(column = 1, row = 6, sticky = E)

#dos2 = Label(group_pval, text= 'P-value: ', font=("Arial", 10))
#dos2.grid(column = 1, row = 6, sticky= S+W)

#ccpval = StringVar()#--------------------------------------------------------------
#tres2 = ttk.Combobox(group_pval, textvariable = ccpval, font=("Arial", 10),
#                     values = Pvalues_valores, width=5)
#tres2.grid(column=2, row=6, sticky= S+W)
#tres2.current(14)

#vacia2 = Label(group_pval, text = ' ', font = ("Arial", 10))
#vacia2.grid(column=3, row=6, sticky= S+W)


cuatro2 = Label(group_adj, text= 'FDR: ', font=("Arial", 10))
cuatro2.grid(column=1, row=6, sticky= S+W)

ccz = StringVar()#--------------------------------------------------------------
cinco2 = ttk.Combobox(group_adj, textvariable = ccz, font=("Arial", 10),
                     values = valores, width=4)
cinco2.grid(column=2, row=6, sticky= S+W)
cinco2.current(6)

seis2 = Label(group_adj, text=" % ", font=("Arial", 10))
seis2.grid(column=3, row=6, sticky= S+W)

ccgr = IntVar()#--------------------------------------------------------------
siete2 = Checkbutton(group_plots, text=' Plots', font=("Arial", 10),
                    cursor="hand2", variable=ccgr)
siete2.grid(column = 1, row = 6, sticky = W)

######

fila_vacia = Label(root, text=" ")# file vacía
fila_vacia.grid(column=1, row=7)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_adj1 = LabelFrame(root, text="P-value adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_adj1.grid(column=2, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_plots1 = LabelFrame(root, text = "Create") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_plots1.grid(column=3, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_aspect1 = LabelFrame(root, text = "Classification") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect1.grid(column=1, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#group_pval1 = LabelFrame(root, text = "Not adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_pval1.grid(column=2, row=9)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=============================================================================
#oror1 = Label(root,text = " OR ", font = ("Arial", 11, "bold"))
#oror1.grid(column=3, row=9, sticky= W+E+S)

#yyy1 = Label(root, text=" & ", font=("Arial", 12, "bold"))
#yyy1.grid(column=5, row=9, sticky= W+E+S)
#=============================================================================

uno3 = Label(group_aspect1, text = '      KEGG Pathways: ', font=("Arial", 10, "bold"))
uno3.grid(column = 1, row = 9, sticky = E)

#dos3 = Label(group_pval1, text= 'P-value: ', font=("Arial", 10))
#dos3.grid(column = 1, row = 9, sticky= S+W)

#keggpval = StringVar()#--------------------------------------------------------------
#tres3 = ttk.Combobox(group_pval1, textvariable = keggpval, font=("Arial", 10),
#                     values = Pvalues_valores, width=5)
#tres3.grid(column=2, row=9, sticky= S+W)
#tres3.current(14)

#vacia3 = Label(group_pval1, text = ' ', font = ("Arial", 10))
#vacia3.grid(column=3, row=9, sticky= S+W)


cuatro3 = Label(group_adj1, text= 'FDR: ', font=("Arial", 10))
cuatro3.grid(column=1, row=9, sticky= S+W)

keggz = StringVar()#--------------------------------------------------------------
cinco3 = ttk.Combobox(group_adj1, textvariable = keggz, font=("Arial", 10),
                     values = valores, width=4)
cinco3.grid(column=2, row=9, sticky= S+W)
cinco3.current(6)

seis3 = Label(group_adj1, text=" % ", font=("Arial", 10))
seis3.grid(column=3, row=9, sticky= S+W)

kegggr = IntVar()#--------------------------------------------------------------
siete3 = Checkbutton(group_plots1, text=' Plots', font=("Arial", 10),
                    cursor="hand2", variable=kegggr)
siete3.grid(column = 1, row = 9)

###

group_aspect11 = LabelFrame(root, text = "Reference") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect11.grid(column=1, row=10)

group_org = LabelFrame(root, text = "Organism-specific") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_org.grid(column=1, columnspan=3, row=10)


#ocho = Label(group_aspect11, text= '    KEGG Organisms: ', font=("Arial", 10, "bold"))
#ocho.grid(column=1, row=10, sticky= E)    

org_kegg = StringVar()#--------------------------------------------------------------
nueve = ttk.Combobox(group_org, textvariable = org_kegg, font=("Arial", 10),
                     cursor="hand2", values = kegg_organism.Organism.tolist(), height=20, width=50)
nueve.grid(column=1, row=10, sticky= W)
nueve.current(1)


######

fila_vacia1 = Label(root, text=" ")# file vacía
fila_vacia1.grid(column=1, row=11)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_adj2 = LabelFrame(root, text="P-value adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_adj2.grid(column=2, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_plots2 = LabelFrame(root, text = "Create") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_plots2.grid(column=3, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

group_aspect2 = LabelFrame(root, text = "Classification") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect2.grid(column=1, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#group_pval2 = LabelFrame(root, text = "Not adjusted") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#group_pval2.grid(column=2, row=13)            # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=============================================================================
#oror2 = Label(root,text = " OR ", font = ("Arial", 11, "bold"))
#oror2.grid(column=3, row=13, sticky= W+E+S)

#yyy2 = Label(root, text=" & ", font=("Arial", 12, "bold"))
#yyy2.grid(column=5, row=13, sticky= W+E+S)
#=============================================================================

uno4 = Label(group_aspect2, text = '      KEGG Pathways: ', font=("Arial", 10, "bold"))
uno4.grid(column = 1, row = 13, sticky = E)

#dos4 = Label(group_pval2, text= 'P-value: ', font=("Arial", 10))
#dos4.grid(column = 1, row = 13, sticky= S+W)

#keggblpval = StringVar()#--------------------------------------------------------------
#tres4 = ttk.Combobox(group_pval2, textvariable = keggblpval, font=("Arial", 10),
#                     values = Pvalues_valores, width=5)
#tres4.grid(column=2, row=13, sticky= S+W)
#tres4.current(14)

#vacia4 = Label(group_pval2, text = ' ', font = ("Arial", 10))
#vacia4.grid(column=3, row=13, sticky= S+W)


cuatro4 = Label(group_adj2, text= 'FDR: ', font=("Arial", 10))
cuatro4.grid(column=1, row=13, sticky= S+W)

keggblz = StringVar()#--------------------------------------------------------------
cinco4 = ttk.Combobox(group_adj2, textvariable = keggblz, font=("Arial", 10),
                     values = valores, width=4)
cinco4.grid(column=2, row=13, sticky= S+W)
cinco4.current(6)

seis4 = Label(group_adj2, text=" % ", font=("Arial", 10))
seis4.grid(column=3, row=13, sticky= S+W)

keggblgr = IntVar()#--------------------------------------------------------------
siete4 = Checkbutton(group_plots2, text=' Plots', font=("Arial", 10),
                    cursor="hand2", variable=keggblgr)
siete4.grid(column = 1, row = 13)

###

group_aspect22 = LabelFrame(root, text = "Reference") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect22.grid(column=1, row=14)

group_org1 = LabelFrame(root, text = "Organism-specific") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_org1.grid(column=1, columnspan=3, row=14)


#ocho1 = Label(group_aspect22, text= '    KEGG Organisms: ', font=("Arial", 10, "bold"))
#ocho1.grid(column=1, row=14, sticky= E)    

org_keggbl = StringVar()#--------------------------------------------------------------
nueve1 = ttk.Combobox(group_org1, textvariable = org_keggbl, font=("Arial", 10),
                     cursor="hand2", values = kegg_organism.Organism.tolist(), height=20, width=50)
nueve1.grid(column=1, row=14, sticky= W)
nueve1.current(1)

###

group_aspect22 = LabelFrame(root, text = "Tool") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_aspect22.grid(column=1, row=15)

group_method = LabelFrame(root, text = "Method") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_method.grid(column=2, row=15)



diez = Label(group_aspect22, text= '             Local Blast: ', font=("Arial", 10, "bold"))
diez.grid(column=1, row=15, sticky= E)

keggblmethod = IntVar()#--------------------------------------------------------------
tipos = ['Blastx', 'Blastp']
#mets = {0:2,1:3}
mets = {0:1,1:2}
#span = {0:1,1:1}
for i, tips in enumerate(tipos):
    once = Radiobutton(group_method, text=tips, font=("Arial", 10), cursor="hand2",
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
#             text="2. Threshold level of significance", font=("Arial", 15, "bold"))
#label.grid(column=8, row=0, sticky= W+S)

#umbralelegido = IntVar()#--------------------------------------------------------------
#umbrales = ['P-value', 'FDR']
#mets = {0:8,1:9}
#for i, tips in enumerate(umbrales):
#    radio10 = Radiobutton(group_label, text = tips, font=("Arial", 10), cursor="hand2",
#                         variable = umbralelegido, value=i)
#    radio10.grid(column=mets[i], row=1, sticky= W)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fila_vacia1 = Label(root, text=" ")# file vacía
fila_vacia1.grid(column=4, row=2)

############


def inputfile():
    seleccionado = ''.join(re.findall('[A-Z].*[a-z]', opciones[analysis.get()]))
    if seleccionado == 'KEGG Blast Pathways Enrichment':
        file_path = ''
        while file_path == '':
            file_path = filedialog.askopenfilename()
            if file_path.split('.')[-1] == 'fasta':
                new = "filelocation="+file_path+"\n"
                op = open("NeVOmics_params.txt", "r")
                op = op.read()
                if op == '':
                    with open("NeVOmics_params.txt", "a") as f:
                        f.write(new)
                else:
                    newfile = re.sub('filelocation.*\n', new, op)
                    params = open('NeVOmics_params.txt','w')
                    params.write(newfile)
                    params.close()
            else:
                messagebox.showinfo('Status',
                                    'It must be a file in fasta format and with extension .fasta')
                file_path = ''
    else:
        file_path = ''
        while file_path == '':
            file_path = filedialog.askopenfilename()
            if file_path.split('.')[-1] == 'tsv':
                new = "filelocation="+file_path+"\n"
                op = open("NeVOmics_params.txt", "r")
                op = op.read()
                if op == '':
                    with open("NeVOmics_params.txt", "a") as f:
                        f.write(new)
                else:
                    newfile = re.sub('filelocation.*\n', new, op)
                    params = open('NeVOmics_params.txt','w')
                    params.write(newfile)
                    params.close()
            else:
                messagebox.showinfo('Status',
                                    'It must be a file separated by tabs and with extension .tsv')
                file_path = ''

subir_archivo = Label(root, text="2. Select File with UniProt IDs", font=("Arial", 15, "bold"))
subir_archivo.grid(column=5, row=2, sticky= W+N)
subir_archivo1 = Button(root, text="           Upload           ", bg="lightgray", fg="black",
                activebackground= 'yellow',activeforeground= 'black',
                font=("Arial", 10, "bold"), cursor="hand2", command=inputfile)
subir_archivo1.grid(column=5, row=3, sticky= W+N)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fila_vacia1 = Label(root, text=" ")# file vacía
fila_vacia1.grid(column=5, row=4)

############
group_label = LabelFrame(root, text = "Identifier") # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
group_label.grid(column=5, row=6, sticky= W)
######

label = Label(root,
             text="3. Node label in Networks?", font=("Arial", 15, "bold"))
label.grid(column=5, row=5, sticky= W+N)

node_lab = IntVar()#--------------------------------------------------------------
etiquetas = ['Gene Name', 'UniProt ID']
mets = {0:5,1:6}
#span = {0:1,1:1}
for i, tips in enumerate(etiquetas):
    radio7 = Radiobutton(group_label, text = tips, font=("Arial", 10), cursor="hand2",
                         variable = node_lab, value=i)
    radio7.grid(column=mets[i], row=6, sticky= W)

########
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fila_vacia2 = Label(root, text=" ")# file vacía
fila_vacia2.grid(column=4, row=8)

############

def r_exe_loc():
    if os.path.exists('NeVOmics_locRexe.txt'):
        rexe = open('NeVOmics_locRexe.txt', 'r')
        rexe = rexe.read()
        messagebox.showinfo('Status',
                            'The location of R.exe program is already defined in: NeVOmics_locRexe.txt with location \n'+rexe+'\n\n\
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

locat = Label(root,
             text="4. R.exe program location", font=("Arial", 15, "bold"))
locat.grid(column=5, row=8, sticky= W+N)

boton2 = Button(root, text=" Click to search R.exe ", bg="lightgray", fg="black", font=("Arial", 10, "bold"),
                activebackground= 'yellow',activeforeground= 'black',
                cursor="hand2", command = r_exe_loc)
boton2.grid(column=5, row=9, sticky= W+N)


######
fila_vacia3 = Label(root, text=" ")# file vacía
fila_vacia3.grid(column=8, row=11)

########################################################################################



def parameters():
    oo = open('NeVOmics_params.txt', 'r')
    if re.search('filelocation.*\n', oo.read()):
        oo = open('NeVOmics_params.txt', 'r') # abre el archivo y agrega los parametros de corrida
        savefile = re.search('filelocation.*\n', oo.read()).group()
        params = open('NeVOmics_params.txt','w')
        params.close()
        with open("NeVOmics_params.txt", "a") as f:
            f.write(savefile)
            f.write("#====================\n")
            f.write("analysis="+' '.join(re.findall('\w+', opciones[analysis.get()]))+"\n")
            f.write("#num="+str(analysis.get())+"\n")
            f.write("#====================\n")
            f.write("#**GO ENRICHMENT**\n")
            #f.write("govalforumbral="+val_for_umbral[umbral.get()]+"\n")
            #f.write("bppval="+str(bppval.get())+"\n")
            #f.write("mfpval="+str(mfpval.get())+"\n")
            #f.write("ccpval="+str(ccpval.get())+"\n")
            f.write("bpfdr="+str(bpz.get())+"\n")
            f.write("mffdr="+str(mfz.get())+"\n")
            f.write("ccfdr="+str(ccz.get())+"\n")
            f.write("bpplots="+str(bpgr.get())+"\n")
            f.write("mfplots="+str(mfgr.get())+"\n")
            f.write("ccplots="+str(ccgr.get())+"\n")
            f.write("#====================\n")
            f.write("#**KEGG ENRICHMENT**\n")
            #f.write("keggvalforumbral="+val_for_umbral[umbral.get()]+"\n")
            #f.write("keggpvalue="+str(keggpval.get())+"\n")
            f.write("keggfdr="+str(keggz.get())+"\n")
            f.write("keggplots="+str(kegggr.get())+"\n")
            f.write("keggorganism="+org_kegg.get()+"\n")
            f.write("keggprefix="+dict_org[org_kegg.get()][0]+"\n")
            f.write("keggTnumber="+dict_org[org_kegg.get()][1]+"\n")
            f.write("#====================\n")
            f.write("#**KEGG BLAST ENRICHMENT**\n")
            #f.write("keggblvalforumbral="+val_for_umbral[umbral.get()]+"\n")
            #f.write("keggblastpvalue="+str(keggblpval.get())+"\n")
            f.write("keggblastfdr="+str(keggblz.get())+"\n")
            f.write("keggblastplots="+str(keggblgr.get())+"\n")
            f.write("keggmethodblast="+tipos[keggblmethod.get()]+"\n")
            f.write("keggblastorganism="+org_keggbl.get()+"\n")
            f.write("keggblastprefix="+dict_org[org_keggbl.get()][0]+"\n")
            f.write("keggblastTnumber="+dict_org[org_keggbl.get()][1]+"\n")
            f.write("#====================\n")
            #f.write("umbralelegido="+umbrales[umbralelegido.get()]+"\n")
            f.write("labelnode="+etiquetas[node_lab.get()]+"\n")
            f.close
        
        
        #comando = 'python short_KEGG.py'
        #os.system("start cmd /c "+comando)
        
        run = subprocess.call(['python', 'short_KEGG.py'])
    else:
        messagebox.showinfo('Status',
                            'No file were selected.')
        
        # aquí incluir los scripts a ejecutar, dependiendo del análisis seleccionado


boton3 = Button(root, text="            RUN            ", cursor="hand2",
                activebackground= 'red',activeforeground= 'black',
                bg="green", fg="white",font=("Arial", 20, "bold"), command = parameters) # , command = parameters
boton3.grid(column=5, row=11, rowspan=2, sticky= N+S+E+W)    


def callback(event):
    webbrowser.open_new(r"https://doi.org/10.3390/genes9120569")
link = Label(root, text="If you use NeVOmics, please cite.\n",
             font=("Arial", 10), fg="blue", cursor="hand2")
link.grid(column=5, row=14, sticky= E+S)
link.bind("<Button-1>", callback)


root.mainloop()


# In[ ]:





# In[13]:


kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()
kegg_organism=pd.read_csv(StringIO(kegg_orgs),names=['T_number','Prefix','Organism','Group'],sep='\t')
kegg_organism = kegg_organism.sort_values(by ='Organism',ascending=True)


# In[16]:


kegg_organism


# In[ ]:




