#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append("../NeVOmics_PyMod/")

print("\nParameters\n")
import re
from pandas import Series, DataFrame 
import pandas as pd
import csv
import pathlib
import urllib.request
import webbrowser
import shutil, os
import numpy as np
from urllib.request import urlopen
import requests
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import warnings
from datetime import datetime 
inicio_total = datetime.now()
import os, fnmatch
import tkinter as tk
from tkinter import *
from tkinter import messagebox
import warnings
warnings.filterwarnings("ignore")
from networkx import path_graph, random_layout
import networkx as nx
from matplotlib import cm
import matplotlib
from colormap import Colormap
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib as mpl


# In[2]:


def del_stop_process():
    if os.path.exists("short_KEGG.py"): os.remove("short_KEGG.py")
    if os.path.exists("HD.py"): os.remove("HD.py")
    sys.exit()


# In[127]:


parameters = open('NeVOmics_params.txt', 'r')
parametros = parameters.read()
parameters.close()
#print(parametros)


han = open('../NeVOmics_img/KEGG_Organisms.txt', 'r')
dict_org = {}
for line in han:
    line = line.rstrip()
    if re.search('^#', line):
        infokegg = line
        pass
    else:
        separados = line.split('\t')
        dict_org[separados[2]] = [separados[1], separados[0]]
han.close()


# ## definimos los parametros para  KEGG

# In[129]:


#######################
# parametros elegidos
method_P = 'FDR' # default

# definimos la localizacion y nombre del archivo elegido
file_path = re.search('filelocation.*', parametros).group().split('=')[1]

# colores elegidos por el usuario para los terminos y edges
usercolormap = re.search('edgecolor.*', parametros).group().split('=')[1]

# coloes para rampa, asignacion a valores, numericos relacionados a los genes o proteinas
colormap_definido = re.search('networkcolor.*', parametros).group().split('=')[1]

# titulo para la barra de colormap
barcolortitle = re.search('usertext.*', parametros).group().split('=')[1]

# definicion de nodos si no hay valores asociados a los nodos de la red, default blue
nodecolorsinback = re.search('uniquecolor.*', parametros).group().split('=')[1]

# fdr elegido
FDR = float(re.search('keggfdr.*', parametros).group().split('=')[1]) / 100

# organismo seleccionado
organism = re.search('keggorganism.*', parametros).group().split('=')[1]

# frefijo identificado a partir del organismo
#pref = re.search('keggprefix.*', parametros).group().split('=')[1]
pref = dict_org[organism][0]

# T number identificado a partir del organismo
#t_number = re.search('keggTnumber.*', parametros).group().split('=')[1]
t_number = dict_org[organism][1]

# crear los gráficos?
# si la respuesta es 0 no se crearán, si es 1 se crearán
keggplots = re.search('keggplots.*', parametros).group().split('=')[1]

# crear redes? 0 es no, 1 es si
createnetworks = re.search('networksplots.*', parametros).group().split('=')[1]

# crear circos? 0 es no, 1 es si
createcircos = re.search('circosplots.*', parametros).group().split('=')[1]

# etiqueta de nodos ['Gene Name', 'UniProt ID']
labelnode = re.search('labelnode.*', parametros).group().split('=')[1]

print('method_P =', method_P)
print('file_path =', file_path)
print('usercolormap =', usercolormap)
print('colormap_definido =', colormap_definido)
print('barcolortitle =', barcolortitle)
print('nodecolorsinback =', nodecolorsinback)
print('FDR =', FDR)
print('organism =', organism)
print('pref =', pref)
print('t_number =', t_number)
print('createplots =', keggplots)
print('createnetworks =', createnetworks)
print('createcircos =', createcircos)
print('labelnode =', labelnode)


# In[30]:

print("\n\nWait ...\n")
## read file
inp_file=pd.read_csv(file_path,sep='\t',header=None)   
    
## explore input file
if len(inp_file.columns) == 1:
    hayvalores = 'nohayvalores'
    inp_file['values'] = 1
    ## only gene list
    list_input=inp_file.rename(columns={0:'Entry'},index=str) 
if len(inp_file.columns) == 2:
    hayvalores = 'sihay'
    ## gene list and vales
    list_input=inp_file.rename(columns={0:'Entry',1:'values'},index=str) ## gene list and values
if len(inp_file.columns) == 3:
    hayvalores = 'sihay'
    ## gene list, values and background
    list_input=inp_file.rename(columns={0:'Entry',1:'values',2:'Background'},index=str) 


# In[31]:


"""
extrae el numero taxonomico y el nombre del organismo desde uniprot
para compararlo con el ingresado por el usuario,
si no son iguales el proceso se detiene
"""
id_organism = requests.get("https://www.uniprot.org/uniprot/?query="+list_input.Entry[0]+"&sort=score&columns=organism-id,organism&format=tab&limit=1").content.decode()
Prefix = id_organism.split('\t', 2)[1].split('\n')[1]


k = requests.get("https://www.kegg.jp/dbget-bin/www_bget?gn:"+t_number)
k = k.text.rstrip()
tax = re.findall('TAX:.*', k)[0]
tax2 = re.findall('Info&id=\d+..\d+', tax)[0]
Prefix_user = tax2.split('">')[1]

root = Tk()
root.withdraw()
if Prefix_user != Prefix:
    messagebox.showwarning('Status',
                        'Your selected organism ('+organism+') does not correspond to the \
 organism identified with the UniProt identifiers ('+id_organism.rstrip().split('\t')[-1]+').\n\n\
 !!!Choose the organism correctly!!!')
    del_stop_process()
else:
    pass


## Create a folder
os.makedirs('data',exist_ok=True)




# In[33]:


# descarga el modulo para la estadistica
hd = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/HD.py', './HD.py')


# In[34]:


# all kegg-id and pathway-description
ee=requests.get('http://rest.kegg.jp/list/pathway/'+pref+'').content.decode()
ee = re.sub('path:|- '+organism[0:5]+'.*','',ee)
kegg_pathways = DataFrame([i.split('\t') for i in ee.split('\n')], columns = ['Path','Term'])
kegg_pathways.to_csv('data/Pathways.txt',sep='\t',index=None)


# In[35]:


# info version
print(re.sub('#', '', infokegg))

# In[36]:


######## extraccion de informacion
if Prefix == '9606': # human
    inf1=requests.get('https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:'+Prefix+'&format=tab&columns=id,database(DisGeNET)').content.decode().rstrip()
    guardar = []
    for i in inf1.split('\n')[1:]:
        if i.split('\t')[1] == '':
            uno = np.nan
        else:
            uno = i.split('\t')[1].split(';')[0]
        guardar.append([i.split('\t')[0], uno])
    inf1 = DataFrame(guardar, columns = ['Entry', 'Entry_Kegg'])
    inf2=requests.get('https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:'+Prefix+'&format=tab&columns=id,database(GeneID)').content.decode().rstrip()
    guardar2 = []
    for i in inf2.split('\n')[1:]:
        if i.split('\t')[1] == '':
            uno = np.nan
        else:
            uno = i.split('\t')[1].split(';')[0]
        guardar2.append([i.split('\t')[0], uno])
    inf2 = DataFrame(guardar2, columns = ['Entry', 'Entry_Kegg'])
    Kegg_Uniprot=pd.concat([inf1,inf2], axis=0,sort=False).dropna().reset_index(drop=True).drop_duplicates()
    dd=requests.get('http://rest.kegg.jp/link/pathway/'+pref+'').content.decode().rstrip()
    guardar3 = [(re.sub('^'+pref+':', '', i.split('\t')[0]), re.sub('path:', '', i.split('\t')[1])) for i in dd.split('\n')]
    kegg_path_ID = DataFrame(guardar3, columns = ['Entry_Kegg','Path'])
else:
    # all kegg-id and pathway-id
    dd=requests.get('http://rest.kegg.jp/link/pathway/'+pref+'').content.decode().rstrip()
    guardar3 = [(re.sub('^'+pref+':', '', i.split('\t')[0]), re.sub('path:', '', i.split('\t')[1])) for i in dd.split('\n')]
    kegg_path_ID = DataFrame(guardar3, columns = ['Entry_Kegg','Path'])
    # all kegg-id and pathway-description
    # all kegg-id and uniprot
    ff=requests.get('http://rest.kegg.jp/conv/uniprot/'+pref+'').content.decode().rstrip()
    guardar4 = [(re.sub('^'+pref+':', '', i.split('\t')[0]), re.sub('up:', '', i.split('\t')[1])) for i in ff.split('\n')]
    Kegg_Uniprot = DataFrame(guardar4, columns = ['Entry_Kegg','Entry'])

allanotacion = Kegg_Uniprot.merge(kegg_path_ID, on = 'Entry_Kegg', how = 'left').dropna()


# In[37]:


"""
corroborar que el organismo ingresado sea igual al detectado en el archivo,
si no es igual termina el proceso
"""
id_organism_user = requests.get("https://www.uniprot.org/uniprot/?query="+Kegg_Uniprot.Entry[0]+"&sort=score&columns=organism-id,organism&format=tab&limit=1").content.decode()
Prefix_user = id_organism_user.split('\t', 2)[1].split('\n')[1]


    
## definimos el tipo de etiqueta para los nodos
labnode = re.search('labelnode.*', parametros).group().split('=')[1]


# In[39]:


###############################################


# In[ ]:





# In[40]:


if len(inp_file.columns) == 3:
    print('Data with background column\n')
    provicional = list_input[['Background']].rename(columns={'Background':'Entry'})
    background_info = provicional.merge(allanotacion, on = 'Entry', how = 'inner')
    
    # guardar archivo background, con la columna de genes
    background_info[['Entry_Kegg']].drop_duplicates().to_csv('data/Background.txt',index=None)

    # 2.- Preparation of list with pathways
    ## Protein list mapping against "background_info" and then save this list
    list_input_match = list_input[['Entry']].merge(background_info,how="left", on='Entry').dropna().drop_duplicates()
    list_input_match[['Entry_Kegg']].drop_duplicates().to_csv('data/List.txt',index=None)
    
  
    # 3.- background with: Entry	GO, for association file
    background_info[['Entry_Kegg','Path']].to_csv('data/Association.txt',index=None,sep='\t')
else:
    print('Data without background column\n')
    background_info = allanotacion
    
    # guardar archivo background, con la columna de genes
    background_info[['Entry_Kegg']].drop_duplicates().to_csv('data/Background.txt',index=None)

    # 2.- Preparation of list with pathways
    ## Protein list mapping against "background_info" and then save this list
    list_input_match = list_input[['Entry']].merge(background_info,how="left", on='Entry').dropna().drop_duplicates()
    list_input_match[['Entry_Kegg']].drop_duplicates().to_csv('data/List.txt',index=None)
    
  
    # 3.- background with: Entry	GO, for association file
    background_info[['Entry_Kegg','Path']].to_csv('data/Association.txt',index=None,sep='\t')


# In[ ]:





# In[60]:


#exploratory analysis of P-value in data
analysis = 'Pathways.txt'
subprocess.call(["python", "HD.py", analysis,str(FDR)])


# In[61]:


# abrimos los resultados
enrich_P = pd.read_csv('data/Enrichment_analysis_'+analysis.split('.')[0]+'.tsv',sep='\t')


# #  <font color = red>En este paso si no hay mas de un termino terminar el proceso porque no se pueden realizar redes con un nodo/proteinas<font>

# In[355]:


if enrich_P[enrich_P.Sig == 'T']['FDR'].count() >= 1: # al menos un valor de FDR es significativo      
    results_process_P = enrich_P[enrich_P.Sig == 'T']
    vertice = [(k, j) for i, k in zip(results_process_P.entry.tolist(), results_process_P.base.tolist()) for j in i.split(';')]
    proteins_count_P = len(set([i[1] for i in vertice]))
    GO_count_P = results_process_P.base.count()
else:
    # sin informacion en el kegg
    results_process_P = enrich_P
    no_anotadas = []
    for i in list_input.Entry.drop_duplicates().dropna().tolist():
        if i in list_input_match.Entry.drop_duplicates().tolist():
            continue
        else:
            no_anotadas.append(i)
    
    if len(results_process_P) == 0:
        singleton = 0
    else:
        singleton = int(float(results_process_P.Bonf_corr.iloc[-1:]) / float(results_process_P.P.iloc[-1:]))
    
    
    reporte = {'base':[np.nan,
                       'KEGG DB Last-Modified',
                       'Input file name',
                       'Association file name',
                       'Total number of background',
                       'Total number of list',
                       'Background with Pathways',
                       'List input with Pathways',
                       'Non-singletons value for Bonf_corr',
                       'Correction Method',
                       'Value',
                       np.nan,
                       'Proteins with no information in KEGG Pathways',
                       ';'.join(no_anotadas)],
               'list_count':[np.nan,
                             infokegg,
                             file_path, analysis,
                             background_info['Entry'].drop_duplicates().count(),
                             list_input['Entry'].drop_duplicates().count(),
                             background_info['Entry'].drop_duplicates().count(),
                             list_input_match['Entry'].drop_duplicates().count(),
                             int(float(results_process_P.Bonf_corr.iloc[-1:]) / float(results_process_P.P.iloc[-1:])),
                             'FDR',
                             str(FDR)+' ('+str(FDR * 100)+'%)',
                             np.nan, 
                             len(no_anotadas),
                             np.nan]}
    information = DataFrame(reporte)
    informe_final = pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'base':'Path'})
    informe_final = informe_final[['Path', 'list_count', 'back_count', 'tot_list', 'tot_back', 'P', 'Bonf_corr',
           'Rank', 'FDR', 'Sig', 'Term', 'entry']]
    writer = pd.ExcelWriter('Enrichment_Pathways_Analysis_FDR_'+str(FDR)+'.xlsx')

    informe_final.to_excel(writer,'Significant KEGG Pathways',index=False)

    enrich_P.to_excel(writer,'Enrichment Results',index=False)
    
    writer.save()
    del_stop_process()
    



etiquetas = []
for i in results_process_P.Term:
    i = i.rstrip()
    if len(i.split(' ')) == 1:
        etiquetas.append(i)
    if len(i.split(' ')) == 2:
        etiquetas.append(re.sub(' ', '\n', i))
    if len(i.split(' ')) == 3:
        etiquetas.append(re.sub(' ', '\n', i))
    if len(i.split(' ')) == 4:
        etiquetas.append(' '.join(i.split(' ')[0:2])+'\n'+' '.join(i.split(' ')[2:4]))
    if len(i.split(' ')) > 4:
        etiquetas.append(' '.join(i.split()[0:2])+'\n'+' '.join(i.split()[2:4])+'...')
results_process_P['Short_Term'] = etiquetas



# In[150]:


significativos = []
for x in results_process_P.base.drop_duplicates():
    dff = results_process_P[results_process_P.base == x]
    for index, row in dff.iterrows():
        for i in row.entry.split(';'):
            significativos.append([x, row.P, row.FDR, row.Term, row.Short_Term, i])


keggtabla = DataFrame(significativos, columns = ['Path', 'P', 'FDR', 'Term', 'Short_Term', 'Entry_Kegg'])
keggtabla['LogminFDR'] = -np.log10(keggtabla.FDR)
keggtabla['LogminP'] = -np.log10(keggtabla.P)
n = 0
ranked = []
for i in keggtabla['Entry_Kegg'].drop_duplicates():
    n+=1
    ranked.append([i, str(n)])
rank = DataFrame(ranked, columns = ['Entry_Kegg', 'label'])


keggtabla = keggtabla.merge(rank, on = 'Entry_Kegg', how = 'left')
keggtabla = keggtabla.merge(list_input_match, on = ['Entry_Kegg', 'Path'], how = 'left')

keggtabla = keggtabla.merge(list_input[['Entry', 'values']], on = 'Entry', how = 'left')

edges_frame_excel = keggtabla[['Path','Entry_Kegg','Entry','Term','values']]


# In[151]:


if labelnode == 'Gene Name':
    pass
if labelnode == 'UniProt ID':
    keggtabla = keggtabla.rename({'Entry_Kegg':'Entry', 'Entry':'Entry_Kegg'}, axis='columns')


# In[66]:


###########################################################################################################


# In[67]:


#####################################################################################


# In[68]:


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

# estas rampas tienen una cantidad de colores predefinidos
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
Spectral = [matplotlib.colors.to_hex(i) for i in cm.Spectral(np.arange(12)/12.)]

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

sequentials_colors.update(diverging_colors)
sequentials_colors.update(uniform_sequential)
sequentials_colors.update(qualitative_colors)
sequentials_colors.update(edge_colors)


# In[ ]:



matrix = keggtabla.pivot_table(values='Entry',index=['label'],aggfunc=len,columns=['Path', 'Term', 'Short_Term'])

df_mat = []
for i in list(matrix.columns.values):
    new = DataFrame(matrix[i])
    for x, y in zip(list(new[i].index), list(new[i].values)):
        df_mat.append([x, i, y])

df_mat = DataFrame(df_mat, columns = ['go0', 'go1', 'val']).dropna()

nodos = []
for index, row in df_mat.iterrows():
    if row.go0 == row.go1:
        #print(row.go0, row.go1)
        continue
    else:
        #print(row.go0, row.go1)
        nodos.append([row.go0, row.go1, row.val])
nodos = DataFrame(nodos)
nodos = DataFrame([[i for i in nodos[0]],
                   [i[0] for i in nodos[1]],
                    nodos[2]]).T

G=nx.Graph()
for index, row in nodos.iterrows():
    G.add_edge(str(row[0]), row[1],weight = row[2])
#esmall=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] > 0]


xxx = []
for i in G.nodes():
    xxx.append(str(i))
yyy = DataFrame(xxx, columns = ['label'])

# si hay columna de valores numéricos en el input
if 'values' in list(list_input.select_dtypes('number').columns):
    zzz = yyy.merge(keggtabla[['label', 'values']], on = 'label', how = 'left').drop_duplicates().reset_index(drop = True)
    zzz = zzz.sort_values(by ='values',ascending=False).reset_index(drop=True)
else:
    zzz = yyy.merge(keggtabla[['label']], on = 'label', how = 'left').drop_duplicates().reset_index(drop = True)

G=nx.Graph()
for index, row in nodos.iterrows():
    G.add_edge(str(row[0]), row[1],weight = row[2])
#esmall=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] > 0]


xxx = []
for i in G.nodes():
    xxx.append(str(i))
yyy = DataFrame(xxx, columns = ['label'])
zzz = yyy.merge(keggtabla[['label', 'values']], on = 'label', how = 'left').drop_duplicates().reset_index(drop = True)
zzz = zzz.sort_values(by ='values',ascending=False).reset_index(drop=True)

###################

mycmap = sequentials_colors[colormap_definido].reversed()

nulos = len(zzz['values']) - len(zzz['values'].dropna())
null_col = list(np.repeat('', nulos))
ids = list(np.round(np.linspace(zzz['values'].max(), zzz['values'].min(), len(zzz['values'])*4), 50))

rangoforcolor = []
valor_unico = []
n = ids[0]
for i in ids:
    if i == n:
        valor_unico.append(i)
        continue
    rangoforcolor.append([n, i])
    n = i

#*******************************************************
mycmap0 = mycmap(np.linspace(0, 1, len(rangoforcolor)))
##*********************************************************

colores = []
for i in mycmap0:
    colores.append(matplotlib.colors.to_hex(i))
rangos = {}
for k, j in zip(rangoforcolor, colores):
    if len(k) == 2:
        rangos[str(k[0])+','+str(k[1])] = j
    if len(k) == 1:
        rangos[str(k[0])] = j   

        
if len(set(zzz.dropna()['values'].tolist())) == 1:
    zzz['cols'] = list(np.repeat(nodecolorsinback, len(zzz['values'].dropna()))) + null_col
else:
    positivos = []
    for i in rangoforcolor:
        if len(i) == 2:
            for j in [(np.round(x, 50), X) for x, X in zip(zzz['values'], zzz.label) if x > 0]:
                if i[0] >= j[0] >= i[1]:
                    #print(rangos[str(i[0])+','+str(i[1])],  j)
                    positivos.append([j[1], rangos[str(i[0])+','+str(i[1])]])
                    #positivos.append(rangos[str(i[0])+','+str(i[1])])
    positivos = DataFrame(positivos).drop_duplicates()[1].tolist()
    negativos = []
    for i in rangoforcolor:
        if len(i) == 2:
            for j in [np.round(x, 50) for x in zzz['values'] if x < 0]:
                if i[0] >= j >= i[1]:
                    #print(rangos[str(i[0])+','+str(i[1])],  j)
                    negativos.append(rangos[str(i[0])+','+str(i[1])])
    zzz['cols'] = positivos + negativos + null_col
    
######
n = 0
l = 10
k = 15
for i in range(10):
    if n+1 <= len(G.edges()) <= l:
        #print(k)
        valor = k
    #print(n+1, k, l)
    n += 10
    l += 10
    k -= 1
if 101 <= len(G.edges()) <= 200:
    valor = 5
if 201 <= len(G.edges()) <= 300:
    valor = 4
if len(G.edges()) >= 301:
    valor = 3
####
    
if hayvalores == 'nohayvalores':
    colorletra = nodecolorsinback
else:
    colorletra = 'none'   


# In[ ]:





# In[70]:


##############################################################


# In[71]:


# asignacion de colores a cada entry, menos a path
colorder = dict(zip(zzz.label.tolist(), zzz.cols.tolist()))


# In[72]:


# posiciones en el plano cartesiano

#pos = nx.spring_layout(G,iterations=50) #nx.spring_layout(G, k = 0.3, iterations=50, threshold=0.001, seed= 123)
pos = nx.kamada_kawai_layout(G, dist=None, weight='weight', scale=1, center=None, dim=2)
#pos = nx.circular_layout(G)


# In[73]:


# ordenados por FDR
paths = results_process_P.base.drop_duplicates().tolist()


# In[74]:


labnodeterms = dict(zip(paths, edge_colors[usercolormap][0:len(paths)]))


# In[ ]:





# In[77]:


name_term = dict(zip(paths, results_process_P.Term.drop_duplicates().tolist()))
# In[78]:


# size de nodo, amplificado 200 veces
sizenodo = -np.log10(np.array(results_process_P.FDR))


######################################

if labelnode == 'Gene Name':
    genelist = []
    for i in paths:
        df = keggtabla[keggtabla.Path == i]
        geneslist = (df.label.tolist(), df.Entry_Kegg.tolist())
        genelist.append([i, geneslist])
    info_for_url = []
    links_for_entrys = {}
    for i in genelist:
        for j, k in zip(i[1][0], i[1][1]):
            links_for_entrys[j] = [k, 'https://www.genome.jp/dbget-bin/www_bget?'+pref+':'+k]
            info_for_url.append([i[0], k+colorder[j]])
    df_update = DataFrame(info_for_url, columns = ['Path', 'gene_exp'])
    # una lista de los genes para mostrar en los mapas de las vías, , aquí los genes son Entry_Kegg
    genelist_map = []
    for i in paths:
        df = keggtabla[keggtabla.Path == i]
        geneslist = (df.label.tolist(), df.Entry_Kegg.tolist())
        genelist_map.append([i, geneslist])
    info_for_url_map = []
    for i in genelist_map:
        for j, k in zip(i[1][0], i[1][1]):
            info_for_url_map.append([i[0], k+colorder[j]])
    df_update_map = DataFrame(info_for_url_map, columns = ['Path', 'gene_exp'])
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
if labelnode == 'UniProt ID':
    genelist = []
    for i in paths:
        df = keggtabla[keggtabla.Path == i]
        geneslist = (df.label.tolist(), df.Entry_Kegg.tolist())
        genelist.append([i, geneslist])
    info_for_url = []
    links_for_entrys = {}
    for i in genelist:
        for j, k in zip(i[1][0], i[1][1]):
            links_for_entrys[j] = [k, 'https://www.uniprot.org/uniprot/'+k]
            info_for_url.append([i[0], k+colorder[j]])
    df_update = DataFrame(info_for_url, columns = ['Path', 'gene_exp'])
    # una lista de los genes para mostrar en los mapas de las vías, aquí los genes son Entry
    genelist_map = []
    for i in paths:
        df = keggtabla[keggtabla.Path == i]
        geneslist = (df.label.tolist(), df.Entry.tolist())
        genelist_map.append([i, geneslist])
    info_for_url_map = []
    for i in genelist_map:
        for j, k in zip(i[1][0], i[1][1]):
            info_for_url_map.append([i[0], k+colorder[j]])
    df_update_map = DataFrame(info_for_url_map, columns = ['Path', 'gene_exp'])


# In[80]:


url_for_kegg = {}
for i in paths:
    fijo = 'https://www.kegg.jp/kegg-bin/show_pathway?map='+i+'&multi_query='
    #print(fijo)
    df = df_update_map[df_update_map.Path == i]
    uno = ',gainsboro%0D'.join(df.gene_exp.tolist())
    dos = re.sub('#', '+%23', uno)
    tres = re.sub('$', ',gainsboro', dos)
    url_for_kegg[i] = fijo+tres


# In[81]:


valor_maximo = zzz['values'].max().round()
valor_minimo = zzz['values'].min().round()


# In[82]:


verti = []
for i, j in G.edges():
    if i in paths:
        verti.append(labnodeterms[i])
    if j in paths:
        verti.append(labnodeterms[j])


# In[83]:


ccc = []
for i in G.nodes():
    if i in list(labnodeterms.keys()):
        ccc.append(labnodeterms[i])
    else:
        continue


# In[ ]:





# In[ ]:





# # Preparacion y creacion del archivo excel

# ### Datos para el reporte en la hoja de excel

# In[121]:


# sin informacion en el kegg
no_anotadas = []
for i in list_input.Entry.drop_duplicates().dropna().tolist():
    if i in list_input_match.Entry.drop_duplicates().tolist():
        continue
    else:
        no_anotadas.append(i)


# In[122]:

reporte = {'base':[np.nan,
                   'KEGG DB Last-Modified',
                   'Input file name',
                   'Association file name',
                   'Total number of background',
                   'Total number of list',
                   'Background with Pathways',
                   'List input with Pathways',
                   'Non-singletons value for Bonf_corr',
                   'Correction Method',
                   'Value',
                   np.nan,
                   'Proteins with no information in KEGG Pathways',
                   ';'.join(no_anotadas)],
        'list_count':[np.nan,
                      infokegg,
                      file_path, analysis,
                      background_info['Entry'].drop_duplicates().count(),
                      list_input['Entry'].drop_duplicates().count(),
                      background_info['Entry'].drop_duplicates().count(),
                      list_input_match['Entry'].drop_duplicates().count(),
                      int(float(results_process_P.Bonf_corr.iloc[-1:]) / float(results_process_P.P.iloc[-1:])),
                      'FDR',
                      str(FDR)+' ('+str(FDR * 100)+'%)',
                      np.nan,
                      len(no_anotadas),
                      np.nan]}
information = DataFrame(reporte)
informe_final = pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'base':'Path'})
# In[ ]:





# In[139]:


writer = pd.ExcelWriter('Enrichment_Pathways_Analysis_FDR_'+str(FDR)+'.xlsx')

informe_final.to_excel(writer,'Significant KEGG Pathways',index=False)

enrich_P.to_excel(writer,'Enrichment Results',index=False)

edges_frame_excel.to_excel(writer,'Edges Pathways',index=False)
writer.save()


# ### Edges y noded para R

# In[126]:


from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
# Preparacion de colores para bar colormap en R, la informacion la obtengo del 
#diccionario creado preciamente llamado <font color = red>"sequentials_colors"<font>
# con este comando extraigo una lista de colores para la barra colormap en R, la defino desde python

color_for_bar_in_R = [matplotlib.colors.to_hex(i) for i in sequentials_colors[colormap_definido](np.linspace(0, 1, 50))]

colores_bar_R = DataFrame(color_for_bar_in_R, columns = ['bar_color_R'])


# In[ ]:





# In[159]:


if createcircos == '1':
    # agrego los colores de los genes/proteinas para que en R sea leida y los colores ya esten definidos
    colores_entry = []
    for i in keggtabla.label:
        colores_entry.append(colorder[i])
    keggtabla['entry_colors'] = colores_entry

    # agrego el titulo elegido por el usuario
    anex_bartitle = DataFrame([barcolortitle], columns = ['bar_title'])

    # preparacion del archivo de nodos
    edges_file_name = 'for_R_edges_KEGG_Enrich_Analysis_'+''.join(method_P)+'_'+str(FDR)+'.csv'
    edges_frame = keggtabla[['Path','Entry_Kegg','Entry','Term', 'values', 'entry_colors']].drop_duplicates()

    if hayvalores == 'sihay':
        # agregar en el archivo edges la rampa de colores para R
        edges_frame = pd.concat([edges_frame,colores_bar_R], axis=1)
        edges_frame = pd.concat([edges_frame,anex_bartitle], axis=1)
    if hayvalores == 'nohayvalores':
        # agregar en el archivo edges una columna de colores unicos definicos como "nodecolorsinback" 
        edges_frame['bar_color_R'] = nodecolorsinback
        edges_frame = pd.concat([edges_frame,anex_bartitle], axis=1)
    
    edges_frame.to_csv(edges_file_name,index=None)

    # agrego los colores de los terminos a la tabla para que en R sea leida y los colores ya esten definidos
    results_process_P['term_colors'] = list(labnodeterms.values())

    nodes_file_name='for_R_nodes_KEGG_Enrich_Analysis_'+''.join(method_P)+'_'+str(FDR)+'.csv'
    results_process_P.drop_duplicates().to_csv(nodes_file_name,index=None)

if createcircos == '0': # el usuario decició no crear estos gráficos
    pass


# In[160]:


if createnetworks == '1':
    print('\n------------')
    print('Python Plots')
    ## Create a folder
    new_dir_plots = "job_KEGG_plots"
    os.makedirs(new_dir_plots,exist_ok=True)
    
    
    # # Figuras 1 y 2
    
    NEWSIZE = [valor, 0]
    NEWPAD = [0.1, valor*0.7]
    NEWCOLOR = ['white', colorletra]
    IMGLABEL = [1, 2]
    
    acumulacion = ''
    for newsize, newpad, newcolor, imglabel in zip(NEWSIZE, NEWPAD, NEWCOLOR, IMGLABEL):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, .5, 1])
        
        nx.draw_networkx_nodes(G, pos, nodelist = paths,
                               node_size =  list(np.array(sizenodo) * 200),
                               node_shape = 'o',
                               node_color = edge_colors[usercolormap][0:len(paths)],
                               linewidths = 0.5, edgecolors= 'white')
    
        nx.draw_networkx_edges(G, pos, width = 2, alpha= 0.5, edge_color= verti)
        ax.axis('equal')
        
        for nod_entry in links_for_entrys:
            ax.annotate(nod_entry,
                        xy=pos[nod_entry], xycoords='data',
                        url = links_for_entrys[nod_entry][1],
                        xytext=pos[nod_entry], textcoords='data',
                        color = newcolor, #..................................
                        fontweight='bold',
                        size = newsize, #..................................
                        ha='center', va="center",
                        bbox=dict(boxstyle="circle",
                                  pad = newpad, #..................................
                                  alpha=1, fc=colorder[nod_entry],
                        ec="none", url = links_for_entrys[nod_entry][1]))
            
        ax.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax2 = fig.add_axes([0.5, 0.77, 0.015, .18]) # vertical
        
        if len(valor_unico) > 1:
            ax2.axis('off')
            posiciones_leyenda_de_abajo = [.5, 0.2, 0.15, 0.73]
            pass
    
        else:
            norm = mpl.colors.Normalize(vmax = valor_maximo,
                                        vmin = valor_minimo)
            cb1 = mpl.colorbar.ColorbarBase(ax2,
                                            cmap=ListedColormap(mycmap(np.linspace(1, 0, len(zzz['values'])*2))),
                                            norm=norm, spacing='proportional',
                                            orientation='vertical')
            plt.tick_params(axis="y", color="grey")
            cb1.set_label(barcolortitle, labelpad=-30, y=1.2, rotation=0, size=10, fontweight='bold')
            plt.yticks(size = 7)
            cb1.outline.set_linewidth(0)
            
            posiciones_leyenda_de_abajo = [0.5, 0, 0.15, 0.73]
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax3 = fig.add_axes(posiciones_leyenda_de_abajo)
        
        ax3.text(0,0.95, 'KEGG Pathways', size=12,ha='left',color= 'black', fontweight='bold')
        
        n = 0.05
        for t in paths[0:20]:
            n += 0.043
            ax3.annotate(t,
                        xy=np.array([0,1-n]), xycoords='data',
                        url = url_for_kegg[t],
                        xytext=np.array([0,1-n]), textcoords='data',
                        size=8, ha='left', va="center",
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[t], ec="none",
                                 url = url_for_kegg[t]))
        
        ax3.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
        
        
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
        plt.close()
    
    
    # In[ ]:
    
    
    
    
    
    
    # # Figuras 3 y 4
    
    # In[222]:
    
    
    IMGLABEL = [3, 4]
    for newsize, newpad, newcolor, imglabel in zip(NEWSIZE, NEWPAD, NEWCOLOR, IMGLABEL):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, .5, 1])
        
        nx.draw_networkx_nodes(G, pos, nodelist = paths,
                               node_size =  list(np.array(sizenodo) * 200),
                               node_shape = 'o',
                               node_color = edge_colors[usercolormap][0:len(paths)],
                               linewidths = 0.5, edgecolors= 'white')
    
        nx.draw_networkx_edges(G, pos, width = 2, alpha= 0.5, edge_color= verti)
        ax.axis('equal')
        
        for nod_entry in links_for_entrys:
            ax.annotate(nod_entry,
                        xy=pos[nod_entry], xycoords='data',
                        url = links_for_entrys[nod_entry][1],
                        xytext=pos[nod_entry], textcoords='data',
                        color = newcolor, #..................................
                        fontweight='bold',
                        size = newsize, #..................................
                        ha='center', va="center",
                        bbox=dict(boxstyle="circle",
                                  pad = newpad, #..................................
                                  alpha=1, fc=colorder[nod_entry],
                        ec="none", url = links_for_entrys[nod_entry][1]))
            
        ax.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax2 = fig.add_axes([0.5, 0.77, 0.015, .18]) # vertical
        
        if len(valor_unico) > 1:
            ax2.axis('off')
            posiciones_leyenda_de_abajo = [.5, 0.2, 0.15, 0.73]
            pass
    
        else:
            norm = mpl.colors.Normalize(vmax = valor_maximo,
                                        vmin = valor_minimo)
            cb1 = mpl.colorbar.ColorbarBase(ax2,
                                            cmap=ListedColormap(mycmap(np.linspace(1, 0, len(zzz['values'])*2))),
                                            norm=norm, spacing='proportional',
                                            orientation='vertical')
            plt.tick_params(axis="y", color="grey")
            cb1.set_label(barcolortitle, labelpad=-30, y=1.2, rotation=0, size=10, fontweight='bold')
            plt.yticks(size = 7)
            cb1.outline.set_linewidth(0)
            
            posiciones_leyenda_de_abajo = [0.5, 0, 0.15, 0.73]
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax3 = fig.add_axes(posiciones_leyenda_de_abajo)
        
        ax3.text(0,0.95, 'KEGG Pathways', size=12,ha='left',color= 'black', fontweight='bold')
        
        n = 0.05
        for t in paths[0:20]:
            n += 0.043
            ax3.annotate(name_term[t],
                        xy=np.array([0,1-n]), xycoords='data',
                        url = url_for_kegg[t],
                        xytext=np.array([0,1-n]), textcoords='data',
                        size=8, ha='left', va="center",
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[t], ec="none",
                                 url = url_for_kegg[t]))
        
        ax3.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
        
        
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
        plt.close()
    
    
    # # Figura 5 y 6
    
    # In[223]:
    
    
    label_gene = {}
    for i, j in zip(keggtabla.label.tolist(), keggtabla.Entry_Kegg.tolist()):
        label_gene[i] = [i, j]
    
        
    
    IMGLABEL = [5, 6]
    
    for tipo_label, imglabel, sizlab in zip([0, 1], IMGLABEL, [valor, valor * 0.7]):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, .5, 1])
        
        nx.draw_networkx_nodes(G, pos, nodelist = paths,
                               node_size =  list(np.array(sizenodo) * 200),
                               node_shape = 'o',
                               node_color = edge_colors[usercolormap][0:len(paths)],
                               linewidths = 0.5, edgecolors= 'white')
    
        nx.draw_networkx_edges(G, pos, width = 2, alpha= 0.5, edge_color= verti)
        ax.axis('equal')
    
        
        for nod_entry in label_gene:
            ax.annotate(label_gene[nod_entry][tipo_label],
                        xy=pos[nod_entry], xycoords='data',
                        url = links_for_entrys[nod_entry][1], #rotation=45,
                        xytext=pos[nod_entry],textcoords='data',
                        color=colorder[nod_entry],
                        #fontweight='bold',
                        size=sizlab,
                        ha='center', va="center",
                        bbox=dict(boxstyle="Square",
                                  pad=0.1,
                                  alpha=0.1, fc='whitesmoke',
                        ec="none", url = links_for_entrys[nod_entry][1]))
        
        
        
            
        ax.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax2 = fig.add_axes([0.5, 0.77, 0.015, .18]) # vertical
        
        if len(valor_unico) > 1:
            ax2.axis('off')
            posiciones_leyenda_de_abajo = [.5, 0.2, 0.15, 0.73]
            pass
    
        else:
            norm = mpl.colors.Normalize(vmax = valor_maximo,
                                        vmin = valor_minimo)
            cb1 = mpl.colorbar.ColorbarBase(ax2,
                                            cmap=ListedColormap(mycmap(np.linspace(1, 0, len(zzz['values'])*2))),
                                            norm=norm, spacing='proportional',
                                            orientation='vertical')
            plt.tick_params(axis="y", color="grey")
            cb1.set_label(barcolortitle, labelpad=-30, y=1.2, rotation=0, size=10, fontweight='bold')
            plt.yticks(size = 7)
            cb1.outline.set_linewidth(0)
            
            posiciones_leyenda_de_abajo = [0.5, 0, 0.15, 0.73]
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax3 = fig.add_axes(posiciones_leyenda_de_abajo)
        
        ax3.text(0,0.95, 'KEGG Pathways', size=12,ha='left',color= 'black', fontweight='bold')
        
        n = 0.05
        for t in paths[0:20]:
            n += 0.043
            ax3.annotate(name_term[t],
                        xy=np.array([0,1-n]), xycoords='data',
                        url = url_for_kegg[t],
                        xytext=np.array([0,1-n]), textcoords='data',
                        size=8, ha='left', va="center",
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[t], ec="none",
                                 url = url_for_kegg[t]))
        
        ax3.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
        
        
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
        
        plt.close()
    
    
    # In[ ]:
    
    
    
    
    
    # In[236]:
    
    
    columnas = ['Path','LogminP', 'LogminFDR', 'Entry']
    columnas
    
    
    # In[237]:
    
    
    frame1 = DataFrame(keggtabla[columnas].drop_duplicates()          .groupby(['Path', 'LogminP', 'LogminFDR']).Entry.count()).reset_index()
    frame1['constante'] = -np.log10(0.05)
    frame1 = frame1.sort_values(by ='LogminP',ascending=False).reset_index(drop=True)
    frame1
    
    
    # In[238]:
    
    
    col_font_annotate = {}
    for i, j in zip(paths, np.repeat('black', len(paths))):
        col_font_annotate[i] = j
    col_font_annotate
    
    
    # In[239]:
    
    
    # solo se mostraran unicamente las 20 primeras
    if len(frame1) <= 20:
        faltantes = 20 - len(frame1)
        barras_vacias = []
        add_to_labnodeterms = {}
        links_perdidos = {}
        for i in range(faltantes):
            barras_vacias.append([str(i),None, None,None,-np.log10(0.05)])
            labnodeterms.update({str(i):'white'}) # actualizo el dict anterior
            add_to_labnodeterms[str(i)] = 'white'
            links_perdidos[str(i)] = None
            col_font_annotate.update({str(i):'white'}) # actualizo el dict anterior
            empty = DataFrame(barras_vacias, columns = ['Path','LogminP','LogminFDR','Entry', 'constante']) 
        
        
        
        frame2 = frame1.sort_values(by ='LogminP',ascending=True).reset_index(drop=True)
        frame3 = pd.concat([empty, frame2])
        url_for_kegg.update(links_perdidos)
        name_term.update(links_perdidos)
    if len(frame1) >= 21:
        frame2 = frame1.iloc[0:20]
        frame3 = frame2.sort_values(by ='LogminP',ascending=True).reset_index(drop=True)
   
    # una forma de reestablecer este dict
    #for i in add_to_labnodeterms:
    #    labnodeterms.pop(i)
    
    
    # In[243]:
    
    
    def bar_parameters(df = DataFrame([]), colum_val = 1, column_lab = 0):
        ejey = list(df.iloc[0:len(df),colum_val])
        ejex = list(df.iloc[0:len(df),column_lab])
        return ejex, ejey 
    
    
    # In[244]:
    
    
    valores_constantes = frame3.constante.tolist()
    logaritmo_fdr = frame3.LogminFDR.tolist()
    cuentas_entry = frame3.Entry.tolist()
    
    
    # In[ ]:
    
    
    
    
    
    # # Figura 7 y 8
    
    # In[245]:
    
    
    IMGLABEL = [7, 8]
    for newsize, newpad, newcolor, imglabel in zip(NEWSIZE, NEWPAD, NEWCOLOR, IMGLABEL):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, .5, 1])
        
        nx.draw_networkx_nodes(G, pos, nodelist = paths,
                               node_size =  list(np.array(sizenodo) * 200),
                               node_shape = 'o',
                               node_color = edge_colors[usercolormap][0:len(paths)],
                               linewidths = 0.5, edgecolors= 'white')
    
        nx.draw_networkx_edges(G, pos, width = 2, alpha= 0.5, edge_color= verti)
        ax.axis('equal')
        
        for nod_entry in links_for_entrys:
            ax.annotate(nod_entry,
                        xy=pos[nod_entry], xycoords='data',
                        url = links_for_entrys[nod_entry][1],
                        xytext=pos[nod_entry], textcoords='data',
                        color = newcolor, #..................................
                        fontweight='bold',
                        size = newsize, #..................................
                        ha='center', va="center",
                        bbox=dict(boxstyle="circle",
                                  pad = newpad, #..................................
                                  alpha=1, fc=colorder[nod_entry],
                        ec="none", url = links_for_entrys[nod_entry][1]))
            
        ax.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax2 = fig.add_axes([0.5, 0.77, 0.015, .18]) # vertical
        
        if len(valor_unico) > 1:
            ax2.axis('off')
            posiciones_leyenda_de_abajo = [.5, 0.2, 0.5, 0.7]
            pass
    
        else:
            norm = mpl.colors.Normalize(vmax = valor_maximo,
                                        vmin = valor_minimo)
            cb1 = mpl.colorbar.ColorbarBase(ax2,
                                            cmap=ListedColormap(mycmap(np.linspace(1, 0, len(zzz['values'])*2))),
                                            norm=norm, spacing='proportional',
                                            orientation='vertical')
            plt.tick_params(axis="y", color="grey")
            cb1.set_label(barcolortitle, labelpad=-30, y=1.2, rotation=0, size=10, fontweight='bold')
            plt.yticks(size = 7)
            cb1.outline.set_linewidth(0)
            
            posiciones_leyenda_de_abajo = [0.5, 0, 0.5, 0.7]
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax3 = fig.add_axes(posiciones_leyenda_de_abajo)
        
        ejes = bar_parameters(df = frame3)
            
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    
        h = 0
        for i, j, k in zip(ejes[0] , ejes[1], cuentas_entry):
            # barras
            ax3.barh(i, j, height= 0.8,
                        color= labnodeterms[i],
                        align='center',
                        linewidth = 0,
                        alpha = 1)
            # valores
            ax3.annotate(' '+str(k),
                        xy= np.array([0 + j, h]), xycoords='data', #(rect.get_x() + rect.get_width() / 2, height),
                        xytext=np.array([0 + j, h]), textcoords='data',  # 3 points vertical offset
                        size=7, ha='left', va='center')
            ax3.annotate(i,
                        xy=np.array([-0.2, h]), xycoords='data',
                        url = url_for_kegg[i], color = col_font_annotate[i],
                        xytext=np.array([-0.2, h]), textcoords='data',
                        size=7, ha='right', va="center",# fontweight='bold',
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[i],
                                  url = url_for_kegg[i], ec="none"))
            
            h +=1
    
        ax3.plot(np.array(valores_constantes), np.array(list(range(0,20))),
                 color = 'white',linewidth=0.5, zorder=1)
        ax3.scatter(np.array(valores_constantes), np.array(list(range(0,20))),
                 s=5,c='white', marker='s', zorder=2)
        ax3.plot(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 color = 'blue',linewidth=0.5)
        ax3.scatter(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 s=5,c='blue', marker='o', zorder=2)
        
        plt.text(0, 21, "       KEGG Pathways\n\nP-value & Ajusted P-value",size= 8,fontweight='bold')
        # configuracionmanual de la escala
        scala = list(range(0,30+1, 2)) # escala prdeterminada, independiente de los valores
        
        q=0
        w=2
        for i in range(len(scala)): # este bucle encuentra el valor maximo dentro de la escala predeterminada
            if q <= frame2.LogminP.max() <= w:
                val_max_to_scala = w
            if w == 30:
                break
            q+=2
            w+=2
        
        
        ax3.plot([0, val_max_to_scala], [19.8, 19.8], color = 'black',linewidth=0.3, zorder=1)
        
        
        for i in list(range(0, val_max_to_scala+1,2)):
            plt.text(i, 20.2, str(i), size= 6, ha='center')
            plt.text(i, 19.8, '|', size= 4, ha='center')
        
    
        plt.xticks([-2.5] + scala, size=6, color='black') #fontweight='bold'
        
        
        plt.yticks(color='none') # oculta las etiquetas del eje y
        
        
        ax3.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
        
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
        plt.close()
    
    
    # In[ ]:
    
    
    
    
    
    # # Figura 9 y 10
    
    # #  <font color = red>agregar otras redes donde se etiqueten a los genes/proteinas<font>
    
    # In[262]:
    
    
    IMGLABEL = [9, 10]
    for tipo_label, newsize, newpad, newcolor, imglabel in zip([0, 1], [valor, valor * 0.6], [0.1, 0.1],
                                                               ['white', 'white'], IMGLABEL):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, .5, 1])
        
        nx.draw_networkx_nodes(G, pos, nodelist = paths,
                               node_size =  list(np.array(sizenodo) * 200),
                               node_shape = 'o',
                               node_color = edge_colors[usercolormap][0:len(paths)],
                               linewidths = 0.5, edgecolors= 'white')
    
        nx.draw_networkx_edges(G, pos, width = 2, alpha= 0.5, edge_color= verti)
        ax.axis('equal')
        
        for nod_entry in links_for_entrys:
            ax.annotate(label_gene[nod_entry][tipo_label],
                        xy=pos[nod_entry], xycoords='data',
                        url = links_for_entrys[nod_entry][1],
                        xytext=pos[nod_entry], textcoords='data',
                        color = newcolor, #..................................
                        #fontweight='bold',
                        size = newsize, #..................................
                        ha='center', va="center",
                        bbox=dict(boxstyle="round",
                                  #pad = newpad, #..................................
                                  alpha=1, fc=colorder[nod_entry],
                        ec="none", url = links_for_entrys[nod_entry][1]))
            
        ax.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax2 = fig.add_axes([0.5, 0.77, 0.015, .18]) # vertical
        
        if len(valor_unico) > 1:
            ax2.axis('off')
            posiciones_leyenda_de_abajo = [.5, 0.2, 0.5, 0.7]
            pass
    
        else:
            norm = mpl.colors.Normalize(vmax = valor_maximo,
                                        vmin = valor_minimo)
            cb1 = mpl.colorbar.ColorbarBase(ax2,
                                            cmap=ListedColormap(mycmap(np.linspace(1, 0, len(zzz['values'])*2))),
                                            norm=norm, spacing='proportional',
                                            orientation='vertical')
            plt.tick_params(axis="y", color="grey")
            cb1.set_label(barcolortitle, labelpad=-30, y=1.2, rotation=0, size=10, fontweight='bold')
            plt.yticks(size = 7)
            cb1.outline.set_linewidth(0)
            
            posiciones_leyenda_de_abajo = [0.5, 0, 0.5, 0.7]
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
        ax3 = fig.add_axes(posiciones_leyenda_de_abajo)
        
        ejes = bar_parameters(df = frame3)
            
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    
        h = 0
        for i, j, k in zip(ejes[0] , ejes[1], cuentas_entry):
            # barras
            ax3.barh(i, j, height= 0.8,
                        color= labnodeterms[i],
                        align='center',
                        linewidth = 0,
                        alpha = 1)
            # valores
            ax3.annotate(' '+str(k),
                        xy= np.array([0 + j, h]), xycoords='data', #(rect.get_x() + rect.get_width() / 2, height),
                        xytext=np.array([0 + j, h]), textcoords='data',  # 3 points vertical offset
                        size=7, ha='left', va='center')
            ax3.annotate(i,
                        xy=np.array([-0.2, h]), xycoords='data',
                        url = url_for_kegg[i], color = col_font_annotate[i],
                        xytext=np.array([-0.2, h]), textcoords='data',
                        size=7, ha='right', va="center",# fontweight='bold',
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[i],
                                  url = url_for_kegg[i], ec="none"))
            
            h +=1
    
        ax3.plot(np.array(valores_constantes), np.array(list(range(0,20))),
                 color = 'white',linewidth=0.5, zorder=1)
        ax3.scatter(np.array(valores_constantes), np.array(list(range(0,20))),
                 s=5,c='white', marker='s', zorder=2)
        ax3.plot(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 color = 'blue',linewidth=0.5)
        ax3.scatter(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 s=5,c='blue', marker='o', zorder=2)
        
        plt.text(0, 21, "       KEGG Pathways\n\nP-value & Ajusted P-value",size= 8,fontweight='bold')
        # configuracionmanual de la escala
        scala = list(range(0,30+1, 2)) # escala prdeterminada, independiente de los valores
        
        q=0
        w=2
        for i in range(len(scala)): # este bucle encuentra el valor maximo dentro de la escala predeterminada
            if q <= frame2.LogminP.max() <= w:
                val_max_to_scala = w
            if w == 30:
                break
            q+=2
            w+=2
        
        
        ax3.plot([0, val_max_to_scala], [19.8, 19.8], color = 'black',linewidth=0.3, zorder=1)
        
        
        for i in list(range(0, val_max_to_scala+1,2)):
            plt.text(i, 20.2, str(i), size= 6, ha='center')
            plt.text(i, 19.8, '|', size= 4, ha='center')
        
    
        plt.xticks([-2.5] + scala, size=6, color='black') #fontweight='bold'
        
        
        plt.yticks(color='none') # oculta las etiquetas del eje y
        
        
        ax3.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
        
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
        plt.close()
    
    
    
    ###########################################################################
    
    
    # In[263]:
    
    
    
    IMGLABEL = [11, 12]
    for lll, imglabel in zip([0, 1], IMGLABEL):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, 1, 1])
        
        ejes = bar_parameters(df = frame3)
            
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    
        h = 0
        for i, j, k in zip(ejes[0] , ejes[1], cuentas_entry):
            # barras
            ax.barh(i, j, height= 0.8,
                        color= labnodeterms[i],
                        align='center',
                        linewidth = 0,
                        alpha = 1)
            # valores
            ax.annotate(' '+str(k),
                        xy= np.array([0 + j, h]), xycoords='data', #(rect.get_x() + rect.get_width() / 2, height),
                        xytext=np.array([0 + j, h]), textcoords='data',  # 3 points vertical offset
                        size=9, ha='left', va='center')
            if lll == 0:
                ax.annotate(name_term[i],
                            xy=np.array([-0.2, h]), # espacio entre la barra y el path term
                            xycoords='data',
                            url = url_for_kegg[i], color = col_font_annotate[i],
                            xytext=np.array([-0.2, h]), textcoords='data',
                            size=9, ha='right', va="center",# fontweight='bold',
                            bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[i],
                                      url = url_for_kegg[i], ec="none"))
            if lll == 1:
                ax.annotate(name_term[i],
                            xy=np.array([-0.2, h]), # espacio entre la barra y el path term
                            xycoords='data',
                            url = url_for_kegg[i], color = col_font_annotate[i],
                            xytext=np.array([-0.2, h]), textcoords='data',
                            size=9, ha='right', va="center",# fontweight='bold',
                            bbox=dict(boxstyle="round", alpha=0.1, fc='whitesmoke',
                                      url = url_for_kegg[i], ec="none"))
            h +=1
    
        ax.plot(np.array(valores_constantes), np.array(list(range(0,20))),
                 color = 'white',linewidth=0.5, zorder=1)
        ax.scatter(np.array(valores_constantes), np.array(list(range(0,20))),
                 s=5,c='white', marker='s', zorder=2)
        ax.plot(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 color = 'blue',linewidth=0.5)
        ax.scatter(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 s=5,c='blue', marker='o', zorder=2)
        
        
        plt.text(0, 21, "       KEGG Pathways\n\nP-value & Ajusted P-value",size= 9,fontweight='bold')
        # configuracionmanual de la escala
        scala = list(range(0,30+1, 2)) # escala prdeterminada, independiente de los valores
        
        q=0
        w=2
        for i in range(len(scala)): # este bucle encuentra el valor maximo dentro de la escala predeterminada
            if q <= frame2.LogminP.max() <= w:
                val_max_to_scala = w
            if w == 30:
                break
            q+=2
            w+=2
        
        
        ax.plot([0, val_max_to_scala], [19.8, 19.8], color = 'black',linewidth=0.3, zorder=1)
        
        
        for i in list(range(0, val_max_to_scala+1,2)):
            plt.text(i, 20.2, str(i), size= 6, ha='center')
            plt.text(i, 19.8, '|', size= 4, ha='center')
        
    
        plt.xticks([-15] + scala, size=6, color='black') #fontweight='bold'
        
        plt.yticks(list(range(1,25)), color='none') # aumento el margen superior para mostrar el titulo
        #plt.yticks(color='none') # oculta las etiquetas del eje y
        
        
        ax.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
           
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Bar_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Bar_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
    
        plt.close()
    
    
    
    ##################################################################################################
    
    
    
    frame4 = frame3.sort_values(by ='Entry',ascending=False).reset_index(drop=True)
    
    
    # In[269]:
    
    
    for_pie = bar_parameters(df = frame4.dropna(), colum_val = 3, column_lab = 0)
    
    
    
    ######
    if len(paths) < 6:
        sizepielabel = 20
    if 6 <= len(paths) < 8:
        sizepielabel = 18
    if len(paths) == 8:
        sizepielabel = 16
    if 9 <= len(paths) < 11:
        sizepielabel = 14
    if 11 <= len(paths) < 13:
        sizepielabel = 12
    if 13 <= len(paths) < 15:
        sizepielabel = 10
    if len(paths) >= 15:
        sizepielabel = 9
    ####
    
    
    # In[ ]:
    
    
    
    
    
    # In[316]:
    
    
    IMGLABEL = [13, 14]
    for abrir, imglabel in zip([0, 0.5], IMGLABEL):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_axes([0, 0, 1, 1])
    
        ax.pie(for_pie[1], colors = edge_colors[usercolormap][0:len(paths)])
    
        n = 0.9
        for i, j, k in zip(for_pie[0], for_pie[1], edge_colors[usercolormap][0:len(paths)]):
            ax.annotate(name_term[i]+'   ('+str(j)+')',
                        xy=(1.1, n),
                        xytext=(1.1, n),
                        url = url_for_kegg[i],
                        size=sizepielabel,
                        bbox=dict(boxstyle="round", alpha=0.5, fc=k, ec="k", lw=0.72,
                                 url = url_for_kegg[i]))
            n -= sizepielabel * .01 # una décima parte
    
        centre_circle = plt.Circle((0,0),abrir,fc='white', alpha = 1)
        plt.gca().add_artist(centre_circle)
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.3, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
    
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Circle_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Circle_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
          
        plt.close()
    
    
    # In[ ]:
    
    
    
    
    
    # In[317]:
    
    
    short_name_term = dict(zip(paths, results_process_P.Short_Term.drop_duplicates().tolist()))

    
    
    # kt = kegg tabla
    kt = keggtabla[['Path', 'Entry_Kegg']]
    
    df1 = kt.merge(kt, on = 'Entry_Kegg', how = 'left').drop_duplicates()
    matrix = df1.pivot_table(values='Entry_Kegg',index=['Path_x'],aggfunc=len,columns=['Path_y'])
    
    ###
    df_mat = []
    n = -1
    for i in list(matrix.columns.values):
        n += 1
        new = DataFrame(matrix.iloc[n:len(matrix)][i])
        nn = -1
        for index, row in new.iterrows():
            nn += 1
            df_mat.append([index, i, new.iloc[nn][i]])
        nn = 0
    ###
    df_mat = DataFrame(df_mat, columns = ['go0', 'go1', 'val']).dropna()
    ###
    nodos = []
    for index, row in df_mat.iterrows():
        if row.go0 == row.go1:
            #print(row.go0, row.go1)
            continue
        else:
            #print(row.go0, row.go1)
            nodos.append([row.go0, row.go1, row.val])
    nodos = DataFrame(nodos)
    if len(nodos) > 1:
        nodos = DataFrame([[i for i in nodos[0]], [i for i in nodos[1]], nodos[2]]).T
    else:
        pass
    #### >>>>>>>>>>>>>>>>
    # si interacciona con mas uno, eliminar la redundancia, y si no interacciona con ninguno, dejar el nodo
    # y su valor, este se verá en la red como un nodo aislado
    aislado = [i for i in matrix.columns if len(matrix[[i]].dropna()) == 1]
    aislado = [df_mat[df_mat.go0 == i] for i in aislado]
    if len(aislado) > 0:
        aislado = pd.concat(aislado)
        aislado.columns = [0, 1, 2]
        aislado = DataFrame([[i for i in aislado[0]], [i for i in aislado[1]], aislado[2]]).T
        nodos = pd.concat([nodos, aislado])
    else:
        pass
    
    g = nx.Graph()
    for index, row in nodos.iterrows():
        g.add_edge(row[0], row[1],weight = int(row[2]))
    
    ed = [(u,v,d['weight']) for (u,v,d) in g.edges(data=True) if d['weight'] >= 1]
    
    
    # In[ ]:
    
    arc_weight = nx.get_edge_attributes(g,'weight')
    
    NEWCOLOR = ['black', colorletra]
    IMGLABEL = [15, 16]
    for newcolor, imglabel, size_la in zip(NEWCOLOR, IMGLABEL, [sizepielabel * 0.35, sizepielabel * 0.35]):
        acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        
        pos = nx.kamada_kawai_layout(g, dist=None, weight='weight', scale=1, center=None, dim=2)
    
        fig = plt.figure(figsize=(15, 7))
    
        ax = fig.add_axes([0, 0, .5, 1])
    
        nx.draw_networkx_nodes(g,pos,node_list = labnodeterms,
                               node_color= [labnodeterms[i] for i in g.nodes()],
                               alpha= 1,
                               node_size = list(np.array(sizenodo) * 200),
                               zorder = 2)
        nx.draw_networkx_edges(g,pos, edgelist=ed,width = np.array([i[2] for i in ed]),
                                   alpha= 0.6,edge_color=  'grey',style='-')
        
        nx.draw_networkx_edge_labels(g, pos, edge_color= 'grey',
                                     font_size=size_la, edge_labels=arc_weight)
        
        for nod_term in list(g.nodes()):
            ax.annotate(short_name_term[nod_term],
                        xy=pos[nod_term], xycoords='data',
                        #url = url_for_kegg[nod_term], #rotation=45,
                        xytext=pos[nod_term],textcoords='data',
                        color= newcolor,
                        size=size_la,
                        ha='center', va="center",
                        bbox=dict(boxstyle="circle",
                                  pad=0, #url = url_for_kegg[nod_term]
                                  alpha=0.01, fc='whitesmoke',
                                  ec="none", ))
    
        ax.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        ax3 = fig.add_axes([0.5, 0.2, 0.5, 0.7])
        
        ejes = bar_parameters(df = frame3)
            
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    
        h = 0
        for i, j, k in zip(ejes[0] , ejes[1], cuentas_entry):
            # barras
            ax3.barh(i, j, height= 0.8,
                        color= labnodeterms[i],
                        align='center',
                        linewidth = 0,
                        alpha = 1)
            # valores
            ax3.annotate(' '+str(k),
                        xy= np.array([0 + j, h]), xycoords='data', #(rect.get_x() + rect.get_width() / 2, height),
                        xytext=np.array([0 + j, h]), textcoords='data',  # 3 points vertical offset
                        size=7, ha='left', va='center')
            ax3.annotate(i,
                        xy=np.array([-0.2, h]), xycoords='data',
                        url = url_for_kegg[i], color = col_font_annotate[i],
                        xytext=np.array([-0.2, h]), textcoords='data',
                        size=7, ha='right', va="center",# fontweight='bold',
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[i],
                                  url = url_for_kegg[i], ec="none"))
            
            h +=1
    
        ax3.plot(np.array(valores_constantes), np.array(list(range(0,20))),
                 color = 'white',linewidth=0.5, zorder=1)
        ax3.scatter(np.array(valores_constantes), np.array(list(range(0,20))),
                 s=5,c='white', marker='s', zorder=2)
        ax3.plot(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 color = 'blue',linewidth=0.5)
        ax3.scatter(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 s=5,c='blue', marker='o', zorder=2)
        
        plt.text(0, 21, "       KEGG Pathways\n\nP-value & Ajusted P-value",size= 8,fontweight='bold')
        # configuracionmanual de la escala
        scala = list(range(0,30+1, 2)) # escala prdeterminada, independiente de los valores
        
        q=0
        w=2
        for i in range(len(scala)): # este bucle encuentra el valor maximo dentro de la escala predeterminada
            if q <= frame2.LogminP.max() <= w:
                val_max_to_scala = w
            if w == 30:
                break
            q+=2
            w+=2
        
        
        ax3.plot([0, val_max_to_scala], [19.8, 19.8], color = 'black',linewidth=0.3, zorder=1)
        
        
        for i in list(range(0, val_max_to_scala+1,2)):
            plt.text(i, 20.2, str(i), size= 6, ha='center')
            plt.text(i, 19.8, '|', size= 4, ha='center')
        
    
        plt.xticks([-2.5] + scala, size=6, color='black') #fontweight='bold'
        
        
        plt.yticks(color='none') # oculta las etiquetas del eje y
        
        
        ax3.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
    
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
      
    
        plt.close()
    
    
    # In[ ]:
    
    
    arc_weight = nx.get_edge_attributes(g,'weight')
    
    NEWCOLOR = ['black', colorletra]
    IMGLABEL = [17, 18]
    for newcolor, imglabel, size_la in zip(NEWCOLOR, IMGLABEL, [sizepielabel * 0.35, sizepielabel * 0.35]):
        if imglabel == 18:
            acumulacion += str(imglabel)+'.'
        else: 
            acumulacion += str(imglabel)+', '
        sys.stdout.write("\rPlots: %s" % (acumulacion))   
        sys.stdout.flush()
        
        pos = nx.circular_layout(g)
    
        fig = plt.figure(figsize=(15, 7))
    
        ax = fig.add_axes([0, 0, .5, 1])
    
        nx.draw_networkx_nodes(g,pos,node_list = paths,
                               node_color= [labnodeterms[i] for i in g.nodes()],
                               alpha= 1,
                               node_size = list(np.array(sizenodo) * 200),
                               zorder = 2)
        nx.draw_networkx_edges(g,pos, edgelist=ed,width = np.array([i[2] for i in ed]),
                                   alpha= 0.6,edge_color=  'grey',style='-')
        
        nx.draw_networkx_edge_labels(g, pos, edge_color= 'grey',
                                     font_size=size_la, edge_labels=arc_weight)
        
        for nod_term in list(g.nodes()):
            ax.annotate(short_name_term[nod_term],
                        xy=pos[nod_term], xycoords='data',
                        #url = url_for_kegg[nod_term], #rotation=45,
                        xytext=pos[nod_term],textcoords='data',
                        color= newcolor,
                        size=size_la,
                        ha='center', va="center",
                        bbox=dict(boxstyle="circle",
                                  pad=0, #url = url_for_kegg[nod_term]
                                  alpha=0.01, fc='whitesmoke',
                                  ec="none", ))
    
    
        ax.axis('off')
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        ax3 = fig.add_axes([0.5, 0.2, 0.5, 0.7])
        
        ejes = bar_parameters(df = frame3)
            
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    
        h = 0
        for i, j, k in zip(ejes[0] , ejes[1], cuentas_entry):
            # barras
            ax3.barh(i, j, height= 0.8,
                        color= labnodeterms[i],
                        align='center',
                        linewidth = 0,
                        alpha = 1)
            # valores
            ax3.annotate(' '+str(k),
                        xy= np.array([0 + j, h]), xycoords='data', #(rect.get_x() + rect.get_width() / 2, height),
                        xytext=np.array([0 + j, h]), textcoords='data',  # 3 points vertical offset
                        size=7, ha='left', va='center')
            ax3.annotate(i,
                        xy=np.array([-0.2, h]), xycoords='data',
                        url = url_for_kegg[i], color = col_font_annotate[i],
                        xytext=np.array([-0.2, h]), textcoords='data',
                        size=7, ha='right', va="center",# fontweight='bold',
                        bbox=dict(boxstyle="round", alpha=0.5, fc=labnodeterms[i],
                                  url = url_for_kegg[i], ec="none"))
            
            h +=1
    
        ax3.plot(np.array(valores_constantes), np.array(list(range(0,20))),
                 color = 'white',linewidth=0.5, zorder=1)
        ax3.scatter(np.array(valores_constantes), np.array(list(range(0,20))),
                 s=5,c='white', marker='s', zorder=2)
        ax3.plot(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 color = 'blue',linewidth=0.5)
        ax3.scatter(np.array(logaritmo_fdr), np.array(list(range(0,20))),
                 s=5,c='blue', marker='o', zorder=2)
        
        plt.text(0, 21, "       KEGG Pathways\n\nP-value & Ajusted P-value",size= 8,fontweight='bold')
        # configuracionmanual de la escala
        scala = list(range(0,30+1, 2)) # escala prdeterminada, independiente de los valores
        
        q=0
        w=2
        for i in range(len(scala)): # este bucle encuentra el valor maximo dentro de la escala predeterminada
            if q <= frame2.LogminP.max() <= w:
                val_max_to_scala = w
            if w == 30:
                break
            q+=2
            w+=2
        
        
        ax3.plot([0, val_max_to_scala], [19.8, 19.8], color = 'black',linewidth=0.3, zorder=1)
        
        
        for i in list(range(0, val_max_to_scala+1,2)):
            plt.text(i, 20.2, str(i), size= 6, ha='center')
            plt.text(i, 19.8, '|', size= 4, ha='center')
        
    
        plt.xticks([-2.5] + scala, size=6, color='black') #fontweight='bold'
        
        
        plt.yticks(color='none') # oculta las etiquetas del eje y
        
        
        ax3.axis('off')
        
        #■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
    
        captation = fig.add_axes([.39, 0, 0.11, 0.02])
    
        captation.text(0,0.3, 'NeVOmics {}'.format(datetime.now()).split('.')[0], size=7,ha='left',color= 'black')
    
        captation.axis('off')
    
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.svg', bbox_inches='tight')
        plt.savefig(new_dir_plots+'/NeVOmics_Plot_Network_'+str(imglabel)+'.png', dpi = 900, bbox_inches='tight')
    
        plt.close()

if createnetworks == '0':
    pass
    



# In[348]:


# crear circos solo si fue elegido por el usuario
if createcircos == '1':
    ## Create a folder
    new_dir_plots = "job_KEGG_plots"
    os.makedirs(new_dir_plots,exist_ok=True)
    print('\n\n-----------------------')
    print('Building graphics with R ...')
    print('Wait...')
    # localizacion de la libreria
    R_lib = open('../NeVOmics_locRlib.txt', 'r')
    R_lib = R_lib.read()
    r_script = requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
    R_script_enrich = re.sub('rliblocation', R_lib, r_script)
    fr = open('Kegg_Enrichment_Plots.R','w')
    fr.write(R_script_enrich)
    fr.close()
    ######
    R_exe = open('../NeVOmics_locRexe.txt', 'r')
    R_exe = R_exe.read()
    run_uni = subprocess.Popen([R_exe, 'CMD', 'BATCH', '--no-save', 'Kegg_Enrichment_Plots.R'])
    run_uni.wait()
    lapso_total = datetime.now() - inicio_total
    
if createcircos == '0': # el usuario decició no crear estos gráficos
    pass


del_stop_process()



