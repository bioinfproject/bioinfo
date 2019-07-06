#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json
import pandas as pd
from pandas import DataFrame
import urllib.request
from urllib.request import urlopen
import requests
import re
import sys
import os
import numpy as np
# Network
from networkx import path_graph, random_layout
import matplotlib.pyplot as plt
import networkx as nx
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns


# In[ ]:


with urllib.request.urlopen('ftp://ftp.pantherdb.org/pathway/current_release/') as response:
    info = response.read()


# In[ ]:


file = ''.join(re.findall('SequenceAssociation.*txt', info.decode()))
print('\nPanther Association Pathway available in http://www.pantherdb.org/')
print(file)

# In[ ]:


url = 'ftp://ftp.pantherdb.org/pathway/current_release/'+file


# In[ ]:


# salvar en variable
#import urllib.request
#with urllib.request.urlopen(url) as response:
#    mapping = response.read().decode()
#mapping = mapping


# In[ ]:


if os.path.exists(file):
    print('\nThe "Panther Association Pathway" file already exists')
    print(file)
else:
    #https://docs.python.org/3/howto/urllib2.html
    response = urllib.request.urlopen(url)
    total_length = urllib.request.urlopen(url).headers.get('Content-length')
    print('\nWait while downloading "Panther Association Pathway" file: ',
          round(int(total_length)/1000000,2), 'MB\n...')
    with open(file, 'wb') as f:
        for data in response:
            f.write(data)
    print('Downloaded file:', file)
    print('Finished process')


# In[ ]:


association = pd.read_csv(file, sep = '\t', header = None)
association = association[[0 , 1]].drop_duplicates()
association.columns = ['Panther_id', 'label']


# In[ ]:


import tkinter as tk
from tkinter import filedialog
## Control of input file
print('\nSubmit file .json (Only PANTHER Pathways)')
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
root.destroy()
if file_path == ():
    print('\n!!!!!!! File not found !!!!!!!')
    import tkinter as tk
    from tkinter import filedialog
    ## Control of input file
    print('\nSubmit file .json')
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    root.destroy()
    if file_path == ():
        print('\n!!!!!!! File not found !!!!!!!')
        sys.exit()
    else:
        print('=====> : ',file_path)
else:
    print('=====> : ',file_path)


# In[ ]:





# In[ ]:


# read file
with open(file_path, 'r') as myfile:
    data=myfile.read()

# parse file
obj = json.loads(data)


# In[ ]:

print('\nInformation')
organism = obj['overrepresentation']['upload_lists']['input_list']['organism']
print('\norganism:           ', organism)
tool_release_date = obj['overrepresentation']['tool_release_date']
print('tool_release_date:  ', tool_release_date)
data_version_release_date = obj['overrepresentation']['data_version_release_date']
print('version_release:    ', data_version_release_date)
annotation_type = obj['overrepresentation']['annotation_type']
print('annotation_type:    ', annotation_type)
correction = obj['overrepresentation']['correction']
print('correction:         ', correction)
test_type = obj['overrepresentation']['test_type']
print('test_type:          ', test_type)


# In[ ]:


# Fisher's Exact
panther_results = []
for i in obj['overrepresentation']['group']:
    x = i['result']['input_list']['number_in_list']
    if x > 0:
        xx = i['result']['input_list']['number_in_list']
        number_in_list = i['result']['input_list']['number_in_list']
        #print(number_in_list)
        fold_enrichment = i['result']['input_list']['fold_enrichment']
        #print(fold_enrichment)
        expected = i['result']['input_list']['expected']
        #print(expected)
        pValue = i['result']['input_list']['pValue']
        #print(pValue)
        fdr = i['result']['input_list']['fdr']
        #print(fdr)
        number_in_reference = i['result']['number_in_reference']
        #print(number_in_reference)
        label = i['result']['term']['label']
        #print(label)
        ##
        y = i['result']['input_list']
        mapped_id = y['mapped_id_list']['mapped_id']
        #print(mapped_id)
        panther_results.append([label, number_in_reference, number_in_list,
                                fold_enrichment, expected, pValue, fdr, mapped_id])
    else:
        continue
###
names = ['label', 'number_in_reference', 'number_in_list', 'fold_enrich',
         'expected', 'pValue', 'fdr', 'Entry']
table = DataFrame(panther_results, columns = names)
table = table.merge(association, on = 'label', how = 'left').fillna('UNCLASSIFIED')

###

###
#options = ['yes', 'y', 'no', 'n']
#while True:
#    cont = input('\nWould you like include "UNCLASSIFIED" Genes? [yes/y   no/n]\n=====> : ').lower()
#    if cont in options:
#        break
cont = 'no'
###
no = ['no', 'n']
yes = ['yes', 'y']
if cont in no:
    table = table[table.label != 'UNCLASSIFIED']
if cont in no:
    table = table 
###
opciones = ['pValue', 'fdr']
while True:
    sorting = input('\nSort by [pValue/fdr]\n=====> : ')
    if sorting in opciones:
        break
table = table.sort_values(by=sorting, ascending=True).reset_index(drop = True)
table.to_csv('table0.csv',index=None)
table0 = table.sort_values(by=sorting, ascending=True).reset_index(drop = True)
table0.to_csv('table1.csv',index=None)
###
###
# comandos para seleccionar los mejores hits, o un intervalo de 
# datos (filas)
hits = ['fixed','f' , 'interval', 'i']
while True:
    inter = input('\nSelect an option [fixed/f   interval/i]\n=====> : ').lower()
    if inter in hits:
        break
###
fix = ['fixed','f']
if inter in fix:
    intervalo = list(range(1,100001))
    while True:
        interv = input('\nSelect a number [10]\n=====> : ')
        if int(interv) in intervalo:
            break
    table = table.iloc[0:int(interv)]
###
interval = ['interval', 'i']
if inter in interval:
    intervalo = list(range(1,100001))
    while True:
        interv = input('\nSelect a range [5-10]\n=====> : ')
        mini = int(interv.split('-')[0])
        maxi = int(interv.split('-')[1])
        if maxi > mini:
            if mini & maxi in intervalo:
                break
            break
    table = table.iloc[mini - 1:maxi]
###
panther_results1 = []
for index, row in table.iterrows():
    if type(row.Entry) == str:
        panther_results1.append([row.Panther_id, row.label, row.number_in_reference, row.number_in_list,
              row.fold_enrich, row.expected, row.pValue, row.fdr, row.Entry])
    else:
        for j in row.Entry:
            panther_results1.append([row.Panther_id, row.label, row.number_in_reference, row.number_in_list,
                  row.fold_enrich, row.expected, row.pValue, row.fdr, j])
names = ['Panther_id', 'label', 'number_in_reference', 'number_in_list', 'fold_enrich',
         'expected', 'pValue', 'fdr', 'Entry']

table = DataFrame(panther_results1, columns = names)
table['Panther'] = 'Panther'
table.to_csv('table1.csv',index=None)


# In[ ]:


table1 = table[['Panther_id', 'Entry', 'Panther', 'label']]


# In[ ]:


options = ['yes', 'y', 'no', 'n']
while True:
    cont = input('\nWould you like include isolated nodes? [yes/y   no/n]\n=====> : ').lower()
    if cont in options:
        break
    ###
no = ['no', 'n']
yes = ['yes', 'y']
#####
if cont in yes:
    respuesta = cont
if cont in no:
    respuesta = cont


# In[ ]:


respuesta


# In[ ]:


def net_plot(df = DataFrame([]), layout = nx.random_layout,label = 'Panther_id',diam_nodos = 10, espe_edges = 0.1, inte = 10,
             color_inter_min = 'k',color_inter_max = 'blue',
             edge_alpha_min = 0.3, edge_alpha_max = 0.3, k_num = 3, color_nodo = 'red', node_alpha = 0.7):
    ##
    # elegir la etiqueta
    if label == 'Panther_id':
        la = 'Panther_id'
        lala = 'label'
    if label == 'label':
        la = 'label'
        lala = 'Panther_id'
        ###
        etiquetas = [] # con este comando fragmento la etiqueta y me quedo con tres palabras unicamente
        for i in df.label:
            if len(i.split(' ')) <= 2:
                etiquetas.append(i)
            else:
                etiquetas.append(''.join(i.split(' ')[0])+'\n'+''.join(i.split(' ')[1])+'...') # +'\n'+i.split(' ')[2]
        df['label'] = etiquetas
        ###
    if label == 'none':
        la = 'Panther_id'
        lala = 'label'
    ##
    df1 = df.groupby(['Panther']).get_group(df.Panther.iloc[0])[[la, 'Entry']].drop_duplicates().merge(df, on = 'Entry', how = 'left').drop_duplicates()
    df2 = DataFrame(df[[la, lala, 'Entry']].drop_duplicates().groupby([la, lala]).Entry.count()).reset_index()
    
    #### >>>>>>>>>>>>>>>>
    #### A partir de una matriz de datos extrae valores no redundantes
    matrix = df1.pivot_table(values='Entry',index=la+'_x',aggfunc=len,columns=la+'_y')
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
    df_mat = DataFrame(df_mat, columns = ['pan0', 'pan1', 'val'])
    ###
    no = ['no', 'n']
    yes = ['yes', 'y']
    #####
    nodos2 = []
    nodos = []
    for index, row in df_mat.iterrows():
        if row.pan0 == row.pan1:
            if respuesta in yes:
                nodos.append([row.pan0, row.pan1, row.val])
            if respuesta in no:
                continue
        else:
            #print(row.pan0, row.pan1)
            nodos.append([row.pan0, row.pan1, row.val])
    #### >>>>>>>>>>>>>>>>
    ###
    nodos = DataFrame(nodos).dropna()
    ####################
    # https://networkx.github.io/documentation/networkx-2.3/auto_examples/drawing/plot_weighted_graph.html#sphx-glr-auto-examples-drawing-plot-weighted-graph-py
    G=nx.Graph()
    for index, row in nodos.iterrows():
        G.add_edge(row[0], row[1],weight = row[2])
    elarge=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] >= inte]
    esmall=[(u,v,d['weight']) for (u,v,d) in G.edges(data=True) if d['weight'] < inte]
    ###
    #circular_layout
    #random_layout
    #shell_layout
    #spring_layout
    #spectral_layout
    pos=layout(G) # positions for all nodes k = k_num
    #print(pos)
    #>>>>>>>>>>>>>>>>>>>>
    #pos = procesamiento(pos)
    #>>>>>>>>>>>>>>>>>>>>
    # nodes
    #------------------------------------------------------------------------------
    # ordenar los valores para representarlos en el tama;o del nodo
    order = []
    for index, row in nodos.iterrows():
        order.append(row[0])
        order.append(row[1])
    orden3 = DataFrame(order).drop_duplicates(keep = 'first').reset_index(drop = True)
    orden3.columns = [la]
    orden4 = pd.merge(orden3, df2, on = la, how = 'left')
    #------------------------------------------------------------------------------
    # https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.drawing.nx_pylab.draw_networkx_edges.html
    nx.draw_networkx_nodes(G,pos,node_size= np.array(orden4.Entry) * diam_nodos,
                           node_color= color_nodo,alpha= node_alpha)
    # edges
    nx.draw_networkx_edges(G,pos,edgelist=esmall,
                           width = np.array([i[2] for i in esmall]) * espe_edges,
                           alpha= edge_alpha_min,edge_color= color_inter_min,style='-')
    nx.draw_networkx_edges(G,pos, edgelist=elarge,
                           width = np.array([i[2] for i in elarge]) * espe_edges,
                           alpha= edge_alpha_max,edge_color= color_inter_max,style= '-')
    # labels
    posicion = {} ## posicion de las etiquetas, ligeramente arriba
    for key, value in pos.items():
        posicion[key] = value + 0.05
    ###
    # arreglo de las posiciones de los nodos en el plano cartesiano
    arr = np.array([[i for i in value] for key, value in pos.items()])
    ###
    if label == 'label':
        nx.draw_networkx_labels(G, posicion, font_size=10, font_weight='bold') # ,font_weight='bold'
        if label == 'label':
            plt.axis([arr[:,0].min() - 0.3, arr[:,0].max() + 0.3,
                      arr[:,1].min() - 0.3, arr[:,1].max() + 0.3])
    if label == 'Panther_id':
        nx.draw_networkx_labels(G, posicion, font_size=10, font_weight='bold') # ,font_weight='bold'
        if label == 'Panther_id':
            plt.axis([arr[:,0].min() - 0.175, arr[:,0].max() + 0.175,
                      arr[:,1].min() - 0.175, arr[:,1].max() + 0.175])
        
    if label == 'none':
        plt.axis([arr[:,0].min() - 0.1, arr[:,0].max() + 0.1,
                  arr[:,1].min() - 0.1, arr[:,1].max() + 0.1])
    ########
    #ax2 = plt.axes([0.62, 0.2, 0.27, 0.3])
    return pos


# In[ ]:


#nx.circular_layout
#nx.random_layout
#nx.shell_layout
#nx.spring_layout
#nx.spectral_layout


# In[ ]:


layouts = [nx.circular_layout, nx.random_layout, nx.shell_layout,
           nx.spring_layout, nx.spectral_layout]
modelos = ['circular_layout', 'random_layout', 'shell_layout',
           'spring_layout', 'spectral_layout']


# In[ ]:

print('')
import warnings
warnings.filterwarnings("ignore")
n = 0
for i, j in zip(layouts, modelos):
    print('Network_Pathways_'+j+'_python.png')
    n += 1
    plt.subplots(figsize=(8,8))
    net_plot(table1, layout = i , diam_nodos = 50, inte = 5,
             # color = plt.cm.viridis(np.linspace(0,1,20))
             color_nodo = 'red', node_alpha = 0.9,
             label = 'label', espe_edges = 3, k_num = 0.2, color_inter_min = 'k',
             color_inter_max = 'b')
    plt.savefig('Network_Pathways_'+j+'_python.png', dpi = 600, bbox_inches='tight')


# In[ ]:


def reject_outliers_2(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return data[s<m]
##
def procesamiento(pos):
    arr = np.array([[i for i in value] for key, value in pos.items()])
    uno0 = list(reject_outliers_2(abs(arr[:,0])))
    uno1 = list(reject_outliers_2(abs(arr[:,1])))
    unos = []
    for i, j in zip(uno0, uno1):
        unos.append([i,j])
    unos = np.array(unos)
    min0 = min(unos[:,0])
    max0 = max(unos[:,0])
    min1 = min(unos[:,1])
    max1 = max(unos[:,1])
    mean_min = np.mean([min0, min1]) + 0.05
    mean_max = np.mean([max0, max1]) + 0.05
    dos0 = list(abs(arr[:,0]))
    out0 = []
    n = 0
    for i in dos0:
        n += 1
        if i in unos[:,0]:
            continue
        else:
            out0.append(i)
    dos1 = list(abs(arr[:,1]))
    out1 = []
    n = 0
    for i in dos1:
        n += 1
        if i in unos[:,1]:
            continue
        else:
            out1.append(i)
    dos = []
    for i, j in zip(out0, out1):
        dos.append([i,j])
    dos = np.array(dos)
    minmax_scale = preprocessing.MinMaxScaler(feature_range=(mean_min, mean_max))
    processing = minmax_scale.fit_transform(dos)
    arreglo = np.concatenate((unos, processing), axis=0)
    keys = [key for key, value, in pos.items()]
    posi = {}
    for i, j in zip(keys, arreglo):
        posi[i] = j
    #return posi


# In[ ]:





# In[ ]:


table2 = table[['Panther_id', 'Entry', 'Panther', 'label', 'pValue', 'fdr']]


def heatmap_plot(df = DataFrame([]), colors = 'Spectral', label_x = 'GO',
                 label_y = 'Term', size_plot = 10, xticks_size = 12,
                 yticks_size = 12, ylabel_size = 12):
    df1 = table2.groupby(['Panther']).get_group(table2.Panther.iloc[0])[['Panther_id', 'Entry']].drop_duplicates().merge(table2, on = 'Entry', how = 'left').drop_duplicates()
    etiquetas = [] # con este comando fragmento la etiqueta y me quedo con tres palabras unicamente
    for i in df1.label:
        if len(i.split(' ')) <= 2:
            etiquetas.append(i)
        else:
            etiquetas.append(' '.join(i.split(' ')[0:2])+' ...')
    df1['label'] = etiquetas
    new_df1 = []
    for index, row in df1.iterrows():
        if row.Panther_id_x == row.Panther_id_y:
            continue
        else:
            new_df1.append([row.Panther_id_x, row.Entry, row.Panther_id_y, row.Panther, row.label])
    df1 = DataFrame(new_df1, columns =['GO_x','Entry','GO_y','Aspect','Term'])
    matrix = df1.pivot_table(values='Entry',index=['GO_x'],aggfunc=len,columns='GO_y')
    matrix2 = df1.pivot_table(values='Entry',index=['Term'],aggfunc=len,columns='GO_y')
    plt.subplots(figsize=(size_plot,size_plot))
    plt.gca().set_facecolor('whitesmoke')
    plt.xticks(size = xticks_size) # fontweight='bold'
    plt.yticks(size = yticks_size) # fontweight='bold'

    if label_x == 'GO':
        xlab = list(matrix.columns)
    if label_x == 'Term':
        xlab = list(matrix2.index)
    if label_y == 'GO':
        ylab = list(matrix.columns)
    if label_y == 'Term':
        ylab = list(matrix2.index)
        
    sns.heatmap(np.array(matrix),  square=True, annot = False, cmap = colors,
            #annot_kws={"size": 13, 'fontweight': 'bold'},fmt= '.0f',
            linewidths=0.005, linecolor='w',
            cbar=True, cbar_kws={"shrink": 0.7}, # 'label': '
            xticklabels= xlab,
            yticklabels= ylab)
    plt.gca().figure.axes[-1].set_ylabel('Interaction degrees\n(Proteins)', size= ylabel_size)


# In[ ]:


import warnings
warnings.filterwarnings("ignore")
heatmap_plot(table2, size_plot = 5)
plt.savefig('Heatmap_Pathways_python.png', dpi = 600, bbox_inches='tight')
print('Heatmap_Pathways_python.png')

# In[ ]:


# si desean incluir algún valor asociado a cada gen introducir una lista


# # con esta función se crean los dos archivos requeridos por el script de R de NeVOmics

# In[ ]:


def nodes_edges_R(df):
    graf = DataFrame(df.groupby(['Panther_id', 'label', 'pValue', 'fdr']).Entry.count()).reset_index()
    graf = graf.sort_values(by ='Entry',ascending=False).reset_index(drop=True)
    graf['num'] = list(graf.index.values + 1)
    gr =graf[['Panther_id','Entry','num','label', 'pValue', 'fdr']]
    gr.columns = ['Entry','Freq','num','Term', 'P', 'FDR']
    v1 = max(gr.num + 1)
    v2 = df[['Entry']].drop_duplicates().count()
    v3 = min(gr.Freq) * 0.5
    pant = df[['Entry']].drop_duplicates()
    pant['num'] = list(np.repeat(v1,v2))
    pant['Freq'] = list(np.repeat(v3,v2))
    nodes = pd.concat([gr,pant],sort=False)
    nodes = nodes.sort_values(by = 'P', ascending = True).reset_index(drop = True)
    nodes.to_csv('nodes_Panther_'+test_type+'.csv',index=None)
    print('\nFile:  nodes_Panther_'+test_type+'.csv')
    edges = df[['Panther_id','Entry']]
    edges.columns = ['GO','Entry']
    edges.to_csv('edges_Panther_'+test_type+'.csv',index=None)
    print('File:  edges_Panther_'+test_type+'.csv')
    #return nodes


# In[ ]:


import warnings
warnings.filterwarnings("ignore")
nodes_edges_R(table2)


# In[ ]:

import subprocess
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import os
import re
import urllib.request
from urllib.request import urlopen

print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nNow you can create the graphics using the NeVOmics R script\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
import tkinter as tk
from tkinter import filedialog
## Control of input file
print('R.exe location')
root = tk.Tk()
root.withdraw()
R_exe = filedialog.askopenfilename()
root.destroy()
print('=====> : ',R_exe)

print('\nFirst select the nodes file and then select the edges file')

file_R = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots_Modified.R', 'Plots_Modified.R')
#

run_uni = subprocess.Popen([R_exe, 'CMD', 'BATCH', 'Plots_Modified.R'])
run_uni.wait()


print('\nFinished process')

if os.path.exists(".RData"): os.remove(".RData")
if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




