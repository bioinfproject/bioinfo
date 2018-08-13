
# coding: utf-8

# In[ ]:


print('\n▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬ Loading data ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n..........\n')
## Packages import
from pandas import Series, DataFrame 
import pandas as pd
from pandas.compat import StringIO
import csv
import pandas
import pathlib
pd.set_option('max_rows',100000)
pd.set_option('max_colwidth',100000)
import urllib.request
import webbrowser
import re
import shutil, os
import numpy as np
from urllib.request import urlopen
from bs4 import BeautifulSoup
import requests
#from tqdm import tqdm
#tqdm.monitor_interval = 0
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")


# In[ ]:


## KEGG organisms list
kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()


# In[ ]:


## KEGG organisms list
#kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()

## Control of inputs= organism, prefix and file name
content_dir_list=os.listdir("./")
content_dir_str =" ".join(content_dir_list)

organism=input('▬▬▬▬▬▬▬▬\nStep 1: Submit a gender (e.g., Penicillium)\n▬▬▬▬▬▬▬▬\n-----> : ')
if re.search(organism,kegg_orgs):
    species=DataFrame(re.findall('T..........'+organism+' [a-z]{1,20}',
                                 kegg_orgs)).rename(columns={0:'Organism'}).replace({'\t':'   '},
                                                                                    regex=True) 
    print('\n',species)
    species1=DataFrame.to_string(species,header=None,index=None)
    inf=DataFrame(re.findall('T..........'+organism+' [a-z]{1,20}',
                                 kegg_orgs)).rename(columns={0:'Organism'}).replace({'\t':'   '},
                                                                                    regex=True)
    info=DataFrame.to_string(inf,header=False,index=False).split('\n')
    pref=DataFrame(info).replace({'T[0-9]{1,5}   ':'','   '+organism+'.*':''},regex=True)
    prefi=DataFrame.to_string(pref,header=None,index=None)
    
    
    Prefix=input('\n▬▬▬▬▬▬▬▬\nStep 2: Select Prefix (e.g., pcs)\n▬▬▬▬▬▬▬▬\n-----> : ')
else:
    print('\n!!!!!!! Organism not found !!!!!!!\n')
    organism=input('\n▬▬▬▬▬▬▬▬\nStep 1: Submit a gender (e.g., Penicillium)\n▬▬▬▬▬▬▬▬\n-----> : ')
    if re.search(organism,kegg_orgs):
        species=DataFrame(re.findall('T..........'+organism+' [a-z]{1,20}',
                                 kegg_orgs)).rename(columns={0:'Organism'}).replace({'\t':'   '},
                                                                                    regex=True) 
        print('\n',species)
        species1=DataFrame.to_string(species,header=None,index=None)
        inf=DataFrame(re.findall('T..........'+organism+' [a-z]{1,20}',
                                 kegg_orgs)).rename(columns={0:'Organism'}).replace({'\t':'   '},
                                                                                    regex=True)
        info=DataFrame.to_string(inf,header=False,index=False).split('\n')
        pref=DataFrame(info).replace({'T[0-9]{1,5}   ':'','   '+organism+'.*':''},regex=True)
        prefi=DataFrame.to_string(pref,header=None,index=None)
        Prefix=input('\n▬▬▬▬▬▬▬▬\nStep 2: Select Prefix (e.g., pcs)\n▬▬▬▬▬▬▬▬\n-----> : ')
    else:
        print('\n!!!!!!! Organism not found !!!!!!!\n')
        sys.exit()
        
if re.search(Prefix,prefi):
    input_file=input('\n▬▬▬▬▬▬▬▬\nStep 3: Submit gene list\n▬▬▬▬▬▬▬▬\n-----> : ')
else:
    print('\n!!!!!!! Prefix not found !!!!!!!\n')
    Prefix=input('\n▬▬▬▬▬▬▬▬\nStep 2: Select Prefix (e.g., pcs)\n▬▬▬▬▬▬▬▬\n-----> : ')
    if re.search(Prefix,prefi):
        input_file=input('\n▬▬▬▬▬▬▬▬\nStep 3: Submit gene list\n▬▬▬▬▬▬▬▬\n-----> : ')
    else:
        print('\n!!!!!!! Prefix not found !!!!!!!\n')
        sys.exit()
        
if re.search(input_file,content_dir_str):
    #print('')
    inp_file=pd.read_csv(input_file,sep='\t',header=None)
else:
    print('\n!!!!!!! File not found !!!!!!!\n')
    input_file=input('\n▬▬▬▬▬▬▬▬\nStep 3: Submit gene list\n▬▬▬▬▬▬▬▬\n-----> : ')
    if re.search(input_file,content_dir_str):
        #print('')
        inp_file=pd.read_csv(input_file,sep='\t',header=None)
    else:
        print('\n!!!!!!! File not found !!!!!!!\n')
        sys.exit()

## Create a folder
os.makedirs('data',exist_ok=True) 
## Download GeneMerge software
outputx=urllib.request.urlretrieve('https://raw.githubusercontent.com/eduardo1011/Programas/master/GeneMerge1.4.pl', 'GeneMerge1.4.pl')


## All KEGG pathways 
kegg_pathways = requests.get('http://rest.kegg.jp/list/pathway').content.decode()

## Pathways for mapping
Pathways_Mapping=(DataFrame(re.findall('[0-9]{5}.*',kegg_pathways)).replace({'^':Prefix},regex=True))[0].str.split('\t',
                                                                   expand=True).rename(columns={0:'Entry', 1:'Pathway'})
## Entry pathways species
path_species=requests.get('http://rest.kegg.jp/link/pathway/'+Prefix).content.decode()

## Data Frame of pathways species
Pathway_Species=(DataFrame(re.findall(Prefix+':.*',path_species)).replace({Prefix+':':'','path:':''},regex=True))[0].str.split('\t',
                                                                   expand=True).rename(columns={0:'Gene', 1:'Entry'})


## Mapping for obtain Specifit pathways
Specific_sp=pd.merge(Pathway_Species,Pathways_Mapping,on='Entry',how='left')

## Save pathways of species, without duplicates (Entry and Pathway), is used for Enrichment
Specific_sp[['Entry','Pathway']].drop_duplicates().to_csv('./data/Pathways_Species.txt',sep='\t',header=None,index=None)

concert_csv=DataFrame.to_csv(Specific_sp,index=None)
Specific_Path_Species=pd.read_csv(StringIO(concert_csv))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Open user's gene list 
## explore input file
if len(inp_file.columns) == 1:
    ## only gene list
    list_input=inp_file.rename(columns={0:'Gene'},index=str) 
else:
    print('')
if len(inp_file.columns) == 2:
    ## gene list and fold change
    list_input=inp_file.rename(columns={0:'Gene',1:'Fold Change'},index=str) ## gene list and fold change
else:
    print('')
if len(inp_file.columns) == 3:
    ## gene list, fold change and background
    list_input=inp_file.rename(columns={0:'Gene',1:'Fold Change',2:'Background'},index=str) 
else:
    print('')
##
## Preparation of background file 
## User's gackground column
if len(inp_file.columns) == 3:
    background=list_input[['Background']].rename(columns={'Background':'Gene'}) ## change column name
    back_with_pathw=pd.DataFrame.merge(background, Specific_Path_Species, how="left", on=list(list_input[['Gene']])).dropna()
    ## User's background list
    back_with_pathw[['Gene']].drop_duplicates().to_csv('./data/Background_Genes.txt',index=None,header=None)
    ## Gene list mapping against "back_with_pathw" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Gene']].dropna(),back_with_pathw, how="left", on=list(list_input[['Gene']])).dropna()
    list_input_match[['Gene']].drop_duplicates().to_csv('./data/Gene_list.txt',header=None,index=False, float_format='%.0f')
    v001="sed -i 's/[.]0//g' ./data/Gene_list.txt"
    subprocess.call(v001,shell=True)
    #output = Popen("sed -i 's/[.]0//g' ./data/Gene_list.txt", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    ## Build of kegg assotiation file
    aa=back_with_pathw.pivot_table(values='Entry',index=['Gene'],aggfunc=sum).reset_index()
    bb=aa[['Entry']].replace({Prefix:';'+Prefix},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Gene']],bb[['Entry']]],axis=1,join='outer').to_csv('./data/kegg_assotiation',sep='\t',header=None,index=None)
## File without background column
else:
    ## Background complete
    Pathway_Species[['Gene']].drop_duplicates().to_csv('./data/Background_Genes.txt',index=None,header=None)
    ## Mapping bwteen gene input and pathways
    list_input_match=pd.DataFrame.merge(list_input[['Gene']].dropna(), Specific_Path_Species, how="left", on=list(list_input[['Gene']])).dropna()
    ## Save Gene list found in KEGG database
    list_input_match[['Gene']].drop_duplicates().to_csv('./data/Gene_list.txt',header=None,index=False, float_format='%.0f')
    output = Popen("sed -i 's/[.]0//g' ./data/Gene_list.txt", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    ## Build of kegg assotiation file
    aa=Specific_Path_Species.pivot_table(values='Entry',index=['Gene'],aggfunc=sum).reset_index()
    bb=aa[['Entry']].replace({Prefix:';'+Prefix},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Gene']],bb[['Entry']]],axis=1,join='outer').to_csv('./data/kegg_assotiation',sep='\t',header=None,index=None)

## Over-representation test (hyper-geometric distriburion)
orden_columnas=[0,4,5,9,10,11]
output1 = Popen("perl GeneMerge1.4.pl ./data/kegg_assotiation ./data/Pathways_Species.txt ./data/Background_Genes.txt ./data/Gene_list.txt ./data/KEGG_Pathways_enrichment_GeneMerge.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/KEGG_Pathways_enrichment_GeneMerge.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

## Shows number of pathways found
kegg_01=pd.read_csv('./data/KEGG_Pathways_enrichment_GeneMerge.csv',usecols=orden_columnas)
min_val_kegg_adj_pval=kegg_01['P'].min()
max_val_kegg_adj_pval=kegg_01['P'].max()


#########>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if kegg_01[(kegg_01.P < 0.05)]['P'].count() >= 1:
    print('▬▬▬▬▬▬▬▬▬▬▬▬ ',kegg_01[(kegg_01.P < 0.05)]['P'].count(),'Pathways were found with P-value < 0.05')
    print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
    ## Correction of P-value    input_file=input('\nStep 3: Submit gene list\n\n-----> : ')
    ## Choice of method
    Bonferroni='Bonferroni'
    FDR='FDR'

    User_method=input('\n▬▬▬▬▬▬▬▬\nStep 4: Choose a Method (e.g., FDR / Bonferroni)\n▬▬▬▬▬▬▬▬\n-----> : ')
    if User_method == Bonferroni:
        a='b'
    else:
        if User_method == FDR:
            a='b'
        else:       
            #print('\n!!!!! Try again !!!!!')
            User_method=input('\n▬▬▬▬▬▬▬▬\nStep 4: Choose a Method (e.g., FDR / Bonferroni)\n▬▬▬▬▬▬▬▬\n-----> : ')
            if User_method == Bonferroni:
                a='b'
            else:
                if User_method == FDR:
                    a='b'
                else:
                    print('\n!!!!! Incorrect Method !!!!!')
                    sys.exit()

    ##  cut-off for Bonferroni                
    if User_method == Bonferroni:
        print('\n===== Must be in this range: ',min_val_kegg_adj_pval,' - ',max_val_kegg_adj_pval)
        Bon_value=input('\n▬▬▬▬▬▬▬▬\nStep 5: Choose a Value (e.g., 0.1)\n▬▬▬▬▬▬▬▬\n-----> : ')
        Bon_cut_off=float(Bon_value)
        file_name_value=Bon_value
        if min_val_kegg_adj_pval <= Bon_cut_off <= max_val_kegg_adj_pval:
            a='b'
        else:
            print('\nIncorrect value')
            print('\n===== Must be in this range: ',min_val_kegg_adj_pval,' - ',max_val_kegg_adj_pval)
            Bon_value=input('\n▬▬▬▬▬▬▬▬\nStep 5: Choose a Value (e.g., 0.1)\n▬▬▬▬▬▬▬▬\n-----> : ')
            Bon_cut_off=float(Bon_value)
            file_name_value=Bon_value
            if min_val_kegg_adj_pval <= Bon_cut_off <= max_val_kegg_adj_pval:
                a='b'
            else:
                print('\nIncorrect value')
                sys.exit()
            
    ## cut-off for FDR            
    else:
        if User_method == FDR:
            print('\n===== Must be in this range: ',min_val_kegg_adj_pval,' - ',max_val_kegg_adj_pval)
            FDR_value=input('\n▬▬▬▬▬▬▬▬\nStep 5: Choose a Value (e.g., 0.05)\n▬▬▬▬▬▬▬▬\n-----> : ')
            FDR_cut_off=float(FDR_value)*100
            file_name_value=FDR_value
            if min_val_kegg_adj_pval <= FDR_cut_off/100 <= max_val_kegg_adj_pval:
                output3 = Popen("perl GeneMerge1.4.pl ./data/kegg_assotiation ./data/Pathways_Species.txt ./data/Background_Genes.txt ./data/Gene_list.txt ./data/KEGG_Pathways_enrichment_GeneMerge.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                output4 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/KEGG_Pathways_enrichment_GeneMerge.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
            else:
                print('\nIncorrect value')
                print('\n===== Must be in this range: ',min_val_kegg_adj_pval,' - ',max_val_kegg_adj_pval)
                FDR_value=input('\n▬▬▬▬▬▬▬▬\nStep 5: Choose a Value (e.g., 0.05)\n▬▬▬▬▬▬▬▬\n-----> : ')
                FDR_cut_off=float(FDR_value)*100
                file_name_value=FDR_value
                if min_val_kegg_adj_pval <= FDR_cut_off/100 <= max_val_kegg_adj_pval:
                    output3 = Popen("perl GeneMerge1.4.pl ./data/kegg_assotiation ./data/Pathways_Species.txt ./data/Background_Genes.txt ./data/Gene_list.txt ./data/KEGG_Pathways_enrichment_GeneMerge.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    output4 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/KEGG_Pathways_enrichment_GeneMerge.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                else:
                    print('\nIncorrect value')
                    sys.exit()
else:
    print('No significant Pathways were found')
    sys.exit()
                    
## Open raw data of enrichment
kegg=pd.read_csv('./data/KEGG_Pathways_enrichment_GeneMerge.csv',usecols=orden_columnas)

if User_method == Bonferroni:
    kegg_user_cut_off=kegg[(kegg.adj_pval <= Bon_cut_off) == True]
else:
    a='b'
if User_method == FDR:
    filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
    kegg_user_cut_off=kegg[kegg[filter_fdr_T].str.contains("T") == True]
else: 
    a='b'

## Enrichment pathways after cut-off
kegg_user_cut_off['Freq']=kegg_user_cut_off['Entry'].str.split().str.len()

## Data editing
enriched_genes=DataFrame(kegg_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
equalize=str(kegg_user_cut_off.GO*kegg_user_cut_off.Freq)
id_pathways=re.findall(Prefix+'.....',equalize)
enriched_pathways=DataFrame(id_pathways)
frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Gene",1:'Entry'}, index=str)#.reset_index(drop=True) 
Gene_Entry_Pathway_Enriched=gene_kegg.merge(Pathways_Mapping,on='Entry',how='left')

## File name for edges
edges_file_name='Edges_KEGG_'+organism+'_'+User_method+'_'+file_name_value+'.csv'
Gene_Entry_Pathway_Enriched[['Entry','Gene']].to_csv(edges_file_name,index=None)
Gene_Entry_Pathway_Enriched.to_csv('Enrichment_Analysis_'+organism+'_'+User_method+'_'+file_name_value+'.xlsx',index=None)

## Preparing data frame for nodes file
frame=Gene_Entry_Pathway_Enriched[['Entry']].drop_duplicates().reset_index(drop=True)
frame['Freq']=kegg_user_cut_off['Freq']
frame['num']=frame.reset_index().index + 1
df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Gene']],index=None)
df2=pd.read_csv(StringIO(df1))
#
if len(inp_file.columns) >= 2:
    entry_exp=df2.merge(list_input[['Gene','Fold Change']],on='Gene',
                        how='left').rename(columns={'Gene':'Entry','Fold Change':'Exp'},
                                           index=str).drop_duplicates().reset_index(drop=True)
    exp = float(entry_exp[['Exp']].count())
    if  exp > 0:
        print('')
    else:
        entry_exp=entry_exp[['Entry']]
else:
    entry_exp=df2.rename(columns={'Gene':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
#
entry_exp['Freq']=float(kegg_user_cut_off['Freq'].min()*.5)
entry_exp['num']=float(frame['num'].count() + 1)
nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
entry_pathways=Gene_Entry_Pathway_Enriched[['Entry','Pathway']].drop_duplicates()
nodes=nodes.merge(entry_pathways,on='Entry',how='left').rename(columns={'Pathway':'Term'},index=str).reset_index(drop=True)
nodes['P']=kegg_user_cut_off['P']
nodes_file_name='Nodes_KEGG_'+organism+'_'+User_method+'_'+file_name_value+'.csv'
nodes.to_csv(nodes_file_name,index=None)
#
##
max_val_pathways=frame['Entry'].count()
url_webbrowser_path='http://www.kegg.jp/kegg-bin/show_pathway?'
output_google= Popen("export BROWSER=google-chrome",shell=True, stdin=PIPE,stdout=PIPE,stderr=STDOUT,close_fds=True).stdout.read()
output_firefox= Popen("export BROWSER=firefox",shell=True, stdin=PIPE,stdout=PIPE,stderr=STDOUT,close_fds=True).stdout.read()
#
import webbrowser as wb
import webbrowser
##
info_str="\n".join(info)
kegg_genome_id=DataFrame(re.findall('T[0-9]{1,5}   '+Prefix,info_str)).replace({'   '+Prefix+'.*':''},regex=True)[0].iloc[0]
output = Popen("python -m webbrowser https://www.kegg.jp/dbget-bin/www_bget?gn:"+kegg_genome_id+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
#webbrowser.open_new_tab("https://www.kegg.jp/dbget-bin/www_bget?gn:"+kegg_genome_id+"")   
##
print('\n▬▬▬▬▬▬▬▬\nStep 6: Open pathways in the browser\n▬▬▬▬▬▬▬▬\n..........\n')
#
if max_val_pathways > 0:
    aa=kegg_user_cut_off.iloc[[0]]
    first_path=aa['GO'].iloc[0]
    bb=aa[['Entry']].replace({', /':'','^':url_webbrowser_path+first_path,'/$':''},regex=True)
    first=bb['Entry'].iloc[0]
    output = Popen("python -m webbrowser "+first+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #  organism+' '++'\n\n'++'\n'
    brow=[organism+' '+first_path+'\n\n'+first+'\n']
else:
    a='b'
if max_val_pathways > 1:
    cc=kegg_user_cut_off.iloc[[1]]
    second_path=cc['GO'].iloc[0]
    dd=cc[['Entry']].replace({', /':'','^':url_webbrowser_path+second_path,'/$':''},regex=True)
    second=dd['Entry'].iloc[0]
    output = Popen("python -m webbrowser "+second+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n']
else:
    a='b'
if max_val_pathways > 2:
    ee=kegg_user_cut_off.iloc[[2]]
    third_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+third_path,'/$':''},regex=True)
    third=ff['Entry'].iloc[0]
    output = Popen("python -m webbrowser "+third+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n']
else:
    a='b'
if max_val_pathways > 3:
    ee=kegg_user_cut_off.iloc[[3]]
    fourth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+fourth_path,'/$':''},regex=True)
    fourth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+fourth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n']
else:
    a='b'
if max_val_pathways > 4:
    ee=kegg_user_cut_off.iloc[[4]]
    fifth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+fifth_path,'/$':''},regex=True)
    fifth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+fifth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n']
else:
    a='b'
if max_val_pathways > 5:
    ee=kegg_user_cut_off.iloc[[5]]
    sixth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+sixth_path,'/$':''},regex=True)
    sixth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+sixth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n']
else:
    a='b'
if max_val_pathways > 6:
    ee=kegg_user_cut_off.iloc[[6]]
    seventh_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+seventh_path,'/$':''},regex=True)
    seventh=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+seventh+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n']
else:
    a='b'
if max_val_pathways > 7:
    ee=kegg_user_cut_off.iloc[[7]]
    eighth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+eighth_path,'/$':''},regex=True)
    eighth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+eighth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n']
else:
    a='b'
if max_val_pathways > 8:
    ee=kegg_user_cut_off.iloc[[8]]
    ninth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+ninth_path,'/$':''},regex=True)
    ninth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+ninth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n']
else:
    a='b'
if max_val_pathways > 9:
    ee=kegg_user_cut_off.iloc[[9]]
    tenth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+tenth_path,'/$':''},regex=True)
    tenth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+tenth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n',
          organism+' '+tenth_path+'\n\n'+tenth+'\n']
else:
    a='b'
###
if max_val_pathways > 10:
    ee=kegg_user_cut_off.iloc[[10]]
    Eleventh_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+Eleventh_path,'/$':''},regex=True)
    Eleventh=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+Eleventh+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n',
          organism+' '+tenth_path+'\n\n'+tenth+'\n',
          organism+' '+Eleventh_path+'\n\n'+Eleventh+'\n']
else:
    a='b'
if max_val_pathways > 11:
    ee=kegg_user_cut_off.iloc[[11]]
    Twelfth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+Twelfth_path,'/$':''},regex=True)
    Twelfth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+Twelfth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first,
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n',
          organism+' '+tenth_path+'\n\n'+tenth+'\n',
          organism+' '+Eleventh_path+'\n\n'+Eleventh+'\n',
          organism+' '+Twelfth_path+'\n\n'+Twelfth+'\n']
else:
    a='b'
if max_val_pathways > 12:
    ee=kegg_user_cut_off.iloc[[12]]
    Thirteenth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+Thirteenth_path,'/$':''},regex=True)
    Thirteenth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+Thirteenth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n',
          organism+' '+tenth_path+'\n\n'+tenth+'\n',
          organism+' '+Eleventh_path+'\n\n'+Eleventh+'\n',
          organism+' '+Twelfth_path+'\n\n'+Twelfth+'\n',
          organism+' '+Thirteenth_path+'\n\n'+Thirteenth+'\n']
else:
    a='b'
if max_val_pathways > 13:
    ee=kegg_user_cut_off.iloc[[13]]
    Fourteenth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+Fourteenth_path,'/$':''},regex=True)
    Fourteenth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+Fourteenth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n',
          organism+' '+tenth_path+'\n\n'+tenth+'\n',
          organism+' '+Eleventh_path+'\n\n'+Eleventh+'\n',
          organism+' '+Twelfth_path+'\n\n'+Twelfth+'\n',
          organism+' '+Thirteenth_path+'\n\n'+Thirteenth+'\n',
          organism+' '+Fourteenth_path+'\n\n'+Fourteenth+'\n']
else:
    a='b'
if max_val_pathways > 14:
    ee=kegg_user_cut_off.iloc[[14]]
    Fifteenth_path=ee['GO'].iloc[0]
    ff=ee[['Entry']].replace({', /':'','^':url_webbrowser_path+Fifteenth_path,'/$':''},regex=True)
    Fifteenth=ff['Entry'].iloc[0]
    #output = Popen("python -m webbrowser "+Fifteenth+"", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    brow=[organism+' '+first_path+'\n\n'+first+'\n',
          organism+' '+second_path+'\n\n'+second+'\n',
          organism+' '+third_path+'\n\n'+third+'\n',
          organism+' '+fourth_path+'\n\n'+fourth+'\n',
          organism+' '+fifth_path+'\n\n'+fifth+'\n',
          organism+' '+sixth_path+'\n\n'+sixth+'\n',
          organism+' '+seventh_path+'\n\n'+seventh+'\n',
          organism+' '+eighth_path+'\n\n'+eighth+'\n',
          organism+' '+ninth_path+'\n\n'+ninth+'\n',
          organism+' '+tenth_path+'\n\n'+tenth+'\n',
          organism+' '+Eleventh_path+'\n\n'+Eleventh+'\n',
          organism+' '+Twelfth_path+'\n\n'+Twelfth+'\n',
          organism+' '+Thirteenth_path+'\n\n'+Thirteenth+'\n',
          organism+' '+Fourteenth_path+'\n\n'+Fourteenth+'\n',
          organism+' '+Fifteenth_path+'\n\n'+Fifteenth+'\n']
else:
    a='b'


# In[ ]:

## Control of creating of directories
content_dir_list=os.listdir("./")
content_dir_str =" ".join(content_dir_list)
if re.search('job_[0-9]{1,3}_KEGG',content_dir_str):
    plot_dir=DataFrame(re.findall('job_[0-9]{1,3}_KEGG',content_dir_str)).sort_values(by=[0],ascending=False).reset_index(drop=True)
    plot_dir2=plot_dir[0].str.split('_', expand=True)
    plot_dir3=float(plot_dir2[1].iloc[0])
    plot_dir4=str(int(float(plot_dir3+1)))
    if float(plot_dir4) >= 10:
        os.makedirs('job_'+plot_dir4+'_KEGG/'+'job_'+plot_dir4+'_KEGG_plots')
        dir_name_plots='job_'+plot_dir4+'_KEGG'
    else:
        if os.path.exists('job_0'+plot_dir4+'_KEGG'):
            dir_name_plots='job_0'+plot_dir4+'_KEGG'
        else:
            os.makedirs('job_0'+plot_dir4+'_KEGG/'+'job_0'+plot_dir4+'_KEGG_plots')
            dir_name_plots='job_0'+plot_dir4+'_KEGG'
else:
    os.makedirs('job_01_KEGG/job_01_KEGG_plots')
    dir_name_plots='job_01_KEGG'
dir_plots='./'+dir_name_plots+'/'+dir_name_plots+'_plots/'


# In[ ]:


## Open R script from github
r_script=requests.get("https://raw.githubusercontent.com/eduardo1011/Programas/master/Enrichment_Plots.R").content.decode()
R_script_enrich = re.sub("./plots/",dir_plots,r_script)
## Create file with R script
f= open("Enrichment_Plots.R","w")
f.write('#\n#\n# Libraries\n#\n'+
        'library(tidyverse)\n'+
        'library(tidygraph)\n'+
        'library(viridis)\n'+
        'library(ggraph)\n'+
        'library(circlize)\n'+
        'library(RColorBrewer)\n'+
        'library(igraph)\n'+
        'library(cowplot)\n'+
        'library(grid)\n'+
        'library(networkD3)\n'+
        'library(UpSetR)\n'+
        '#\n#\nnodes = read_csv("'+
        nodes_file_name+'")\n#\nlinks = read_csv("'+
        edges_file_name+'")\n#\n#\n'+R_script_enrich)
f.close()


# In[ ]:


print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics\n▬▬▬▬▬▬▬▬\n..........\n')


# In[ ]:

import subprocess
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
process.wait()


# In[ ]:


if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
if os.path.exists("GeneMerge1.4.pl"): os.remove("GeneMerge1.4.pl")
if os.path.exists('data'): shutil.rmtree('data')
if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")


# In[ ]:


## Analysis report
aaa=' '+str(kegg[(kegg.P < 0.05)]['P'].count())+' Pathways were found with P-value < 0.05'
bbb=DataFrame.to_string(kegg[(kegg.P < 0.05)].rename(columns={'GO':'Pathway'},index=str),index=True)
content_job_plots=os.listdir(dir_plots)
content_job_str ="\n".join(content_job_plots)
links="\n".join(brow)
#
f= open("Report_"+dir_name_plots+".txt","w")
f.write('#\n# ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n#'+
        ' Enrichment Analysis Report\n#\n# Program Name: \n#\n#\n#\n#\n#\n# ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬'+
        '\n# Parameters\n#'
        '\n# Organism:           '+organism+
        '\n# KEGG GENOME:        '+kegg_genome_id+
        '\n# KEGG organism code: '+Prefix+
        '\n# Study File:         '+input_file+
        '\n# Correction Method:  '+User_method+
        '\n# Cut-off:            '+file_name_value+
        '\n#\n#\n# ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n#'+
        aaa+'\n#\n'+
        bbb+'#\n\n# ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬'+
        '\n# Links to the KEGG pathways maps\n#\n'+
        links+
        '\n# ▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬'+
        '\n# Graphics stored in '+dir_plots+' directory\n#\n'+
        content_job_str+'\n#\n#\n#')
f.close()


# In[ ]:


report="Report_"+dir_name_plots+".txt"
xlsx='Enrichment_Analysis_'+organism+'_'+User_method+'_'+file_name_value+'.xlsx'
if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, dir_name_plots)
if os.path.exists(edges_file_name): shutil.move(edges_file_name, dir_name_plots)
#if os.path.exists('Enrichment_Plots.R'): shutil.move('Enrichment_Plots.R', dir_name_plots)
if os.path.exists(report): shutil.move(report, dir_name_plots)
if os.path.exists(xlsx): shutil.move(xlsx, dir_name_plots)
if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
output3 = Popen("exit", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
if os.path.exists("enrichment.py"): os.remove("enrichment.py")
