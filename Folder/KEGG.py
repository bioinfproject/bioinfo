
# coding: utf-8

# In[ ]:


## Modules import
import datetime
start = datetime.datetime.now()
from pandas import Series, DataFrame 
import pandas as pd
from pandas.compat import StringIO
import csv
import pandas
import pathlib
#pd.set_option('max_rows',100000)
#pd.set_option('max_colwidth',100000)
import urllib.request
import webbrowser
import re
import shutil, os
import numpy as np
from urllib.request import urlopen
#from bs4 import BeautifulSoup
import requests
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import sys
import warnings
#warnings.filterwarnings("ignore")
from datetime import datetime 
inicio_total = datetime.now()
import os, fnmatch



# #### uniprot_peni_list
# #### Cancer_analysis_platelet.txt

# In[ ]:

import tkinter as tk
from tkinter import filedialog
## Control of input file
print('\n[ Step 1: Submit file (Uniprot IDs) ]')
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
root.destroy()
if file_path == ():
    print('\n!!!!!!! File not found !!!!!!!')
    import tkinter as tk
    from tkinter import filedialog
    ## Control of input file
    print('\n[ Step 1: Submit file (Uniprot IDs) ]')
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

## read file
inp_file=pd.read_csv(file_path,sep='\t',header=None)   
    
## explore input file
if len(inp_file.columns) == 1:
    ## only gene list
    list_input=inp_file.rename(columns={0:'Entry'},index=str) 
else:
    a='b'
if len(inp_file.columns) == 2:
    ## gene list and fold change
    list_input=inp_file.rename(columns={0:'Entry',1:'Fold Change'},index=str) ## gene list and fold change
else:
    a='b'
if len(inp_file.columns) == 3:
    ## gene list, fold change and background
    list_input=inp_file.rename(columns={0:'Entry',1:'Fold Change',2:'Background'},index=str) 
else:
    a='b'

## info of protein
for i in list(list_input['Entry'].head()):
    a=requests.get("http://rest.kegg.jp/conv/genes/uniprot:"+str(i)).content.decode()
    if a != '\n':
        if a == '':
            print('\n!!!!!!! Organism not found, check your identifiers list !!!!!!!')
            sys.exit()
        else:
            pref=re.findall('[a-z]{3,4}',re.findall('\t[a-z]{3,4}:',a)[0])[0]
            break
if a == '\n':
    print('\n!!!!!!! Organism not found, check your identifiers list !!!!!!!')
    sys.exit()
else:
    b=requests.get('http://rest.kegg.jp/info/'+pref+'').content.decode()
    organism=re.sub(' KEGG Genes Database','',re.sub('\n.*','',re.sub('T[0-9]{5}           ','',b)))
    print('\n-- ','Organism identified: ',organism)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Control of Uniprot and GOA directories
ls=[]
for i in os.listdir("./"):
    if re.search('job_KEGG_[0-9]{1,3}',str(i)):
        ls.append(i)
if ls == []:
    new_folder='job_KEGG_1'
    xoxo='1'
    os.makedirs('job_KEGG_1/job_KEGG_1/job_KEGG_plots_1',exist_ok=True)
    level_1_kegg=new_folder+'/job_KEGG_'+xoxo+'/'
    level_2_kegg=new_folder+'/job_KEGG_'+xoxo+'/job_KEGG_plots_'+xoxo+'/'
else:
    n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
    xoxo=str(n+1)
    n=str(n)
    old_folder=''.join(re.findall('job_KEGG_'+n,str(ls)))
    new_folder=re.sub(n,xoxo,old_folder)
    os.makedirs(new_folder+'/job_KEGG_'+xoxo+'/job_KEGG_plots_'+xoxo+'',exist_ok=True)
    level_1_kegg=new_folder+'/job_KEGG_'+xoxo+'/'
    level_2_kegg=new_folder+'/job_KEGG_'+xoxo+'/job_KEGG_plots_'+xoxo+'/'
    
## Create a folder
os.makedirs('data',exist_ok=True)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
######### DataBase
Pathways_database=requests.get('http://rest.kegg.jp/list/pathway').content.decode()
Pathways_database=pd.read_csv(StringIO(Pathways_database),sep='\t',header=None,names=['GO','Term']).replace({'^path:map':pref},regex=True)
Pathways_database.to_csv('data/Pathways.txt',sep='\t',index=None)
#########

## include information to kegg from Uniprot

ses = requests.Session()
ses.max_redirects = 100
first_entry=inp_file[0].iloc[0]
id_organism = ses.get("https://www.uniprot.org/uniprot/?query="+first_entry+"&sort=score&columns=organism-id,organism&format=tab&limit=1").content.decode()
Prefix=DataFrame(re.findall('[0-9]{1,30}',id_organism))[0].iloc[0]

if Prefix == '9606':
    inf1=requests.get('https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:'+Prefix+'&format=tab&columns=id,database(DisGeNET)').content.decode()
    inf1=pd.read_csv(StringIO(inf1),sep='\t').rename(columns={'Cross-reference (DisGeNET)':'Entry_Kegg'}).replace({';.*':''},regex=True)
    inf1=inf1[['Entry_Kegg','Entry']]
    inf2=requests.get('https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:'+Prefix+'&format=tab&columns=id,database(GeneID)').content.decode()
    inf2=pd.read_csv(StringIO(inf2),sep='\t').rename(columns={'Cross-reference (GeneID)':'Entry_Kegg'}).replace({';.*':''},regex=True).dropna()
    inf2=inf2[['Entry_Kegg','Entry']]
    frames=[inf1,inf2]
    Kegg_Uniprot=pd.concat(frames,axis=0,sort=False).dropna().reset_index(drop=True).drop_duplicates()
    # all kegg-id and pathway-id
    dd=requests.get('http://rest.kegg.jp/link/pathway/'+pref+'').content.decode()
    kegg_path_ID=pd.read_csv(StringIO(dd),sep='\t',header=None,names=['Entry_Kegg','GO']).replace({'^'+pref+':|path:':''},regex=True)
    # all kegg-id and pathway-description
    ee=requests.get('http://rest.kegg.jp/list/pathway/'+pref+'').content.decode()
    kegg_pathways=pd.read_csv(StringIO(ee),sep='\t',header=None,names=['GO','Description']).replace({'path:|- '+organism[0:5]+'.*':''},regex=True)
else:
    # all kegg-id and pathway-id
    dd=requests.get('http://rest.kegg.jp/link/pathway/'+pref+'').content.decode()
    kegg_path_ID=pd.read_csv(StringIO(dd),sep='\t',header=None,names=['Entry_Kegg','GO']).replace({'^'+pref+':|path:':''},regex=True)
    # all kegg-id and pathway-description
    ee=requests.get('http://rest.kegg.jp/list/pathway/'+pref+'').content.decode()
    kegg_pathways=pd.read_csv(StringIO(ee),sep='\t',header=None,names=['GO','Description']).replace({'path:|- '+organism[0:5]+'.*':''},regex=True)
    # all kegg-id and uniprot
    ff=requests.get('http://rest.kegg.jp/conv/uniprot/'+pref+'').content.decode()
    Kegg_Uniprot=pd.read_csv(StringIO(ff),sep='\t',header=None,names=['Entry_Kegg','Entry']).replace({'^'+pref+':|up:':''},regex=True)

# In[ ]:


if len(inp_file.columns) == 3:
    #### With background
    # 1.- Preparation of background
    background=list_input[['Background']].rename(columns={'Background':'Entry'})
    background_info=pd.merge(background,Kegg_Uniprot.dropna(),on='Entry',how='left').dropna().drop_duplicates()
    background_info=pd.merge(background_info,kegg_path_ID.dropna(),on='Entry_Kegg',how='left').dropna().drop_duplicates()
    background_info[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/Background.txt',index=None)
    
    # 2.- Preparation of list with pathways
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),background_info,how="left", on='Entry').dropna()
    list_input_match[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/List.txt',index=None)
    
    # 3.- background with: Entry	GO, for association file
    background_info[['Entry_Kegg','GO']].rename(columns={'Entry_Kegg':'Entry'}).to_csv('data/Association.txt',index=None,sep='\t')
    
    #exploratory analysis of P-value in data
    file='Pathways.txt'
    fdr_min=0.05
    subprocess.call(["python","Hypergeometric_distribution.py",file,str(fdr_min),'data/Enrichment_analysis_Path.tsv'])
    enrich_P=pd.read_csv('data/Enrichment_analysis_Path.tsv',sep='\t')
        
    if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
        print('\n----------------------------------------------------\n'+
                '                KEGG Pathways Annotation\n----------------------------------------------------'+
                '\n\n*****  KEGG Pathways:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
                'Pathways were found with P-value < 0.05')
            
        ## Choice of correction corection method
        methods = ['Bonferroni','FDR']
        User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        method_P=[]
        for x in methods:
            if x == User_method:
                method_P.append(x)
        if method_P == []:
            print('\n!!!!! Incorrect Method !!!!!')
            User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
            for x in methods:
                if x == User_method:
                    method_P.append(x)
        if method_P == []:
            print('\n!!!!! Incorrect Method !!!!!')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('data'): shutil.rmtree('data')
            sys.exit()
        ########################################
        ########################################
        ## loop control User value
        a=input('[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
        if a == '':
            a='aaa'
        if re.search(r'[A-Za-z]{1,10}',a):
            print('Incorrect value')
            a=input('[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
            if a == '':
                a='aaa'
            if re.search(r'[A-Za-z]{1,10}',a):
                print('Incorrect value')
                if os.path.exists(new_folder): shutil.rmtree(new_folder)
                if os.path.exists('data'): shutil.rmtree('data')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_P=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_P=float(a)
            else:
                print('Incorrect value')
                a=input('[ Must be in this range: 0.005 - 0.5 ]\n=====> : ')
                if a == '':
                    a='aaa'
                if re.search(r'[A-Za-z]{1,10}',a):
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_P=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        sys.exit()
        ########################################
        ########################################
        analysis='Pathways.txt'
        if ''.join(method_P) == 'Bonferroni':
            subprocess.call(["python","Hypergeometric_distribution.py",
                                analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
            enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
            results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
            proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_P=results_process_P[['GO']].count()[0]

        else:
            if ''.join(method_P) == 'FDR':
                subprocess.call(["python","Hypergeometric_distribution.py",
                                    analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
                enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
                results_process_P=enrich_P[enrich_P.Sig == 'T']
                proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_P=results_process_P[['GO']].count()[0]
        
         # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_P['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),background_info[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            lis['entry']='entry'
            lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
            lis=lis.groupby('entry')['ent'].sum().reset_index()
            www=pd.merge(list_input_match[['Entry_Kegg']],background_info,on='Entry_Kegg',how='left')
            uuu=www.groupby('GO')['Entry'].count().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True).drop_duplicates()
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+file_path+
                      '\nAssociation file name\t'+analysis+
                      '\nTotal number of background\t'+str(background['Entry'].drop_duplicates().count())+
                      '\nTotal number of list\t'+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with Pathways\t'+str(background_info['Entry'].drop_duplicates().count())+
                      '\nList input with Pathways\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+
                      '\nCorrection Method\t'+''.join(method_P)+
                      '\nValue\t'+str(User_value_P)+
                      '\n\t\n'+
                      '\nProteins with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]
            
            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})

            a01=[]
            for index, row in results_process_P.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            process_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_KEGG_Enrichment_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            process_goa.to_csv(level_1_kegg+edges_file_name,index=None)

                
            enrich_P=enrich_P.rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
            process_goa=process_goa.rename(columns={'GO':'Path'})
            ##
            process_goa=process_goa.merge(list_input_match[['Entry','Entry_Kegg']].rename(columns={'Entry':'UniProt','Entry_Kegg':'Entry'}),on='Entry',how='left')
            df2=process_goa[['Entry']]
            process_goa=pd.merge(process_goa,kegg_pathways[['GO','Description']].rename(columns={'GO':'Path'}),on='Path',how='left').drop_duplicates()
            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Pathways_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant KEGG Pathways',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_goa.to_excel(writer,'Edges Pathways',index=False)
            writer.save()

            frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            #df2=process_goa[['Entry']] ## fue definido antes de convertir la columna entry a formato numérico

            #
            if len(inp_file.columns) >= 2:
                list_in=list_input_match.merge(list_input[['Entry','Fold Change']],on='Entry',
                                                  how='left').rename(columns={'Fold Change':'Exp'},
                                                                     index=str).drop_duplicates().reset_index(drop=True)
                
                entry_exp=df2.merge(list_in[['Entry_Kegg','Exp']].rename(columns={'Entry_Kegg':'Entry'}),
                                    on='Entry',how='left')
                exp = float(entry_exp[['Exp']].count())
                if  exp > 0:
                    a='b'
                else:
                    entry_exp=entry_exp[['Entry']]
            else:
                entry_exp=df2

            entry_exp['Freq']=float(frame['Freq'].min()*.5)
            entry_exp['num']=float(frame['num'].count() + 1)
            nodes=pd.concat([frame,entry_exp],sort=False)
            nodes_file_name='nodes_KEGG_Enrichment_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_kegg+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/eduardo1011/Programas/master/Enrichment_Plots.R').content.decode()
            R_script_enrich = re.sub('qwertyuiop',level_1_kegg+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_kegg+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_kegg,R_script_enrich) # store plots
            R_script_enrich = re.sub('ASPECT','Pathways',R_script_enrich)
            R_script_enrich = re.sub('poiuytrewq','Pathways Annotation (BP)',R_script_enrich) # legend
            f= open(level_1_kegg+'/Kegg_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\n*****  KEGG Pathways: No significant terms were found')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),background_info[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            lis['entry']='entry'
            lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
            lis=lis.groupby('entry')['ent'].sum().reset_index()
            www=pd.merge(list_input_match[['Entry_Kegg']],background_info,on='Entry_Kegg',how='left')
            uuu=www.groupby('GO')['Entry'].count().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True).drop_duplicates()
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+file_path+
                      '\nAssociation file name\t'+analysis+
                      '\nTotal number of background\t'+str(background['Entry'].drop_duplicates().count())+
                      '\nTotal number of list\t'+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with Pathways\t'+str(background_info['Entry'].drop_duplicates().count())+
                      '\nList input with Pathways\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+  # ****
                      '\nCorrection Method\t'+''.join(method_P)+ # ***
                      '\nValue\t'+str(User_value_P)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_P, information], axis=0, sort=False) # ***
            process_goa=results_process_P[['GO','entry']].rename(columns={'entry':'Entry'})

            enrich_P=enrich_P.rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
            process_goa=process_goa.rename(columns={'GO':'Path'})
            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Pathways_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant KEGG Pathways',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            #process_goa.to_excel(writer,'Edges Pathways',index=False)
            writer.save()
    
    else:
        print('\n*****  KEGG Pathways: No significant terms were found')
        results_process_P=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_P=results_process_P[['entry']].count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        method_P='0'
        User_value_P='0'
            
## Without background           
else:
    print('\n* [ The analysis will be done using the whole information ]')
    #### With background
    # 1.- Preparation of background
    background_info=kegg_path_ID.drop_duplicates()
    background_info=pd.merge(background_info,Kegg_Uniprot.dropna(),on='Entry_Kegg',how='left').dropna().drop_duplicates()
    background_info[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/Background.txt',index=None)
    # 2.- Preparation of list with pathways
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),background_info,how="left", on='Entry').dropna()
    list_input_match[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/List.txt',index=None)
    # 3.- background with: Entry	GO, for association file
    background_info[['Entry_Kegg','GO']].rename(columns={'Entry_Kegg':'Entry'}).to_csv('data/Association.txt',index=None,sep='\t')
     
    #exploratory analysis of P-value in data
    file='Pathways.txt'
    fdr_min=0.05
    subprocess.call(["python","Hypergeometric_distribution.py",file,str(fdr_min),'data/Enrichment_analysis_Path.tsv'])
    enrich_P=pd.read_csv('data/Enrichment_analysis_Path.tsv',sep='\t')
        
    if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
        print('\n----------------------------------------------------\n'+
                '                KEGG Pathways Annotation\n----------------------------------------------------'+
                '\n\n*****  KEGG Pathways:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
                'Pathways were found with P-value < 0.05')
            
        ## Choice of correction corection method
        methods = ['Bonferroni','FDR']
        User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        method_P=[]
        for x in methods:
            if x == User_method:
                method_P.append(x)
        if method_P == []:
            print('\n!!!!! Incorrect Method !!!!!')
            User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
            for x in methods:
                if x == User_method:
                    method_P.append(x)
        if method_P == []:
            print('\n!!!!! Incorrect Method !!!!!')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('data'): shutil.rmtree('data')
            sys.exit()
        ########################################
        ########################################
        ## loop control User value
        a=input('[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
        if a == '':
            a='aaa'
        if re.search(r'[A-Za-z]{1,10}',a):
            print('Incorrect value')
            a=input('[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
            if a == '':
                a='aaa'
            if re.search(r'[A-Za-z]{1,10}',a):
                print('Incorrect value')
                if os.path.exists(new_folder): shutil.rmtree(new_folder)
                if os.path.exists('data'): shutil.rmtree('data')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_P=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_P=float(a)
            else:
                print('Incorrect value')
                a=input('[ Must be in this range: 0.005 - 0.5 ]\n=====> : ')
                if a == '':
                    a='aaa'
                if re.search(r'[A-Za-z]{1,10}',a):
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_P=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        sys.exit()
        ########################################
        ########################################
        analysis='Pathways.txt'
        if ''.join(method_P) == 'Bonferroni':
            subprocess.call(["python","Hypergeometric_distribution.py",
                                analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
            enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
            results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
            proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_P=results_process_P[['GO']].count()[0]

        else:
            if ''.join(method_P) == 'FDR':
                subprocess.call(["python","Hypergeometric_distribution.py",
                                    analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
                enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
                results_process_P=enrich_P[enrich_P.Sig == 'T']
                proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_P=results_process_P[['GO']].count()[0]
        
         # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_P['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),background_info[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            lis['entry']='entry'
            lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
            lis=lis.groupby('entry')['ent'].sum().reset_index()
            www=pd.merge(list_input_match[['Entry_Kegg']],background_info,on='Entry_Kegg',how='left')
            uuu=www.groupby('GO')['Entry'].count().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True).drop_duplicates()
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+file_path+
                      '\nAssociation file name\t'+analysis+
                      '\nTotal number of background\t'+str(background_info['Entry'].drop_duplicates().count())+
                      '\nTotal number of list\t'+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with Pathways\t'+str(background_info['Entry'].drop_duplicates().count())+
                      '\nList input with Pathways\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+
                      '\nCorrection Method\t'+''.join(method_P)+
                      '\nValue\t'+str(User_value_P)+
                      '\n\t\n'+
                      '\nProteins with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]
            
            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})

            a01=[]
            for index, row in results_process_P.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            process_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_KEGG_Enrichment_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            process_goa.to_csv(level_1_kegg+edges_file_name,index=None)

                
            enrich_P=enrich_P.rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
            process_goa=process_goa.rename(columns={'GO':'Path'})
            ##
            process_goa=process_goa.merge(list_input_match[['Entry','Entry_Kegg']].rename(columns={'Entry':'UniProt','Entry_Kegg':'Entry'}),on='Entry',how='left')
            df2=process_goa[['Entry']]
            process_goa=pd.merge(process_goa,kegg_pathways[['GO','Description']].rename(columns={'GO':'Path'}),on='Path',how='left').drop_duplicates()
            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Pathways_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant KEGG Pathways',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_goa.to_excel(writer,'Edges Pathways',index=False)
            writer.save()

            frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            #df2=process_goa[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                list_in=list_input_match.merge(list_input[['Entry','Fold Change']],on='Entry',
                                               how='left').rename(columns={'Fold Change':'Exp'},
                                                                  index=str).drop_duplicates().reset_index(drop=True)
                entry_exp=df2.merge(list_in[['Entry_Kegg','Exp']].rename(columns={'Entry_Kegg':'Entry'}),
                                    on='Entry',how='left')
                exp = float(entry_exp[['Exp']].count())
                if  exp > 0:
                    a='b'
                else:
                    entry_exp=entry_exp[['Entry']]
            else:
                entry_exp=df2

            entry_exp['Freq']=float(frame['Freq'].min()*.5)
            entry_exp['num']=float(frame['num'].count() + 1)
            nodes=pd.concat([frame,entry_exp],sort=False)
            nodes_file_name='nodes_KEGG_Enrichment_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_kegg+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/eduardo1011/Programas/master/Enrichment_Plots.R').content.decode()
            R_script_enrich = re.sub('qwertyuiop',level_1_kegg+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_kegg+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_kegg,R_script_enrich) # store plots
            R_script_enrich = re.sub('ASPECT','Pathways',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','Pathways Annotation (BP)',R_script_enrich) # legend
            f= open(level_1_kegg+'/Kegg_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\n*****  KEGG Pathways: No significant terms were found')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),background_info[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            lis['entry']='entry'
            lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
            lis=lis.groupby('entry')['ent'].sum().reset_index()
            www=pd.merge(list_input_match[['Entry_Kegg']],background_info,on='Entry_Kegg',how='left')
            uuu=www.groupby('GO')['Entry'].count().reset_index().sort_values(by ='Entry',ascending=False).reset_index(drop=True).drop_duplicates()
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+file_path+
                      '\nAssociation file name\t'+analysis+
                      '\nTotal number of background\t'+str(background_info['Entry'].drop_duplicates().count())+
                      '\nTotal number of list\t'+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with Pathways\t'+str(background_info['Entry'].drop_duplicates().count())+
                      '\nList input with Pathways\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+  # ****
                      '\nCorrection Method\t'+''.join(method_P)+ # ***
                      '\nValue\t'+str(User_value_P)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_P, information], axis=0, sort=False) # ***
            process_goa=results_process_P[['GO','entry']].rename(columns={'entry':'Entry'})

            enrich_P=enrich_P.rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
            process_goa=process_goa.rename(columns={'GO':'Path'})
            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Pathways_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant KEGG Pathways',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_goa.to_excel(writer,'Edges Pathways',index=False)
            writer.save()
    
    else:
        print('\n*****  KEGG Pathways: No significant terms were found')
        results_process_P=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_P=results_process_P[['entry']].count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        method_P='0'
        User_value_P='0'
            


# In[ ]:


## summary
database='KEGG'
aspect='Pathway'

information=[[database],[aspect],[''.join(method_P)],[str(User_value_P)],[GO_count_P],[proteins_count_P]]
summary = DataFrame(information)
summary=summary.T.reset_index(drop=True).rename(columns={0:'DATABASE',1:'TYPE',2:'US_METH',3:'US_VAL',4:'PATHWAYS',5:'UNIQ_PROT'})

print('\nENRICHMENT SUMMARY')
print('\n',summary)


# In[ ]:


#  API REST KEGG
# extrae id y vías de un organismo
    #http://rest.kegg.jp/list/pathway/pcs
#    path:pcs00010	Glycolysis / Gluconeogenesis - Penicillium rubens
#    path:pcs00020	Citrate cycle (TCA cycle) - Penicillium rubens
#    path:pcs00030	Pentose phosphate pathway - Penicillium rubens

# extrae todos los identificadores del kegg y su respectivo id de la via
    #http://rest.kegg.jp/link/pathway/pcs
#    pcs:Pc06g00180	path:pcs00010
#    pcs:Pc12g05620	path:pcs00010
#    pcs:Pc12g07100	path:pcs00010

# extrae identificador kegg y id del gen de NCBI
    #http://rest.kegg.jp/conv/pcs/ncbi-geneid
#    pcs:Pc01g00010	ncbi-geneid:8305880
#    pcs:Pc01g00020	ncbi-geneid:8305881
#    pcs:Pc01g00030	ncbi-geneid:8305882

# extrae identificadores kegg y id uniprot
    #http://rest.kegg.jp/conv/uniprot/pcs
#    pcs:Pc01g00010	up:B6GVP5
#    pcs:Pc01g00020	up:B6GVP6
#    pcs:Pc01g00030	up:B6GVP7

# extrae identificador kegg y id de la proteina de NCBI
    #http://rest.kegg.jp/conv/ncbi-proteinid/pcs
#    pcs:Pc01g00010	ncbi-proteinid:XP_002556565
#    pcs:Pc01g00020	ncbi-proteinid:XP_002556566
#    pcs:Pc01g00030	ncbi-proteinid:XP_002556567


from tkinter import * 
from tkinter.ttk import *
from tkinter import Tk, Label, Button, Radiobutton, IntVar
import tkinter as tk
import io
import base64
try:
    # Python2
    import Tkinter as tk
    from urllib2 import urlopen
except ImportError:
    # Python3
    import tkinter as tk
    from urllib.request import urlopen

def ask_multiple_choice_question(prompt, options):
    root = Tk()
    root.title("NeVOmics")
    root.geometry("850x500")
    root.overrideredirect(1)

    if prompt:
        Label(root, text="\nIf you use NeVOmics in your research, please cite:      ",font=("Arial", 11)).pack()
        Label(root, text="NeVOmics: an enrichment tool for gene ontology and functional\n"+
        "network analysis of data from OMICs technologies",font=("Arial", 11),bg="pale green").pack()
        Label(root, text="\n[ You want to create the network visualizations ]\n",font=("Arial", 15)).pack()
    v = IntVar()
    for i, option in enumerate(options):
        Radiobutton(root, text=option,font=("Arial", 13,"bold"),bg="gold2", variable=v, value=i).pack(anchor="n")
    Label(root, text="!!!!!   It may take several minutes   !!!!!",font=("Arial", 12)).pack()
    Button(root,text="          Submit          ",bg="black", fg="white",font=("Arial", 11), command=root.destroy).pack(pady=15) #anchor must be n, ne, e, se, s, sw, w, nw, or center
    image_url = "https://raw.githubusercontent.com/eduardo1011/Programas/master/NeVOmicsPlots.gif"
    image_byt = urlopen(image_url).read()
    image_b64 = base64.encodestring(image_byt)
    logo = tk.PhotoImage(data=image_b64)
    w1 = tk.Label(root, image=logo).pack()
    root.mainloop()						             #pady  separación entre la tercer linea y submit
    if v.get() == 0: return None
    return options[v.get()]

result = ask_multiple_choice_question(" ",
    [
        "   None     ",
        "   KEGG PATHWAYS       "
    ]
)

if format(repr(result)) == 'None':
    print('\n!!!!! Graphics not generated !!!!!')
else:
    folders = [level_1_kegg]
    # Biological process
    if ''.join(re.findall('KEGG',format(repr(result)))) == 'KEGG':
        # find R scripst process
        plots_selection=[]
        for i in folders:
            plots_selection.append([i+''.join(fnmatch.filter(os.listdir((i)), '*.R'))]) 
        # run R scripts
        for i in plots_selection:
            run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
            run_uni.wait()


# In[ ]:
if os.path.exists('data'): shutil.rmtree('data')
if os.path.exists("./*.RData"): os.remove("./*.RData")
if os.path.exists("./Rplots.pdf"): os.remove("./Rplots.pdf")
# print total time of analysis
lapso_total = datetime.now() - inicio_total
print('\n'+new_folder+': Analysis Time (hh:mm:ss.ms) {}'.format(lapso_total),'\n')

