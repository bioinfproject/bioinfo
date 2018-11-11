
# coding: utf-8

# In[432]:


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


import tkinter as tk
from tkinter import filedialog
## Control of input file
print('\n[ Step 1: Submit file (Uniprot IDs) ]')
root = tk.Tk()
root.withdraw()
input_file = filedialog.askopenfilename()
root.destroy()
if input_file == ():
    print('\n!!!!!!! File not found !!!!!!!')
    import tkinter as tk
    from tkinter import filedialog
    ## Control of input file
    print('\n[ Step 1: Submit file (Uniprot IDs) ]')
    root = tk.Tk()
    root.withdraw()
    input_file = filedialog.askopenfilename()
    root.destroy()
    if input_file == ():
        print('\n!!!!!!! File not found !!!!!!!')
        sys.exit()
    else:
        print('=====> : ',input_file)
else:
    print('=====> : ',input_file)
   
## Extract first id of protein
inp_file=pd.read_csv(''.join(input_file),sep='\t',header=None)
first_entry=inp_file[0].iloc[0]

##############################################################
################           Uniprot         ###################
##############################################################
## exttract id-organism
id_organism = requests.get("https://www.uniprot.org/uniprot/?query="+first_entry+"&sort=score&columns=organism-id,organism&format=tab&limit=1").content.decode()
## Kill process if not found id-organism
if id_organism == '':
    print('\n!!!!!!! ID-Organism not found, Check your identifiers list !!!!!!!')
    if os.path.exists(new_folder): shutil.rmtree(new_folder)
    if os.path.exists('data'): shutil.rmtree('data')
    if os.path.exists('GO.py'): os.remove('GO.py')
    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
    if os.path.exists('HD.py'): os.remove('HD.py')
    sys.exit()
else:
    ## id-organism
    Prefix=DataFrame(re.findall('[0-9]{1,30}',id_organism))[0].iloc[0]
    ## Organism identified
    strain=DataFrame(re.findall('[A-Za-z].*',id_organism))[0].iloc[1]
    print('\n-- ','Organism identified: ',strain)
    print('')
    
    ## Control of Uniprot and GOA directories
    ls=[]
    for i in os.listdir("./"):
        if re.search('job_GO_[0-9]{1,3}',str(i)):
            ls.append(i)
    if ls == []:
        new_folder='job_GO_1'
        xoxo='1'
        os.makedirs('job_GO_1/job_Uniprot_1/job_Uniprot_plots_1',exist_ok=True)
        os.makedirs('job_GO_1/job_GOA_1/job_GOA_plots_1',exist_ok=True)
        level_1_uniprot=new_folder+'/job_Uniprot_'+xoxo+'/'
        level_2_uniprot=new_folder+'/job_Uniprot_'+xoxo+'/job_Uniprot_plots_'+xoxo+'/'
        level_1_goa=new_folder+'/job_GOA_'+xoxo+'/'
        level_2_goa=new_folder+'/job_GOA_'+xoxo+'/job_GOA_plots_'+xoxo+'/'
    else:
        n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
        xoxo=str(n+1)
        n=str(n)
        old_folder=''.join(re.findall('job_GO_'+n,str(ls)))
        new_folder=re.sub(n,xoxo,old_folder)
        os.makedirs(new_folder+'/job_Uniprot_'+xoxo+'/job_Uniprot_plots_'+xoxo+'',exist_ok=True)
        os.makedirs(new_folder+'/job_GOA_'+xoxo+'/job_GOA_plots_'+xoxo+'',exist_ok=True)
        level_1_uniprot=new_folder+'/job_Uniprot_'+xoxo+'/'
        level_2_uniprot=new_folder+'/job_Uniprot_'+xoxo+'/job_Uniprot_plots_'+xoxo+'/'
        level_1_goa=new_folder+'/job_GOA_'+xoxo+'/'
        level_2_goa=new_folder+'/job_GOA_'+xoxo+'/job_GOA_plots_'+xoxo+'/'
    
## Create a folder
os.makedirs('data',exist_ok=True)

# save all ontology from go.obo file from GOC
#go0=BeautifulSoup(urlopen("http://purl.obolibrary.org/obo/go.obo"),"html.parser" ).decode()
import os, fnmatch
def find(pattern,path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
file_obo=find('go.obo', './')
if file_obo == []:
    url = 'http://purl.obolibrary.org/obo/go.obo'
    go_file = new_folder+'/go.obo'
    with open(go_file, 'wb') as f:
        #print ("Downloading %s" % file_name)
        response = requests.get(url, stream=True)
        total_length = response.headers.get('content-length')
        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\rLoading [%s%s]" % ('>' * done, ' ' * (50-done)) )    
                sys.stdout.flush()
    with open(new_folder+'/go.obo', 'r') as g:
        go0 = g.read()
else:
    if os.path.exists(file_obo[0]):
        with open(file_obo[0], 'r') as g:
            go0 = g.read()

go1=DataFrame(re.findall('\nid: GO:[0-9]{7}',go0),columns=['GO']).replace({'\nid: ':''},regex=True)
go2=DataFrame(re.findall('\nname:.*',go0),columns=['Term']).replace({'\nname: ':''},regex=True)
go3=DataFrame(re.findall('\nnamespace:.*',go0),columns=['Aspect']).replace({'\nnamespace:.':'',
                                                                            'biological_process':'P',
                                                                            'molecular_function':'F',
                                                                            'cellular_component':'C'},regex=True)
frames=[go1,go2,go3]
all_GO=pd.concat(frames,axis=1).dropna().reset_index(drop=True)
ww=all_GO['Term'].str.contains('molecular_function|cellular_component|biological_process')
all_annotations = all_GO[~ww].reset_index(drop=True) # without root aspects
#all_annotations.to_csv('./data/All_annotation.txt',sep='\t',index=None)
#
## GO description for all Process category
# Extract GO terms by category (C, F, P), description files
pro=[]
fun=[]
com=[]
for index, row in all_annotations.iterrows():
    if row['Aspect'] == 'P':
        a=(str(row['GO']),str(row['Term']),str(row['Aspect']))
        pro.append(a)
    if  row['Aspect'] == 'F':
        b=(str(row['GO']),str(row['Term']),str(row['Aspect']))
        fun.append(b)
    if  row['Aspect'] == 'C':
        c=(str(row['GO']),str(row['Term']),str(row['Aspect']))
        com.append(c)
GO_BP=pd.DataFrame(list(pro)).rename(columns={0:'GO',1:'Term',2:'Aspect'})
GO_BP[['GO','Term']].to_csv('./data/GO_BP.txt',sep='\t',index=None)
GO_MF=pd.DataFrame(list(fun)).rename(columns={0:'GO',1:'Term',2:'Aspect'})
GO_MF[['GO','Term']].to_csv('./data/GO_MF.txt',sep='\t',index=None)
GO_CC=pd.DataFrame(list(com)).rename(columns={0:'GO',1:'Term',2:'Aspect'})
GO_CC[['GO','Term']].to_csv('./data/GO_CC.txt',sep='\t',index=None)
##########

##############################################################
################           Uniprot         ###################
##############################################################

## Download all acc uniprot and GO IDs
import os, fnmatch
def find(pattern,path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
file_uniprot=find('annotation_'+Prefix, './')
if file_uniprot == []:
    annotation_go_uniprot=urllib.request.urlretrieve('https://www.uniprot.org/uniprot/?query=organism:'+Prefix+'&format=tab&columns=id,go-id', new_folder+'/annotation_'+Prefix)
    acc_uniprot_GO_id=pd.read_csv(new_folder+'/annotation_'+Prefix,sep='\t').rename(columns={'Gene Ontology IDs':'GO'}).dropna().reset_index(drop=True)
else:
    if os.path.exists(file_uniprot[0]):
        #print('ya existe')
        acc_uniprot_GO_id=pd.read_csv(file_uniprot[0],sep='\t').rename(columns={'Gene Ontology IDs':'GO'}).dropna().reset_index(drop=True)

## Edit information uniprot see "acc_uniprot_GO_id"
# convert data frame acc_uniprot_GO_id to two columns

#########   option 1
#aa=[]
#counter=0
#for go in acc_uniprot_GO_id['GO']:
#    ss=re.findall('GO:[0-9]{7}',go) # find all go by row
#    dd='\n\t'.join(ss) # convert list to str
#    ff='; '.join(ss)
#    gg=acc_uniprot_GO_id['Entry'][counter]
#    if dd == '': # counter of cycles
#        print('not content')
#    else:
#        counter += 1
#    aa.append(gg+'\t'+dd)
#hh='\n'.join(aa)
#Entry_GO_id=pd.read_csv(StringIO(hh),sep='\t',header=None,names=['Entry','GO']).fillna(method='pad')
#Entry_GO_id

#########   option 2 short version
b01=[]
for index, row in acc_uniprot_GO_id.iterrows():
    b02=re.findall('GO:[0-9]{7}',row['GO'])
    for i in b02:
        b01.append([row['Entry'],i])
Entry_GO_id=DataFrame(b01)[DataFrame(b01)[1] != ''].reset_index(drop=True).rename(columns={0:'Entry',1:'GO'})

#
## Open user's gene list 
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


# In[ ]:


####################  Enrichment Analysis  ####################

##############################################################
################       Uniprot GOA         ###################
##############################################################
## GOA Proteome Sets
## Download Uniprot GOA file
#print('\n----------------------------------------------------\n               GOA-Uniprot Annotation\n----------------------------------------------------\n\n               ■ Biological Process\n')
#all_goa=BeautifulSoup(urlopen("https://www.ebi.ac.uk/inc/drupal/goa/proteomes_release.html"),"html.parser" ).decode()
#goa_all=re.sub(r'<.*?>',' ',all_goa)
#x=str(re.findall('\n  '+Prefix+'.*.goa',goa_all))
#y = ''.join(re.findall('[0-9]{0,20}.[A-Z]_[a-z].*goa',x)) #convierte lista a string

################################
link = 'https://www.ebi.ac.uk/inc/drupal/goa/proteomes_release.html'
file_name = 'data/proteomes'
with open(file_name, 'wb') as f:
        #print ("Downloading %s" % file_name)
        response = requests.get(link, stream=True)
        total_length = response.headers.get('content-length')

        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                sys.stdout.write("\rLoading [%s%s]" % ('>' * done, ' ' * (50-done)) )    
                sys.stdout.flush()
with open('data/proteomes', 'r') as g:
    h = g.read()
a=re.sub(r'<.*?>',' ',h)
x=str(re.findall('  '+Prefix+'  '+strain[0]+'.*.goa',a))
y = ''.join(re.findall('[0-9]{0,20}.[A-Z]_[a-z].*goa',x)) #convierte lista a string
################################


if y == '':
    print('\n!!!!!!! ID-Organism not found in GOA-Uniprot (Complete Annotation). Enrichment will be done only with the Annotation of Uniprot Database !!!!!!!\n')
## descarga
else:
    # download goa file from database
    #out2 =urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/'+y, 'data/goa_file')
    import os, fnmatch
    def find(pattern,path):
        result = []
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result
    file_goa=find(y, './')
    if file_goa == []:
        #print('no está, entonces descargarlo')

        link = 'http://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/'+y+''
        file_name = new_folder+'/'+y
        with open(file_name, 'wb') as f:
            #print ("Downloading %s" % file_name)
            response = requests.get(link, stream=True)
            total_length = response.headers.get('content-length')
            if total_length is None: # no content length header
                f.write(response.content)
            else:
                dl = 0
                total_length = int(total_length)
                for data in response.iter_content(chunk_size=4096):
                    dl += len(data)
                    f.write(data)
                    done = int(50 * dl / total_length)
                    sys.stdout.write("\rLoading [%s%s]" % ('>' * done, ' ' * (50-done)) )  
                    sys.stdout.flush()
        ######################################
        orden_columnas_goa=[1,4,8]
        names=['Entry','GO','Aspect']
        with open(new_folder+'/'+y, 'r') as f3:
            my_gaf = f3.read()
            goa_file= pd.read_csv(StringIO(re.sub('!.*|\n!.*','',my_gaf)),sep='\t',header=None,usecols=orden_columnas_goa,names=names)

    else:
        if os.path.exists(file_goa[0]):
            #print('si existe')
            orden_columnas_goa=[1,4,8]
            names=['Entry','GO','Aspect']
            with open(file_goa[0], 'r') as f3:
                my_gaf = f3.read()
                goa_file= pd.read_csv(StringIO(re.sub('!.*|\n!.*','',my_gaf)),sep='\t',header=None,usecols=orden_columnas_goa,names=names)
    ######################################


# In[ ]:

## File with background column
## User's gackground column
if len(inp_file.columns) == 3:
    #>>>>>
    # 1.- Preparation of background
    background=list_input[['Background']].rename(columns={'Background':'Entry'}) ## change column name
    back_with_GO=pd.DataFrame.merge(background, goa_file[['Entry','GO']].drop_duplicates(), how="left", on='Entry').dropna()    ## User's background list
    back_with_GO[['Entry']].drop_duplicates().to_csv('data/Background.txt',index=None)

    # 2.- Preparation of list with GO terms
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),back_with_GO,how="left", on='Entry').dropna()
    list_input_match[['Entry']].drop_duplicates().to_csv('data/List.txt',index=None,)
    
    # 3.- background with: Entry	GO, for association file
    entry_go_term=pd.merge(back_with_GO,all_annotations,on='GO',how='left').dropna()
    entry_go_term[['Entry','GO']].to_csv('data/Association.txt',index=None,sep='\t')
    #>>>>>
    
    ########################################
    ## Enrichment analysis
    ########################################
    
    ########################################
    #*************************************** Biological Process
    ########################################
    
    #exploratory analysis of P-value in data
    file='GO_BP.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_BP.tsv'])
    enrich_P=pd.read_csv('data/Enrichment_analysis_GO_BP.tsv',sep='\t')
    
    if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
        print('\n\n----------------------------------------------------\n'+
              '                Uniprot-GOA Annotation\n----------------------------------------------------'+
              '\n\n*****  Biological Process:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        
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
            if os.path.exists('GO.py'): os.remove('GO.py')
            if os.path.exists('KEGG.py'): os.remove('KEGG.py')
            if os.path.exists('HD.py'): os.remove('HD.py')
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
                if os.path.exists('GO.py'): os.remove('GO.py')
                if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                if os.path.exists('HD.py'): os.remove('HD.py')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_P=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
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
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_P=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        if os.path.exists('GO.py'): os.remove('GO.py')
                        if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                        if os.path.exists('HD.py'): os.remove('HD.py')
                        sys.exit()
        ########################################
        ########################################
        analysis='GO_BP.txt'
        if ''.join(method_P) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
            enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
            results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
            proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_P=results_process_P[['GO']].count()[0]

        else:
            if ''.join(method_P) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
                enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
                results_process_P=enrich_P[enrich_P.Sig == 'T']
                proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_P=results_process_P[['GO']].count()[0]
        
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_P['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_P=results_process_P.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_P['Bonf_corr'][0]/results_process_P['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_P)+
                      '\nValue\t'+str(User_value_P)+
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_P, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_P.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            process_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_GOA_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            process_goa.to_csv(level_1_goa+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_goa.to_excel(writer,'Edges Process',index=False)
            writer.save()

            frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=process_goa[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_GOA_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_goa+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_goa+'Process')
            R_script_enrich = re.sub('qwertyuiop',level_1_goa+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_goa+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_goa+'Process/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (BP)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Biological Process',R_script_enrich)
            f= open(level_1_goa+'/Process_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Biological Process (GOA)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_P['Bonf_corr'][0]/enrich_P['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_P)+ # ***
                      '\nValue\t'+str(User_value_P)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    
            combine=pd.concat([results_process_P, information], axis=0, sort=False) # ***
            
            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
 
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°                
    else:
        print('\n*****  Biological Process: No significant terms were found')
        results_process_P=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_P=results_process_P[['entry']].count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        User_value_P='0'
        method_P='0'
     
    ########################################
    #*************************************** Molecular Function
    ########################################    
    
    #exploratory analysis of P-value in data
    file='GO_MF.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_MF.tsv'])
    enrich_F=pd.read_csv('data/Enrichment_analysis_GO_MF.tsv',sep='\t')
    
    if enrich_F[(enrich_F.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Molecular Function:  ',enrich_F[(enrich_F.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        
        ## Choice of correction corection method
        methods = ['Bonferroni','FDR']
        User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        method_F=[]
        for x in methods:
            if x == User_method:
                method_F.append(x)
        if method_F == []:
            print('\n!!!!! Incorrect Method !!!!!')
            User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
            for x in methods:
                if x == User_method:
                    method_F.append(x)
        if method_F == []:
            print('\n!!!!! Incorrect Method !!!!!')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('data'): shutil.rmtree('data')
            if os.path.exists('GO.py'): os.remove('GO.py')
            if os.path.exists('KEGG.py'): os.remove('KEGG.py')
            if os.path.exists('HD.py'): os.remove('HD.py')
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
                if os.path.exists('GO.py'): os.remove('GO.py')
                if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                if os.path.exists('HD.py'): os.remove('HD.py')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_F=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_F=float(a)
            else:
                print('Incorrect value')
                a=input('[ Must be in this range: 0.005 - 0.5 ]\n=====> : ')
                if a == '':
                    a='aaa'
                if re.search(r'[A-Za-z]{1,10}',a):
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_F=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        if os.path.exists('GO.py'): os.remove('GO.py')
                        if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                        if os.path.exists('HD.py'): os.remove('HD.py')
                        sys.exit()
        ########################################
        ########################################
        analysis='GO_MF.txt'
        if ''.join(method_F) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
            enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
            results_process_F=enrich_F[enrich_F.Bonf_corr < User_value_F]
            proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_F=results_process_F[['GO']].count()[0]
        else:
            if ''.join(method_F) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
                enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
                results_process_F=enrich_F[enrich_F.Sig == 'T']
                proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_F=results_process_F[['GO']].count()[0]
        
                    
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_F['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_F=results_process_F.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_F['Bonf_corr'][0]/results_process_F['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_F)+
                      '\nValue\t'+str(User_value_F)+
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_F, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_F.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            function_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_GOA_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            function_goa.to_csv(level_1_goa+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False)
            function_goa.to_excel(writer,'Edges Function',index=False)
            writer.save()

            frame=results_process_F[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=function_goa[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_GOA_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_goa+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_goa+'Function')
            R_script_enrich = re.sub('qwertyuiop',level_1_goa+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_goa+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_goa+'Function/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (MF)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Molecular Function',R_script_enrich)
            f= open(level_1_goa+'/Function_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Molecular Function (GOA)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(enrich_F['Bonf_corr'][0]/enrich_F['P'][0]))+  # ****
                    '\nCorrection Method\t'+''.join(method_F)+ # ***
                    '\nValue\t'+str(User_value_F)+ # ***
                    '\n\t\n'+
                    '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    
            combine=pd.concat([results_process_F, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
                   
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

    else:
        print('\n*****  Molecular Function: No significant terms were found')
        results_process_F=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_F=results_process_F[['entry']].count()[0]
        GO_count_F=results_process_F[['GO']].count()[0]
        User_value_F='0'
        method_F='0'
        
    ########################################
    #*************************************** Cellular Component
    ########################################
                     
    #exploratory analysis of P-value in data
    file='GO_CC.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_CC.tsv'])
    enrich_C=pd.read_csv('data/Enrichment_analysis_GO_CC.tsv',sep='\t')
    
    if enrich_C[(enrich_C.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Cellular Component:  ',enrich_C[(enrich_C.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        
        ## Choice of correction corection method
        methods = ['Bonferroni','FDR']
        User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        method_C=[]
        for x in methods:
            if x == User_method:
                method_C.append(x)
        if method_C == []:
            print('\n!!!!! Incorrect Method !!!!!')
            User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
            for x in methods:
                if x == User_method:
                    method_C.append(x)
        if method_C == []:
            print('\n!!!!! Incorrect Method !!!!!')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('data'): shutil.rmtree('data')
            if os.path.exists('GO.py'): os.remove('GO.py')
            if os.path.exists('KEGG.py'): os.remove('KEGG.py')
            if os.path.exists('HD.py'): os.remove('HD.py')
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
                if os.path.exists('GO.py'): os.remove('GO.py')
                if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                if os.path.exists('HD.py'): os.remove('HD.py')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_C=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_C=float(a)
            else:
                print('Incorrect value')
                a=input('[ Must be in this range: 0.005 - 0.5 ]\n=====> : ')
                if a == '':
                    a='aaa'
                if re.search(r'[A-Za-z]{1,10}',a):
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_C=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        if os.path.exists('GO.py'): os.remove('GO.py')
                        if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                        if os.path.exists('HD.py'): os.remove('HD.py')
                        sys.exit()
        ########################################
        ########################################
        analysis='GO_CC.txt'
        if ''.join(method_C) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_C),'data/Enrichment_Analysis_CC.csv'])
            enrich_C=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
            results_process_C=enrich_C[enrich_C.Bonf_corr < User_value_C]
            proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_C=results_process_C[['GO']].count()[0]

        else:
            if ''.join(method_C) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_C),'data/Enrichment_Analysis_CC.csv'])
                enrich_C=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
                results_process_C=enrich_C[enrich_C.Sig == 'T']
                proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_C=results_process_C[['GO']].count()[0]
        
                    
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_C['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_C=results_process_C.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_C['Bonf_corr'][0]/results_process_C['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_C)+
                      '\nValue\t'+str(User_value_C)+
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_C, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_C.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            component_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_GOA_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            component_goa.to_csv(level_1_goa+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False)
            component_goa.to_excel(writer,'Edges Component',index=False)
            writer.save()

            frame=results_process_C[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=component_goa[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_GOA_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_goa+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_goa+'Component')
            R_script_enrich = re.sub('qwertyuiop',level_1_goa+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_goa+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_goa+'Component/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (CC)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Cellular Component',R_script_enrich)
            f= open(level_1_goa+'/Component_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Cellular Component (GOA)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(enrich_C['Bonf_corr'][0]/enrich_C['P'][0]))+  # ****
                    '\nCorrection Method\t'+''.join(method_C)+ # ***
                    '\nValue\t'+str(User_value_C)+ # ***
                    '\n\t\n'+
                    '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    
            combine=pd.concat([results_process_C, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
 
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

    else:
        print('\n*****  Cellular Component: No significant terms were found')
        results_process_C=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_C=results_process_C[['entry']].count()[0]
        GO_count_C=results_process_C[['GO']].count()[0]
        User_value_C='0'
        method_C='0'

########################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## File without background column
else:
    print('\n\n* [ The analysis will be done using the whole proteome ]')
    # 1.- Background complete (proteome)
    goa_file[['Entry']].drop_duplicates().to_csv('data/Background.txt',index=None)
    background=goa_file[['Entry']]
    back_with_GO=goa_file[['Entry']]
    
    
    # 2.- Preparation of list with GO terms
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file, how="left", on='Entry').dropna()
    list_input_match[['Entry']].drop_duplicates().to_csv('data/List.txt',index=False, float_format='%.0f')    ## Build of kegg assotiation file
    
    # 3.- All proteins with: Entry	GO, for association file
    goa_file[['Entry','GO']].drop_duplicates().to_csv('data/Association.txt',index=None,sep='\t')
    #>>>>> 
                    
    ########################################
    ## Enrichment analysis
    ########################################
    
    ########################################
    #*************************************** Biological Process
    ########################################
    
    #exploratory analysis of P-value in data
    file='GO_BP.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_BP.tsv'])
    enrich_P=pd.read_csv('data/Enrichment_analysis_GO_BP.tsv',sep='\t')
    print('\n----------------------------------------------------\n'+
              '                Uniprot-GOA Annotation\n----------------------------------------------------')
    
    if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Biological Process:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        
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
            if os.path.exists('GO.py'): os.remove('GO.py')
            if os.path.exists('KEGG.py'): os.remove('KEGG.py')
            if os.path.exists('HD.py'): os.remove('HD.py')
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
                if os.path.exists('GO.py'): os.remove('GO.py')
                if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                if os.path.exists('HD.py'): os.remove('HD.py')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_P=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
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
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_P=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        if os.path.exists('GO.py'): os.remove('GO.py')
                        if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                        if os.path.exists('HD.py'): os.remove('HD.py')
                        sys.exit()
        ########################################
        ########################################
        analysis='GO_BP.txt'
        if ''.join(method_P) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
            enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
            results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
            proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_P=results_process_P[['GO']].count()[0]

        else:
            if ''.join(method_P) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
                enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
                results_process_P=enrich_P[enrich_P.Sig == 'T']
                proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_P=results_process_P[['GO']].count()[0]
                
                    
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_P['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_P=results_process_P.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_P['Bonf_corr'][0]/results_process_P['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_P)+
                      '\nValue\t'+str(User_value_P)+
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_P, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_P.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            process_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_GOA_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            process_goa.to_csv(level_1_goa+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_goa.to_excel(writer,'Edges Process',index=False)
            writer.save()

            frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=process_goa[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_GOA_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_goa+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_goa+'Process')
            R_script_enrich = re.sub('qwertyuiop',level_1_goa+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_goa+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_goa+'Process/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (BP)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Biological Process',R_script_enrich)
            f= open(level_1_goa+'/Process_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        
        else:
            print('\nThere are not enrichment terms for Biological Process (GOA)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_P['Bonf_corr'][0]/enrich_P['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_P)+ # ***
                      '\nValue\t'+str(User_value_P)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    
            combine=pd.concat([results_process_P, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                        
    else:
        print('\n*****  Biological Process: No significant terms were found')
        results_process_P=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_P=results_process_P[['entry']].count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        User_value_P='0'
        method_P='0'
     
    ########################################
    #*************************************** Molecular Function
    ########################################    
    
    #exploratory analysis of P-value in data
    file='GO_MF.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_MF.tsv'])
    enrich_F=pd.read_csv('data/Enrichment_analysis_GO_MF.tsv',sep='\t')
    
    if enrich_F[(enrich_F.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Molecular Function:  ',enrich_F[(enrich_F.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        
        ## Choice of correction corection method
        methods = ['Bonferroni','FDR']
        User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        method_F=[]
        for x in methods:
            if x == User_method:
                method_F.append(x)
        if method_F == []:
            print('\n!!!!! Incorrect Method !!!!!')
            User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
            for x in methods:
                if x == User_method:
                    method_F.append(x)
        if method_F == []:
            print('\n!!!!! Incorrect Method !!!!!')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('data'): shutil.rmtree('data')
            if os.path.exists('GO.py'): os.remove('GO.py')
            if os.path.exists('KEGG.py'): os.remove('KEGG.py')
            if os.path.exists('HD.py'): os.remove('HD.py')
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
                if os.path.exists('GO.py'): os.remove('GO.py')
                if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                if os.path.exists('HD.py'): os.remove('HD.py')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_F=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_F=float(a)
            else:
                print('Incorrect value')
                a=input('[ Must be in this range: 0.005 - 0.5 ]\n=====> : ')
                if a == '':
                    a='aaa'
                if re.search(r'[A-Za-z]{1,10}',a):
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_F=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        if os.path.exists('GO.py'): os.remove('GO.py')
                        if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                        if os.path.exists('HD.py'): os.remove('HD.py')
                        sys.exit()
        ########################################
        ########################################
        analysis='GO_MF.txt'
        if ''.join(method_F) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
            enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
            results_process_F=enrich_F[enrich_F.Bonf_corr < User_value_F]
            proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_F=results_process_F[['GO']].count()[0]
        else:
            if ''.join(method_F) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
                enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
                results_process_F=enrich_F[enrich_F.Sig == 'T']
                proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_F=results_process_F[['GO']].count()[0]
        
                    
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_F['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_F=results_process_F.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_F['Bonf_corr'][0]/results_process_F['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_F)+
                      '\nValue\t'+str(User_value_F)+
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_F, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_F.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            function_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_GOA_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            function_goa.to_csv(level_1_goa+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False)
            function_goa.to_excel(writer,'Edges Function',index=False)
            writer.save()

            frame=results_process_F[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=function_goa[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_GOA_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_goa+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_goa+'Function')
            R_script_enrich = re.sub('qwertyuiop',level_1_goa+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_goa+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_goa+'Function/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (MF)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Molecular Function',R_script_enrich)
            f= open(level_1_goa+'/Function_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Molecular Function (GOA)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(enrich_F['Bonf_corr'][0]/enrich_F['P'][0]))+  # ****
                    '\nCorrection Method\t'+''.join(method_F)+ # ***
                    '\nValue\t'+str(User_value_F)+ # ***
                    '\n\t\n'+
                    '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    
            combine=pd.concat([results_process_F, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

    else:
        print('\n*****  Molecular Function: No significant terms were found')
        results_process_F=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_F=results_process_F[['entry']].count()[0]
        GO_count_F=results_process_F[['GO']].count()[0]
        User_value_F='0'
        method_F='0'
        
    ########################################
    #*************************************** Cellular Component
    ########################################
                     
    #exploratory analysis of P-value in data
    file='GO_CC.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_CC.tsv'])
    enrich_C=pd.read_csv('data/Enrichment_analysis_GO_CC.tsv',sep='\t')
    
    if enrich_C[(enrich_C.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Cellular Component:  ',enrich_C[(enrich_C.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        
        ## Choice of correction corection method
        methods = ['Bonferroni','FDR']
        User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        method_C=[]
        for x in methods:
            if x == User_method:
                method_C.append(x)
        if method_C == []:
            print('\n!!!!! Incorrect Method !!!!!')
            User_method=input('\n[ Step 2: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
            for x in methods:
                if x == User_method:
                    method_C.append(x)
        if method_C == []:
            print('\n!!!!! Incorrect Method !!!!!')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('data'): shutil.rmtree('data')
            if os.path.exists('GO.py'): os.remove('GO.py')
            if os.path.exists('KEGG.py'): os.remove('KEGG.py')
            if os.path.exists('HD.py'): os.remove('HD.py')
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
                if os.path.exists('GO.py'): os.remove('GO.py')
                if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                if os.path.exists('HD.py'): os.remove('HD.py')
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_C=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_C=float(a)
            else:
                print('Incorrect value')
                a=input('[ Must be in this range: 0.005 - 0.5 ]\n=====> : ')
                if a == '':
                    a='aaa'
                if re.search(r'[A-Za-z]{1,10}',a):
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists('GO.py'): os.remove('GO.py')
                    if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                    if os.path.exists('HD.py'): os.remove('HD.py')
                    sys.exit()
                else:
                    if 0.005 <= float(a) <= 0.5:
                        User_value_C=float(a)
                    else:
                        print('Incorrect value')
                        if os.path.exists(new_folder): shutil.rmtree(new_folder)
                        if os.path.exists('data'): shutil.rmtree('data')
                        if os.path.exists('GO.py'): os.remove('GO.py')
                        if os.path.exists('KEGG.py'): os.remove('KEGG.py')
                        if os.path.exists('HD.py'): os.remove('HD.py')
                        sys.exit()
        ########################################
        ########################################
        analysis='GO_CC.txt'
        if ''.join(method_C) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_C),'data/Enrichment_Analysis_CC.csv'])
            enrich_C=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
            results_process_C=enrich_C[enrich_C.Bonf_corr < User_value_C]
            proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_C=results_process_C[['GO']].count()[0]

        else:
            if ''.join(method_C) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_C),'data/Enrichment_Analysis_CC.csv'])
                enrich_C=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
                results_process_C=enrich_C[enrich_C.Sig == 'T']
                proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_C=results_process_C[['GO']].count()[0]
        
                    
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_C['GO'].count() >= 1:
            ## Information about analysis for GOA
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_C=results_process_C.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_C['Bonf_corr'][0]/results_process_C['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_C)+
                      '\nValue\t'+str(User_value_C)+
                      '\n\t\n'+
                      '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_C, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_C.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            component_goa=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_GOA_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            component_goa.to_csv(level_1_goa+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False)
            component_goa.to_excel(writer,'Edges Component',index=False)
            writer.save()

            frame=results_process_C[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=component_goa[['Entry']]

             #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_GOA_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_goa+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_goa+'Component')
            R_script_enrich = re.sub('qwertyuiop',level_1_goa+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_goa+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_goa+'Component/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (CC)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Cellular Component',R_script_enrich)
            f= open(level_1_goa+'/Component_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Cellular Component (GOA)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(enrich_C['Bonf_corr'][0]/enrich_C['P'][0]))+  # ****
                    '\nCorrection Method\t'+''.join(method_C)+ # ***
                    '\nValue\t'+str(User_value_C)+ # ***
                    '\n\t\n'+
                    '\nProteins with no information in Uniprot-GOA\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    
            combine=pd.concat([results_process_C, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_GOA_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
 
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

    else:
        print('\n*****  Cellular Component: No significant terms were found')
        results_process_C=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_C=results_process_C[['entry']].count()[0]
        GO_count_C=results_process_C[['GO']].count()[0]
        User_value_C='0'
        method_C='0'


# In[ ]:


## loop for goa
database=['Uniprot-GOA','Uniprot-GOA','Uniprot-GOA']
aspect=['P','F','C']
protein_count=[proteins_count_P,proteins_count_F,proteins_count_C]
user_methods=[''.join(method_P),''.join(method_F),''.join(method_C)]
descriptions=['GO_BP.txt','GO_MF.txt','GO_CC.txt']
user_values=[str(User_value_P),str(User_value_F),str(User_value_C)]
GOs_count=[GO_count_P,GO_count_F,GO_count_C]

summary_GOA=[]
for i , j , k , l , m , q in zip(database, aspect, user_methods, user_values, GOs_count, protein_count):
    summary_GOA.append([i ,j , k ,l ,m, q])
    #print(i , j , k , l , m)


# In[ ]:


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##### Uniprot
#
if len(inp_file.columns) == 3:
    
    # first clean variables
    background = back_with_GO =  list_input_match = entry_go_term = enrich_P = enrich_F = enrich_C =     results_process_P = proteins_count_P  = GO_count_P = results_process_F =     proteins_count_F  = GO_count_F = results_process_C = proteins_count_C  = GO_count_C = ''

    # 1.- Preparation of background
    background=list_input[['Background']].rename(columns={'Background':'Entry'}) ## change column name
    back_with_GO=pd.DataFrame.merge(background, Entry_GO_id[['Entry','GO']].drop_duplicates(), how="left", on='Entry').dropna()    ## User's background list
    back_with_GO[['Entry']].drop_duplicates().to_csv('data/Background.txt',index=None)
    
    # 2.- Preparation of list with GO terms
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),back_with_GO,how="left", on='Entry').dropna()
    list_input_match[['Entry']].drop_duplicates().to_csv('data/List.txt',index=None,)
    
    # 3.- background with: Entry	GO, for association file
    entry_go_term=pd.merge(back_with_GO,all_annotations,on='GO',how='left').dropna()
    entry_go_term[['Entry','GO']].to_csv('data/Association.txt',index=None,sep='\t')
    #>>>>>
    
    ########################################
    ## Enrichment analysis
    ########################################
        
    ########################################
    #*************************************** Biological Process
    ########################################
        
    #exploratory analysis of P-value in data
    file='GO_BP.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_BP.tsv'])
    enrich_P=pd.read_csv('data/Enrichment_analysis_GO_BP.tsv',sep='\t')
    print('\n----------------------------------------------------\n'+
    '                UniprotKB Annotation\n----------------------------------------------------')
        
    if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Biological Process:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        analysis='GO_BP.txt'
        if ''.join(method_P) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
            enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
            results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
            proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_P=results_process_P[['GO']].count()[0]

        else:
            if ''.join(method_P) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
                enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
                results_process_P=enrich_P[enrich_P.Sig == 'T']
                proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_P=results_process_P[['GO']].count()[0]
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_P['GO'].count() >= 1:
            ## Information about analysis for Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_P=results_process_P.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_P['Bonf_corr'][0]/results_process_P['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_P)+
                      '\nValue\t'+str(User_value_P)+
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+str(lis['ent'][0])]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_P, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_P.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            process_uniprot=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_Uniprot_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            process_uniprot.to_csv(level_1_uniprot+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_uniprot.to_excel(writer,'Edges Process',index=False)
            writer.save()

            frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=process_uniprot[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_uniprot+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_uniprot+'Process')
            R_script_enrich = re.sub('qwertyuiop',level_1_uniprot+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_uniprot+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_uniprot+'Process/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (BP)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Biological Process',R_script_enrich)
            f= open(level_1_uniprot+'/Process_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Biological Process (Uniprot)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_P['Bonf_corr'][0]/enrich_P['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_P)+ # ***
                      '\nValue\t'+str(User_value_P)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_P, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                            
    else:
        print('\n*****  Biological Process: No significant terms were found')
        results_process_P=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_P=results_process_P[['entry']].count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        User_value_P='0'
        method_P='0'

    ########################################
    #*************************************** Molecular Function
    ########################################
        
    #exploratory analysis of P-value in data
    file='GO_MF.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_MF.tsv'])
    enrich_F=pd.read_csv('data/Enrichment_analysis_GO_MF.tsv',sep='\t')

    if enrich_F[(enrich_F.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Molecular Function:  ',enrich_F[(enrich_F.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        analysis='GO_MF.txt'
        if ''.join(method_F) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
            enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
            results_process_F=enrich_F[enrich_F.Bonf_corr < User_value_F]
            proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_F=results_process_F[['GO']].count()[0]

        else:
            if ''.join(method_F) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
                enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
                results_process_F=enrich_F[enrich_F.Sig == 'T']
                proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_F=results_process_F[['GO']].count()[0]
                
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_F['GO'].count() >= 1:
            ## Information about analysis for Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_F=results_process_F.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_F['Bonf_corr'][0]/results_process_F['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_F)+
                      '\nValue\t'+str(User_value_F)+
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_F, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_F.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            function_uniprot=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})
    
            ## save file with edges for graph
            edges_file_name='edges_Uniprot_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            function_uniprot.to_csv(level_1_uniprot+edges_file_name,index=None)
    
            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False)
            function_uniprot.to_excel(writer,'Edges Function',index=False)
            writer.save()

            frame=results_process_F[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=function_uniprot[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_uniprot+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_uniprot+'Function')
            R_script_enrich = re.sub('qwertyuiop',level_1_uniprot+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_uniprot+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_uniprot+'Function/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (MF)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Molecular Function',R_script_enrich)
            f= open(level_1_uniprot+'/Function_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Molecular Function (Uniprot)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_F['Bonf_corr'][0]/enrich_F['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_F)+ # ***
                      '\nValue\t'+str(User_value_F)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_F, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                            
    else:
        print('\n*****  Molecular Function: No significant terms were found')
        results_process_F=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_F=results_process_F[['entry']].count()[0]
        GO_count_F=results_process_F[['GO']].count()[0]
        User_value_F='0'
        method_F='0'


    ########################################
    #*************************************** Cellular Component
    ########################################
        
    #exploratory analysis of P-value in data
    file='GO_CC.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_CC.tsv'])
    enrich_C=pd.read_csv('data/Enrichment_analysis_GO_CC.tsv',sep='\t')

    if enrich_C[(enrich_C.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Celular Component:  ',enrich_C[(enrich_C.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        analysis='GO_CC.txt'
        if ''.join(method_C) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_C),'data/Enrichment_Analysis_CC.csv'])
            enrich_CC=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
            results_process_C=enrich_C[enrich_C.Bonf_corr < User_value_C]
            proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_C=results_process_C[['GO']].count()[0]

        else:
            if ''.join(method_C) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_F),'data/Enrichment_Analysis_CC.csv'])
                enrich_C=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
                results_process_C=enrich_C[enrich_C.Sig == 'T']
                proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_C=results_process_C[['GO']].count()[0]
                
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_C['GO'].count() >= 1:
            ## Information about analysis for Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_C=results_process_C.reset_index(drop=True)
            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(results_process_C['Bonf_corr'][0]/results_process_C['P'][0]))+
                    '\nCorrection Method\t'+''.join(method_C)+
                    '\nValue\t'+str(User_value_C)+
                    '\n\t\n'+
                    '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])

            combine=pd.concat([results_process_C, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_C.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            component_uniprot=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_Uniprot_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            component_uniprot.to_csv(level_1_uniprot+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False)
            component_uniprot.to_excel(writer,'Edges Component',index=False)
            writer.save()

            frame=results_process_C[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=component_uniprot[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_uniprot+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_uniprot+'Component')
            R_script_enrich = re.sub('qwertyuiop',level_1_uniprot+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_uniprot+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_uniprot+'Component/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (CC)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Cellular Component',R_script_enrich)
            f= open(level_1_uniprot+'/Component_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Cellular Component (Uniprot)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_C['Bonf_corr'][0]/enrich_C['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_C)+ # ***
                      '\nValue\t'+str(User_value_C)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_C, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        
    else:
        print('\n*****  Cellular Component: No significant terms were found')
        results_process_C=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_C=results_process_C[['entry']].count()[0]
        GO_count_C=results_process_C[['GO']].count()[0]
        User_value_C='0'
        method_C = '0'

## File without background column
else:
    # first clean variables
    background =  list_input_match = entry_go_term = enrich_P = enrich_F = enrich_C =     results_process_P = proteins_count_P  = GO_count_P = results_process_F =     proteins_count_F  = GO_count_F = results_process_C = proteins_count_C  = GO_count_C = ''
    
    print('\n* [ The analysis will be done using the whole proteome ]')
    
    ## File without background column
    # 1.- Background complete (proteome)
    background=Entry_GO_id[['Entry']].drop_duplicates()
    back_with_GO=Entry_GO_id[['Entry']].drop_duplicates()
    back_with_GO.to_csv('data/Background.txt',index=None)
        
    # 2.- Preparation of list with GO terms
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id, how="left", on='Entry').dropna()
    list_input_match[['Entry']].drop_duplicates().to_csv('data/List.txt',index=False, float_format='%.0f')    ## Build of kegg assotiation file
        
    # 3.- All proteins with: Entry	GO, for association file
    Entry_GO_id[['Entry','GO']].drop_duplicates().to_csv('data/Association.txt',index=None,sep='\t')
    #>>>>>
    
    ########################################
    ## Enrichment analysis
    ########################################
        
    ########################################
    #*************************************** Biological Process
    ########################################
        
    #exploratory analysis of P-value in data
    file='GO_BP.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_BP.tsv'])
    enrich_P=pd.read_csv('data/Enrichment_analysis_GO_BP.tsv',sep='\t')
    print('\n----------------------------------------------------\n'+
    '                UniprotKB Annotation\n----------------------------------------------------')
        
    if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Biological Process:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        analysis='GO_BP.txt'
        if ''.join(method_P) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
            enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
            results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
            proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_P=results_process_P[['GO']].count()[0]

        else:
            if ''.join(method_P) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_P),'data/Enrichment_Analysis_BP.csv'])
                enrich_P=pd.read_csv('data/Enrichment_Analysis_BP.csv',sep='\t')
                results_process_P=enrich_P[enrich_P.Sig == 'T']
                proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_P=results_process_P[['GO']].count()[0]
        
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_P['GO'].count() >= 1:
            ## Information about analysis for Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_P=results_process_P.reset_index(drop=True)
            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(results_process_P['Bonf_corr'][0]/results_process_P['P'][0]))+
                      '\nCorrection Method\t'+''.join(method_P)+
                      '\nValue\t'+str(User_value_P)+
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+str(lis['ent'][0])]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_P, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_P.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            process_uniprot=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_Uniprot_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            process_uniprot.to_csv(level_1_uniprot+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False)
            process_uniprot.to_excel(writer,'Edges Process',index=False)
            writer.save()

            frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=process_uniprot[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_uniprot+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_uniprot+'Process')
            R_script_enrich = re.sub('qwertyuiop',level_1_uniprot+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_uniprot+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_uniprot+'Process/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (BP)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Biological Process',R_script_enrich)
            f= open(level_1_uniprot+'/Process_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('There are not enrichment terms for Biological Process (Uniprot)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_P['Bonf_corr'][0]/enrich_P['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_P)+ # ***
                      '\nValue\t'+str(User_value_P)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_P, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Process_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_P.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                            
    else:
        print('\n*****  Biological Process: No significant terms were found')
        results_process_P=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_P=results_process_P[['entry']].count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        User_value_P='0'
        method_P='0'

    ########################################
    #*************************************** Molecular Function
    ########################################
        
    #exploratory analysis of P-value in data
    file='GO_MF.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_MF.tsv'])
    enrich_F=pd.read_csv('data/Enrichment_analysis_GO_MF.tsv',sep='\t')

    if enrich_F[(enrich_F.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Molecular Function:  ',enrich_F[(enrich_F.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        analysis='GO_MF.txt'
        if ''.join(method_F) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
            enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
            results_process_F=enrich_F[enrich_F.Bonf_corr < User_value_F]
            proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_F=results_process_F[['GO']].count()[0]

        else:
            if ''.join(method_F) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_F),'data/Enrichment_Analysis_MF.csv'])
                enrich_F=pd.read_csv('data/Enrichment_Analysis_MF.csv',sep='\t')
                results_process_F=enrich_F[enrich_F.Sig == 'T']
                proteins_count_F=DataFrame(results_process_F['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_F=results_process_F[['GO']].count()[0]
        
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_F['GO'].count() >= 1:
            ## Information about analysis for Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_F=results_process_F.reset_index(drop=True)
            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(results_process_F['Bonf_corr'][0]/results_process_F['P'][0]))+
                    '\nCorrection Method\t'+''.join(method_F)+
                    '\nValue\t'+str(User_value_F)+
                    '\n\t\n'+
                    '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_F, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_F.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            function_uniprot=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_Uniprot_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            function_uniprot.to_csv(level_1_uniprot+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False)
            function_uniprot.to_excel(writer,'Edges Function',index=False)
            writer.save()

            frame=results_process_F[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=function_uniprot[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_uniprot+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_uniprot+'Function')
            R_script_enrich = re.sub('qwertyuiop',level_1_uniprot+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_uniprot+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_uniprot+'Function/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (MF)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Molecular Function',R_script_enrich)
            f= open(level_1_uniprot+'/Function_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
            
            #-------------------
        else:
            print('\nThere are not enrichment terms for Molecular Function (Uniprot)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_F['Bonf_corr'][0]/enrich_F['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_F)+ # ***
                      '\nValue\t'+str(User_value_F)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_F, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Function_'+''.join(method_F)+'_'+str(User_value_F)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_F.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
                            
    else:
        print('\n*****  Molecular Function: No significant terms were found')
        results_process_F=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_F=results_process_F[['entry']].count()[0]
        GO_count_F=results_process_F[['GO']].count()[0]
        User_value_F='0'
        method_F='0'

    ########################################
    #*************************************** Cellular Component
    ########################################
        
    #exploratory analysis of P-value in data
    file='GO_CC.txt'
    fdr_min=0.05
    subprocess.call(["python",'HD.py',file,str(fdr_min),'data/Enrichment_analysis_GO_CC.tsv'])
    enrich_C=pd.read_csv('data/Enrichment_analysis_GO_CC.tsv',sep='\t')

    if enrich_C[(enrich_C.P < 0.05)]['P'].count() >= 1:
        print('\n*****  Cellular Component:  ',enrich_C[(enrich_C.P < 0.05)]['P'].count(),
              'GO terms were found with P-value < 0.05')
        analysis='GO_CC.txt'
        if ''.join(method_C) == 'Bonferroni':
            subprocess.call(["python",'HD.py',
                             analysis,str(User_value_C),'data/Enrichment_Analysis_CC.csv'])
            enrich_CC=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
            results_process_C=enrich_C[enrich_C.Bonf_corr < User_value_C]
            proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
            GO_count_C=results_process_C[['GO']].count()[0]

        else:
            if ''.join(method_C) == 'FDR':
                subprocess.call(["python",'HD.py',
                                 analysis,str(User_value_F),'data/Enrichment_Analysis_CC.csv'])
                enrich_C=pd.read_csv('data/Enrichment_Analysis_CC.csv',sep='\t')
                results_process_C=enrich_C[enrich_C.Sig == 'T']
                proteins_count_C=DataFrame(results_process_C['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
                GO_count_C=results_process_C[['GO']].count()[0]
        
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        if results_process_C['GO'].count() >= 1:
            ## Information about analysis for Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()
            results_process_C=results_process_C.reset_index(drop=True)
            report = ['\n\t\n'+
                    'NeVOmics\t'+new_folder+
                    '\n\nInput file name\t'+''.join(input_file)+
                    '\nAssociation file name\t'+analysis+
                    "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                    "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                    '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                    '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                    '\nNon-singletons value\t'+str(int(results_process_C['Bonf_corr'][0]/results_process_C['P'][0]))+
                    '\nCorrection Method\t'+''.join(method_C)+
                    '\nValue\t'+str(User_value_C)+
                    '\n\t\n'+
                    '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                    '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_C, information], axis=0, sort=False)

            a01=[]
            for index, row in results_process_C.iterrows():
                a02=re.findall('[A-Za-z0-9-_]{0,20}',row['entry'])
                for i in a02:
                    a01.append([row['GO'],i])
            component_uniprot=DataFrame(a01)[DataFrame(a01)[1] != ''].reset_index(drop=True).rename(columns={0:'GO',1:'Entry'})

            ## save file with edges for graph
            edges_file_name='edges_Uniprot_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            component_uniprot.to_csv(level_1_uniprot+edges_file_name,index=None)

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx')
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False)
            component_uniprot.to_excel(writer,'Edges Component',index=False)
            writer.save()

            frame=results_process_C[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
            df2=component_uniprot[['Entry']]

            #
            if len(inp_file.columns) >= 2:
                entry_exp=df2.merge(list_input[['Entry','Fold Change']],on='Entry',
                                    how='left').rename(columns={'Fold Change':'Exp'},
                                                       index=str).drop_duplicates().reset_index(drop=True)
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
            nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.csv'
            nodes.drop_duplicates().to_csv(level_1_uniprot+nodes_file_name,index=None)

            ## Open R script from github and run
            r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
            os.makedirs(level_2_uniprot+'Component')
            R_script_enrich = re.sub('qwertyuiop',level_1_uniprot+nodes_file_name,r_script) # name edges file
            R_script_enrich = re.sub('asdfghjkl',level_1_uniprot+edges_file_name,R_script_enrich) # name nodes file
            R_script_enrich = re.sub('zxcvbnm',level_2_uniprot+'Component/',R_script_enrich) # store plots
            R_script_enrich = re.sub('poiuytrewq','GO Annotation (CC)',R_script_enrich) # legend
            R_script_enrich = re.sub('ASPECT','Cellular Component',R_script_enrich)
            f= open(level_1_uniprot+'/Component_Enrichment_Plots.R','w')
            f.write(R_script_enrich)
            f.close()
        else:
            print('\nThere are not enrichment terms for Cellular Component (Uniprot)')
            ## Information Uniprot
            non_annoted=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id[['Entry','GO']],how="left", on='Entry').fillna('N')

            lis=non_annoted[non_annoted.GO == 'N'][['Entry']] 
            if lis.count()[0] == 0:
                lis['ent']=['0']
            else:
                lis['entry']='entry'
                lis['ent']=lis[['Entry']].replace({'$':'; '},regex=True)
                lis=lis.groupby('entry')['ent'].sum().reset_index()

            report = ['\n\t\n'+
                      'NeVOmics\t'+new_folder+
                      '\n\nInput file name\t'+''.join(input_file)+
                      '\nAssociation file name\t'+analysis+
                      "\nTotal number of background\t"+str(background['Entry'].drop_duplicates().count())+
                      "\nTotal number of list\t"+str(list_input['Entry'].drop_duplicates().count())+
                      '\n\nBackground with GO terms\t'+str(back_with_GO['Entry'].drop_duplicates().count())+
                      '\nList input with GO terms\t'+str(list_input_match['Entry'].drop_duplicates().count())+
                      '\nNon-singletons value\t'+str(int(enrich_C['Bonf_corr'][0]/enrich_C['P'][0]))+  # ****
                      '\nCorrection Method\t'+''.join(method_C)+ # ***
                      '\nValue\t'+str(User_value_C)+ # ***
                      '\n\t\n'+
                      '\nProteins with no information in UniprotKB\t'+str(non_annoted[non_annoted.GO == 'N'][['Entry']].count()[0])+
                      '\n'+lis['ent'][0]]

            rep=''.join(report)
            information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
        
            combine=pd.concat([results_process_C, information], axis=0, sort=False) # ***

            ## Excel file with all information
            writer = pd.ExcelWriter(new_folder+'/Enrichment_Uniprot_Analysis_Component_'+''.join(method_C)+'_'+str(User_value_C)+'.xlsx') # ***
            combine.to_excel(writer,'Significant GO terms',index=False)
            enrich_C.to_excel(writer,'Enrichment Results',index=False) # ***
            writer.save()
        
        # °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
        
    else:
        print('\n*****  Cellular Component: No significant terms were found')
        results_process_C=pd.DataFrame(data={'GO':[],'entry':[]})
        proteins_count_C=results_process_C[['entry']].count()[0]
        GO_count_C=results_process_C[['GO']].count()[0]
        User_value_C='0'
        method_C='0'
            


# In[ ]:


## loop for uniprot
database2=['UniprotKB','UniprotKB','UniprotKB']
aspect=['P','F','C']
protein_count=[proteins_count_P,proteins_count_F,proteins_count_C]
user_methods=[''.join(method_P),''.join(method_F),''.join(method_C)]
descriptions=['GO_BP.txt','GO_MF.txt','GO_CC.txt']
user_values=[str(User_value_P),str(User_value_F),str(User_value_C)]
GOs_count=[GO_count_P,GO_count_F,GO_count_C]

summary_Uniprot=[]
for i , j , k , l , m , q in zip(database2, aspect, user_methods, user_values, GOs_count, protein_count):
    summary_Uniprot.append([i ,j , k ,l ,m, q])
    #print(i , j , k , l , m)
summary = summary_GOA + summary_Uniprot
summary = DataFrame(summary).rename(columns={0:'DATABASE',1:'ASPECT',2:'US_METH',3:'US_VAL',4:'GO_TERMS',5:'UNIQ_PROT'})
print('\nENRICHMENT SUMMARY')
print('\n',summary)




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
    root.geometry("850x600")
    root.overrideredirect(1)

    if prompt:
        Label(root, text="\nIf you use NeVOmics in your research, please cite:      ",font=("Arial", 11)).pack()
        Label(root, text="NeVOmics: an enrichment tool for GO and functional\n"+
        "network analysis of data from OMICs technologies",font=("Arial", 11),bg="pale green").pack()
        Label(root, text="\n[ You want to create the network visualizations ]\n"+
              "Select any category\n",font=("Arial", 15)).pack()
    v = IntVar()
    for i, option in enumerate(options):
        Radiobutton(root, text=option,font=("Arial", 13,"bold"),bg="gold2", variable=v, value=i).pack()
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
        "   1    Biological Process       ",
        "   2    Molecular Function       ",
        "   3    Cellular Component       ",
        "   4    All Categories              "
    ]
)
if format(repr(result)) == 'None':
    print('\n!!!!! Graphics not generated !!!!!')
else:
    folders = [level_1_goa,level_1_uniprot]
    # Biological process
    if float(re.findall('[0-9]{1}',format(repr(result)))[0]) == 1:
        # find R scripst process
        plots_selection=[]
        for i in folders:
            plots_selection.append([i+''.join(fnmatch.filter(os.listdir((i)), 'Process*.R'))]) 
        # run R scripts
        for i in plots_selection:
            run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
            run_uni.wait()
    else:
        if float(re.findall('[0-9]{1}',format(repr(result)))[0]) == 2:
            # find R scripst function
            plots_selection=[]
            for i in folders:
                plots_selection.append([i+''.join(fnmatch.filter(os.listdir((i)), 'Function*.R'))])
            # run R scripts
            for i in plots_selection:
                run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
                run_uni.wait()
        else:
            if float(re.findall('[0-9]{1}',format(repr(result)))[0]) == 3:
                # find R scripst component
                plots_selection=[]
                for i in folders:
                    plots_selection.append([i+''.join(fnmatch.filter(os.listdir((i)), 'Component*.R'))])
                 # run R scripts
                for i in plots_selection:
                    run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
                    run_uni.wait()
            else:
                if float(re.findall('[0-9]{1}',format(repr(result)))[0]) == 4:
                    # find R scripst for all categories
                    plots_selection=[]
                    for i in folders:
                        for x in fnmatch.filter(os.listdir((i)), '*.R'):
                            plots_selection.append(i+x)
                    # run R scripts
                    for i in plots_selection:
                        run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
                        run_uni.wait()


# In[ ]:
if os.path.exists('data'): shutil.rmtree('data')
if os.path.exists('GO.py'): os.remove('GO.py')
if os.path.exists('KEGG.py'): os.remove('KEGG.py')
if os.path.exists('HD.py'): os.remove('HD.py')
if os.path.exists("./*.RData"): os.remove("./*.RData")
if os.path.exists("./Rplots.pdf"): os.remove("./Rplots.pdf")

# print total time of analysis
lapso_total = datetime.now() - inicio_total
print('\n'+new_folder+': Analysis Time (hh:mm:ss.ms) {}'.format(lapso_total),'\n')

