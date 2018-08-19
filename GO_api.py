
# coding: utf-8

# In[1]:


## Packages import
import datetime
start = datetime.datetime.now()
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
from time import sleep
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")


# In[ ]:


#from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:95% !important; }</style>"))


# In[2]:


## Control of inputs= organism, id-organism and file name
content_dir_list=os.listdir("./")
content_dir_str =" ".join(content_dir_list)

input_file=input('\n[ Step 1: Submit file (Uniprot IDs) ]\n=====> : ')
if input_file == '':
    print('\n!!!!!!! File not found !!!!!!!')
    input_file=input('\n[ Step 1: Submit file (Uniprot IDs) ]\n=====> : ')
    if input_file == '':
        print('\n!!!!!!! File not found !!!!!!!')
        sys.exit()
else:
    if re.search(input_file,content_dir_str):
        inp_file=pd.read_csv(input_file,sep='\t',header=None)
    else:
        print('\n!!!!!!! File not found !!!!!!!')
        input_file=input('\n[ Step 1: Submit file (Uniprot IDs) ]\n=====> : ')
        if re.search(input_file,content_dir_str):
            inp_file=pd.read_csv(input_file,sep='\t',header=None)
        else:
            print('\n!!!!!!! File not found !!!!!!!\n')
            sys.exit()
    
## Extract first id of protein
inp_file=pd.read_csv(input_file,sep='\t',header=None)
first_entry=inp_file[0].iloc[0]

##############################################################
################           Uniprot         ###################
##############################################################
## exttract id-organism
id_organism = requests.get("https://www.uniprot.org/uniprot/?query="+first_entry+"&sort=score&columns=organism-id,organism&format=tab&limit=1").content.decode()
## Kill process if not found id-organism
if id_organism == '':
    print('\n!!!!!!! ID-Organism not found, Check your identifiers list !!!!!!!')
    sys.exit()
else:
    ## id-organism
    Prefix=DataFrame(re.findall('[0-9]{1,30}',id_organism))[0].iloc[0]
    ## Organism identified
    strain=DataFrame(re.findall('[A-Za-z].*',id_organism))[0].iloc[1]
    print('\n▬ ','Organism identified')
    print('\n▬ ',strain)
    print('..........')
    ## Control of Uniprot directories
    content_dir_list=os.listdir("./")
    content_dir_str =" ".join(content_dir_list)
    if re.search('job_[0-9]{1,3}_GO',content_dir_str):
        plot_dir=DataFrame(re.findall('job_[0-9]{1,3}_GO',content_dir_str)).sort_values(by=[0],ascending=False).reset_index(drop=True)
        plot_dir2=plot_dir[0].str.split('_', expand=True)
        plot_dir3=float(plot_dir2[1].iloc[0])
        plot_dir4=str(int(float(plot_dir3+1)))
        if float(plot_dir4) >= 10:
            os.makedirs('job_'+plot_dir4+'_GO/job_'+plot_dir4+'_Uniprot/'+'job_'+plot_dir4+'_Uniprot_plots')
            dir_name_plots='job_'+plot_dir4+'_GO'
            u='job_0'+plot_dir4
            uu=dir_name_plots+'/'+u+'_Uniprot'
            dir_uniprot=dir_name_plots+'/'+u+'_Uniprot/'+u+'_Uniprot_plots'
        else:
            if os.path.exists('job_'+plot_dir4+'_GO'):
                dir_name_plots='job_0'+plot_dir4+'_GO'
                u='job_0'+plot_dir4
                uu=dir_name_plots+'/'+u+'_Uniprot'
                dir_uniprot=dir_name_plots+'/'+u+'_Uniprot/'+u+'_Uniprot_plots'
            else:
                os.makedirs('job_0'+plot_dir4+'_GO/job_0'+plot_dir4+'_Uniprot/'+'job_0'+plot_dir4+'_Uniprot_plots')
                dir_name_plots='job_0'+plot_dir4+'_GO'
                u='job_0'+plot_dir4
                uu=dir_name_plots+'/'+u+'_Uniprot'
                dir_uniprot=dir_name_plots+'/'+u+'_Uniprot/'+u+'_Uniprot_plots'
    else:
        os.makedirs('job_01_GO/job_01_Uniprot/job_01_Uniprot_plots')
        dir_name_plots='job_01_GO'
        u='job_01'
        uu=dir_name_plots+'/'+u+'_Uniprot'
        dir_uniprot=dir_name_plots+'/'+u+'_Uniprot/'+u+'_Uniprot_plots' 
    
## Create a folder
os.makedirs('data',exist_ok=True)

## Download GeneMerge software from github
output1=urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/GeneMerge1.4.pl', './data/GeneMerge1.4.pl')

## Download all terms from github
output1=urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Terms_IDs.txt', './data/Terms_IDs.txt')

## Open >>>>>>>>>>>>
go_obo_columns=pd.read_csv('./data/Terms_IDs.txt',sep='\t',header=None, names=['GO','Term','Aspect'])

## GO description for all Process category
GO_BP=go_obo_columns[go_obo_columns.Aspect.str.contains("P")==True].reset_index(drop=True)
GO_BP[['GO','Term']].to_csv('./data/GO_BP.txt',sep='\t',index=None,header=None)
## GO description for all Functions category
GO_MF=go_obo_columns[go_obo_columns.Aspect.str.contains("F")==True].reset_index(drop=True)
GO_MF[['GO','Term']].to_csv('./data/GO_MF.txt',sep='\t',index=None,header=None)
## GO description for all Components category
GO_CC=go_obo_columns[go_obo_columns.Aspect.str.contains("C")==True].reset_index(drop=True)
GO_CC[['GO','Term']].to_csv('./data/GO_CC.txt',sep='\t',index=None,header=None)
##########

##############################################################
################           Uniprot         ###################
##############################################################

## Download all acc uniprot and GO IDs
output2=urllib.request.urlretrieve('https://www.uniprot.org/uniprot/?query=organism:'+Prefix+'&format=tab&columns=id,go-id', './data/acc_uniprot_GO_id')

## Edit downloaded file
output3 = Popen("sed -i 's/ //g; s/$/;/g; s/^[A-Z0-9].*\t;//g; /^$/d' ./data/acc_uniprot_GO_id",
                shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
output4 = Popen("sed -i '1d' ./data/acc_uniprot_GO_id",
                shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()

## Open entry and GO-IDs obtained from Uniprot [[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]>>>>>>>>>>>>
acc_uniprot_GO_id=pd.read_csv('./data/acc_uniprot_GO_id',sep='\t',header=None).rename(columns={0:'Entry',1:'GO'})

## Editions in file
acc_uniprot_GO_id=pd.read_csv('./data/acc_uniprot_GO_id',sep='\t',header=None).rename(columns={0:'Entry',1:'GO'})
acc_uniprot_GO_id['Entry']=acc_uniprot_GO_id[['Entry']].replace({'$':';'},regex=True)
acc_uniprot_GO_id['frec']=acc_uniprot_GO_id['GO'].str.split(';').str.len() - 1
acc_uniprot_GO_id['Entry']=acc_uniprot_GO_id.Entry*acc_uniprot_GO_id.frec
## Extract all Entry of pforeins
entrys=acc_uniprot_GO_id['Entry'].str.extractall('([A-Z0-9]{1,20};)').replace({';':''}, regex=True).reset_index(drop=True)
## Extract all GO IDs
gos=acc_uniprot_GO_id['GO'].str.extractall('([A-Z]{2}\W\d{7})').replace({';':''}, regex=True).reset_index(drop=True)
df=[entrys,gos]
## Two columns Entry and GO IDs for mapping [[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]>>>>>>>>>>>>
## ## For background complete (proteome uniprot)
Entry_GO_id=pd.concat(df, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 

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
    
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■

####################  Enrichment Analysis

#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■

##############################################################
################       Uniprot GOA         ###################
##############################################################
## GOA Proteome Sets
## Download Uniprot GOA file
#print('\n▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n               GOA-Uniprot Annotation\n▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n\n               ■ Biological Process\n')

goa_proteomes = urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/20180622/proteome2taxid', './data/goa_proteomes')

## Find Prefix (id-organism)
goa_id_organism = Popen('cut -f2,3 ./data/goa_proteomes | grep '+'^'+Prefix+' | cut -f2',
                shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read().decode()
goa_id_organism=re.sub("\n",'',goa_id_organism)
##
## loop for ID-organism in GOA proteomes
if goa_id_organism == '':
    print('\n!!!!!!! ID-Organism not found in GOA-Uniprot (Complete Annotation). Enrichment will be done only with the Annotation of Uniprot Database !!!!!!!\n')
    if os.path.exists("./data/goa_proteomes"): os.remove("./data/goa_proteomes")
else:
    ## Control of GOA directories
    if os.path.exists(dir_uniprot): os.makedirs(dir_name_plots+'/'+u+'_GOA/'+u+'_GOA_plots')
    oa=dir_name_plots+'/'+u+'_GOA/'
    dir_goa=dir_name_plots+'/'+u+'_GOA/'+u+'_GOA_plots'
    
    #Download Uniprot GOA file
    out2 =urllib.request.urlretrieve("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/20180622/"+goa_id_organism+'.gz', './data/'+goa_id_organism+'.gz')
    out3 = Popen("gzip -d ./data/"+goa_id_organism+".gz",
                 shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    out3 = Popen("sed -i 's/^!.*//g;/^$/d' ./data/"+goa_id_organism,
                shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    orden_columnas_goa=[1,4,6,8,9]
    names=['Entry','GO','Evidence','Aspect','Name']
    ## Background goa file
    goa_file=pd.read_table('./data/'+goa_id_organism,header=None,usecols=orden_columnas_goa,names=names)
    if os.path.exists("./data/goa_proteomes"): os.remove("./data/goa_proteomes")
    
    ## >>>>>>>>>>>> Biological Process
    ##
    ## Preparation of background file
    ## File with background column
    ## User's gackground column
    if len(inp_file.columns) == 3:
        background=list_input[['Background']].rename(columns={'Background':'Entry'}) ## change column name
        back_with_GO=pd.DataFrame.merge(background, goa_file[['Entry','GO']].drop_duplicates(), how="left", on='Entry').dropna()    ## User's background list
        ## Background with GO terms
        back_with_GO[['Entry']].drop_duplicates().to_csv('./data/Background_Proteins.txt',index=None,header=None)
        ## Protein list mapping against "back_with_GO" and then save this list
        list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),back_with_GO,how="left", on='Entry').dropna()
        list_input_match[['Entry']].drop_duplicates().to_csv('./data/Protein_list.txt',header=None,index=False, float_format='%.0f')
        ##
        ## >>>>>>>>>>>> Biological Process    
        ## Build of Process assotiation file
        ## Biological process in list
        entry_go_term=pd.merge(back_with_GO,go_obo_columns,on='GO',how='left').dropna()
        category_P=entry_go_term[entry_go_term.Aspect.str.contains("P")==True]
        ##
        aa=category_P.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
        bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
        pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Process_association.txt',sep='\t',header=None,index=None)
        ##
        ## >>>>>>>>>>>> Molecular Function    
        ## Build of Function assotiation file
        ## Molecular Function in list
        category_F=entry_go_term[entry_go_term.Aspect.str.contains("F")==True]
        ##
        aa=category_F.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
        bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
        pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Function_association.txt',sep='\t',header=None,index=None)
        ##
        ## >>>>>>>>>>>> Cellular Component    
        ## Build of Component assotiation file
        ## Cellular Component in list
        category_C=entry_go_term[entry_go_term.Aspect.str.contains("C")==True]
        ##
        aa=category_C.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
        bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
        pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Component_association.txt',sep='\t',header=None,index=None)
        ##
        ## °°°°°°°°>>>>>>>>>>>> Enrichment for Biological Process ▬▬▬▬▬▬▬▬▬▬▬▬▬▬   
        ##
        orden_columnas=[0,4,5,9,10,11]
        ##
        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        #
        ## Open file for count terms enriched
        pro=pd.read_csv('./data/Process_Enrichment_GOA.csv',usecols=orden_columnas)
        min_val_pro_P=pro['P'].min()
        max_val_pro_P=pro['P'].max()
        ##
        ## °°°°°°°°>>>>>>>>>>>> Enrichment for Molecular Function ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
        output3 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        output4 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        #
        ## Open file for count terms enriched
        fun=pd.read_csv('./data/Function_Enrichment_GOA.csv',usecols=orden_columnas)
        min_val_fun_F=fun['P'].min()
        max_val_fun_F=fun['P'].max()
        ##
        ## °°°°°°°°>>>>>>>>>>>> Enrichment for Cellular Component ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
        output5 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        output6 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        #
        ## Open file for count terms enriched
        com=pd.read_csv('./data/Component_Enrichment_GOA.csv',usecols=orden_columnas)
        min_val_com_C=com['P'].min()
        max_val_com_C=com['P'].max()
        ###

    ## File without background column
    else:
        ## Background complete (proteome)
        goa_file[['Entry']].drop_duplicates().to_csv('./data/Background_Proteins.txt',index=None,header=None)
        ## Mapping bwteen gene input and background
        list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),goa_file, how="left", on='Entry').dropna()
        ## Save Proteins list found in background database
        list_input_match[['Entry']].drop_duplicates().to_csv('./data/Protein_list.txt',header=None,index=False, float_format='%.0f')    ## Build of kegg assotiation file
        ##
        ## >>>>>>>>>>>> Biological Process    
        ## Build of Process assotiation file
        ## Biological process in list
        entry_go_term=pd.merge(goa_file[['Entry','GO']],go_obo_columns,on='GO',how='left').dropna()
        category_P=entry_go_term[entry_go_term.Aspect.str.contains("P")==True]
        #
        aa=category_P.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
        bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
        pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Process_association.txt',sep='\t',header=None,index=None)
        ##
        ## >>>>>>>>>>>> Molecular Function    
        ## Build of Function assotiation file
        ## Molecular Function in list
        ##
        entry_go_term=pd.merge(goa_file[['Entry','GO']],go_obo_columns,on='GO',how='left').dropna()
        category_F=entry_go_term[entry_go_term.Aspect.str.contains("F")==True]
        #
        aa=category_F.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
        bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
        pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Function_association.txt',sep='\t',header=None,index=None)
        ##
        ## >>>>>>>>>>>> Cellular Component    
        ## Build of Component assotiation file
        ## Cellular Component in list
        ##
        entry_go_term=pd.merge(goa_file[['Entry','GO']],go_obo_columns,on='GO',how='left').dropna()
        category_C=entry_go_term[entry_go_term.Aspect.str.contains("C")==True]
        #
        aa=category_C.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
        bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
        pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Component_association.txt',sep='\t',header=None,index=None)
        ##
        ## °°°°°°°°>>>>>>>>>>>> Enrichment for Biological Process ▬▬▬▬▬▬▬▬▬▬▬▬▬▬   
        ##
        orden_columnas=[0,4,5,9,10,11]
        ##
        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        #
        ## Open file for count terms enriched
        pro=pd.read_csv('./data/Process_Enrichment_GOA.csv',usecols=orden_columnas)
        min_val_pro_P=pro['P'].min()
        max_val_pro_P=pro['P'].max()
        ##
        ## °°°°°°°°>>>>>>>>>>>> Enrichment for Molecular Function ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
        output3 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        output4 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        #
        ## Open file for count terms enriched
        fun=pd.read_csv('./data/Function_Enrichment_GOA.csv',usecols=orden_columnas)
        min_val_fun_F=fun['P'].min()
        max_val_fun_F=fun['P'].max()
        ##
        ## °°°°°°°°>>>>>>>>>>>> Enrichment for Cellular Component ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
        output5 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        output6 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
        #
        ## Open file for count terms enriched
        com=pd.read_csv('./data/Component_Enrichment_GOA.csv',usecols=orden_columnas)
        min_val_com_C=com['P'].min()
        max_val_com_C=com['P'].max()
        ###
        ################################# Process
        
    if pro[(pro.P < 0.05)]['P'].count() >= 1:
        print('▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n                 Uniprot-GOA Annotation\n▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n\n              ■  Biological Process\n')
        print('======== ',pro[(pro.P < 0.05)]['P'].count(),'Biological Processes were found with P-value < 0.05')
        print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
        ## Choice of method
        Bonferroni='Bonferroni'
        FDR='FDR'
        #
        User_method=input('\n[ Step 2: Choose a Method (e.g., FDR / Bonferroni) ]\n=====> : ')
        if User_method == Bonferroni:
            print('')
        else:
            if User_method == FDR:
                print('')
            else:       
                #print('\n!!!!! Try again !!!!!')
                User_method=input('\n[ Step 2: Choose a method (e.g., FDR / Bonferroni) ]\n=====> : ')
                if User_method == Bonferroni:
                    print('')
                else:
                    if User_method == FDR:
                        print('')
                    else:
                        print('\n!!!!! Incorrect Method !!!!!')
                        sys.exit()
        ##  cut-off for Bonferroni                
        if User_method == Bonferroni:
            print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
            Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',Bon_value)
            if match:
                print('\nEnter a numeric value\n')
                print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',Bon_value)
                if match:
                    print('\nIt is not a numerical value')
                    sys.exit()
                else:
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                        null=''
                    else:
                        print('\nIncorrect value')
                        print('\n===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                        Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
                        Bon_cut_off=float(Bon_value)
                        file_name_value=Bon_value
                        if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                            null=''
                        else:
                            print('\nIncorrect value')
                            sys.exit()
            else:
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                    null=''
                else:
                    print('\nIncorrect value')
                    print('\n===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                    Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                        null=''
                    else:
                        print('\nIncorrect value')
                        sys.exit()
            
        ## cut-off for FDR            
        else:
            if User_method == FDR:
                print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                FDR_cut_off=float(FDR_value)*100
                file_name_value=FDR_value
                match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                if match:
                    print('\nEnter a numeric value\n')
                    print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                    FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    #FDR_cut_off=float(FDR_value)*100
                    #file_name_value=FDR_value
                    match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                    if match:
                        print('\nIt is not a numerical value')
                        sys.exit()
                    else:
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        else:
                            print('\nIncorrect value')
                            print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                            FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                            FDR_cut_off=float(FDR_value)*100
                            file_name_value=FDR_value
                            if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                                output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                                output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                            else:
                                print('\nIncorrect value')
                                sys.exit()
                else:
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    else:
                        print('\nIncorrect value')
                        print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                        FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                        else:
                            print('\nIncorrect value')
                            sys.exit()
        ## Open raw data of enrichment GOA
        goaP=pd.read_csv('./data/Process_Enrichment_GOA.csv',usecols=orden_columnas)
        if User_method == Bonferroni:
            goaP_user_cut_off=goaP[(goaP.adj_pval <= Bon_cut_off) == True]
            goaP_method=User_method
            goaP_cutoff=Bon_cut_off
        else:
            if User_method == FDR:
                filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
                goaP_user_cut_off=goaP[goaP[filter_fdr_T].str.contains("T") == True]
                goaP_method=User_method
                goaP_cutoff=FDR_cut_off/100
            else: 
                null=''
        ## Enrichment pathways after cut-off
        goaP_user_cut_off['Freq']=goaP_user_cut_off['Entry'].str.split().str.len()
        
    else:
        print('\n        ■■■■■  Biological Process\n')
        print('No significant terms were found for Biological Processes')
        goaP_user_cut_off={'GO':[]}
        goaP_user_cut_off=pd.DataFrame(data=goaP_user_cut_off)
        goaP_method='NA'
        goaP_cutoff='NA'

#################################### Function

    if fun[(fun.P < 0.05)]['P'].count() >= 1:
        print('\n              ■  Molecular Function\n')
        print('======== ',fun[(fun.P < 0.05)]['P'].count(),'Molecular Functions were found with P-value < 0.05')
        print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
        ## Choice of method
        Bonferroni='Bonferroni'
        FDR='FDR'
        #
        User_method=input('\n[ Step 4: Choose a Method (e.g., FDR / Bonferroni) ]\n=====> : ')
        if User_method == Bonferroni:
            print('')
        else:
            if User_method == FDR:
                print('')
            else:       
                #print('\n!!!!! Try again !!!!!')
                User_method=input('\n[ Step 4: Choose a method (e.g., FDR / Bonferroni) ]\n=====> : ')
                if User_method == Bonferroni:
                    print('')
                else:
                    if User_method == FDR:
                        print('')
                    else:
                        print('\n!!!!! Incorrect Method !!!!!')
                        sys.exit()

        ##  cut-off for Bonferroni                
        if User_method == Bonferroni:
            print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
            Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',Bon_value)
            if match:
                print('\nEnter a numeric value\n')
                print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',Bon_value)
                if match:
                    print('\nIt is not a numerical value')
                    sys.exit()
                else:
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                        null=''
                    else:
                        print('\nIncorrect value')
                        print('\n===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                        Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
                        Bon_cut_off=float(Bon_value)
                        file_name_value=Bon_value
                        if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                            null=''
                        else:
                            print('\nIncorrect value')
                        sys.exit()
            else:
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                    null=''
                else:
                    print('\nIncorrect value')
                    print('\n===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                    Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                        null=''
                    else:
                        print('\nIncorrect value')
                        sys.exit()
            
        ## cut-off for FDR            
        else:
            if User_method == FDR:
                print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                if match:
                    print('\nEnter a numeric value\n')
                    print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                    FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                    if match:
                        print('\nIt is not a numerical value')
                        sys.exit()
                    else:
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        else:
                            print('\nIncorrect value')
                            print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                            FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                            FDR_cut_off=float(FDR_value)*100
                            file_name_value=FDR_value
                            if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                                output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                                output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                            else:
                                print('\nIncorrect value')
                                sys.exit()
                else:
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    else:
                        print('\nIncorrect value')
                        print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                        FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                        else:
                            print('\nIncorrect value')
                            sys.exit()
        ## Open raw data of enrichment GOA
        goaF=pd.read_csv('./data/Function_Enrichment_GOA.csv',usecols=orden_columnas)
        if User_method == Bonferroni:
            goaF_user_cut_off=goaF[(goaF.adj_pval <= Bon_cut_off) == True]
            goaF_method=User_method
            goaF_cutoff=Bon_cut_off 
        else:
            if User_method == FDR:
                filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
                goaF_user_cut_off=goaF[goaF[filter_fdr_T].str.contains("T") == True]
                goaF_method=User_method
                goaF_cutoff=FDR_cut_off/100
            else: 
                null=''
        ## Enrichment pathways after cut-off
        goaF_user_cut_off['Freq']=goaF_user_cut_off['Entry'].str.split().str.len()
       
    else:
        print('\n        ■■■■■  Molecular Function\n')
        print('No significant terms were found for Molecular Functions')
        goaF_user_cut_off={'GO':[]}
        goaF_user_cut_off=pd.DataFrame(data=goaF_user_cut_off)
        goaF_method='NA'
        goaF_cutoff='NA'

#################################### Component

    if com[(com.P < 0.05)]['P'].count() >= 1:
        print('\n              ■  Cellular Component\n')
        print('======== ',com[(pro.P < 0.05)]['P'].count(),'Cellular Components were found with P-value < 0.05')
        print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
        ## Choice of method
        Bonferroni='Bonferroni'
        FDR='FDR'
        #
        User_method=input('\n[ Step 6: Choose a Method (e.g., FDR / Bonferroni) ]\n=====> : ')
        if User_method == Bonferroni:
            print('')
        else:
            if User_method == FDR:
                print('')
            else:       
                #print('\n!!!!! Try again !!!!!')
                User_method=input('\n[ Step 6: Choose a method (e.g., FDR / Bonferroni) ]\n=====> : ')
                if User_method == Bonferroni:
                    print('')
                else:
                    if User_method == FDR:
                        print('')
                    else:
                        print('\n!!!!! Incorrect Method !!!!!')
                        sys.exit()

        ##  cut-off for Bonferroni                
        if User_method == Bonferroni:
            print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
            Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',Bon_value)
            if match:
                print('\nEnter a numeric value\n')
                print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',Bon_value)
                if match:
                    print('\nIt is not a numerical value')
                    sys.exit()
                else:
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                        null=''
                    else:
                        print('\nIncorrect value')
                        print('\n===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                        Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
                        Bon_cut_off=float(Bon_value)
                        file_name_value=Bon_value
                        if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                            null=''
                        else:
                            print('\nIncorrect value')
                            sys.exit()
            else:
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                    null=''
                else:
                    print('\nIncorrect value')
                    print('\n===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                    Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                        null=''
                    else:
                        print('\nIncorrect value')
                        sys.exit()
            
        ## cut-off for FDR            
        else:
            if User_method == FDR:
                print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                if match:
                    print('\nEnter a numeric value\n')
                    print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                    FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                    if match:
                        print('\nIt is not a numerical value')
                        sys.exit()
                    else:
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        else:
                            print('\nIncorrect value')
                            print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                            FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                            FDR_cut_off=float(FDR_value)*100
                            file_name_value=FDR_value
                            if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                                output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                                output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                            else:
                                print('\nIncorrect value')
                                sys.exit()
                else:
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    else:
                        print('\nIncorrect value')
                        print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                        FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_GOA.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_GOA.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                        else:
                            print('\nIncorrect value')
                            sys.exit()
        ## Open raw data of enrichment GOA
        goaC=pd.read_csv('./data/Component_Enrichment_GOA.csv',usecols=orden_columnas)
        if User_method == Bonferroni:
            goaC_user_cut_off=goaC[(goaC.adj_pval <= Bon_cut_off) == True]
            goaC_method=User_method
            goaC_cutoff=Bon_cut_off 
        else:
            if User_method == FDR:
                filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
                goaC_user_cut_off=goaC[goaC[filter_fdr_T].str.contains("T") == True]
                goaC_method=User_method
                goaC_cutoff=FDR_cut_off/100
            else: 
                null=''
        ## Enrichment pathways after cut-off
        goaC_user_cut_off['Freq']=goaC_user_cut_off['Entry'].str.split().str.len()
        
    else:
        print('\n        ■■■■■  Cellular Component\n')
        print('No significant terms were found for Cellular Components')
        goaC_user_cut_off={'GO':[]}
        goaC_user_cut_off=pd.DataFrame(data=goaC_user_cut_off)
        goaC_method='NA'
        goaC_cutoff='NA'
    
    ## Open raw data of enrichment GOA
    #goaC=pd.read_csv('./data/Component_Enrichment_GOA.csv',usecols=orden_columnas)
    #if User_method == Bonferroni:
        #goaC_user_cut_off=goaC[(goaC.adj_pval <= Bon_cut_off) == True]
        #goaC_method=User_method
        #goaC_cutoff=Bon_cut_off 
    #else:
        #if User_method == FDR:
            #filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
            #goaC_user_cut_off=goaC[goaC[filter_fdr_T].str.contains("T") == True]
            #goaC_method=User_method
            #goaC_cutoff=FDR_cut_off/100
        #else: 
            #null=''
    ## Enrichment pathways after cut-off
    #goaC_user_cut_off['Freq']=goaC_user_cut_off['Entry'].str.split().str.len()

#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
#■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
##############################################################
################           Uniprot         ###################
##############################################################

## >>>>>>>>>>>> Biological Process
##
## Preparation of background file
## File with background column
## User's gackground column
if len(inp_file.columns) == 3:
    background=list_input[['Background']].rename(columns={'Background':'Entry'}) ## change column name
    back_with_GO=pd.DataFrame.merge(background, Entry_GO_id, how="left", on='Entry').dropna()    ## User's background list
    ## Background with GO terms
    back_with_GO[['Entry']].drop_duplicates().to_csv('./data/Background_Proteins.txt',index=None,header=None)
    ## Protein list mapping against "back_with_GO" and then save this list
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),back_with_GO,how="left", on='Entry').dropna()
    list_input_match[['Entry']].drop_duplicates().to_csv('./data/Protein_list.txt',header=None,index=False, float_format='%.0f')
    ##
    ## >>>>>>>>>>>> Biological Process    
    ## Build of Process assotiation file
    ## Biological process in list
    entry_go_term=pd.merge(back_with_GO,go_obo_columns,on='GO',how='left').dropna()
    category_P=entry_go_term[entry_go_term.Aspect.str.contains("P")==True]
    ##
    aa=category_P.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
    bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Process_association.txt',sep='\t',header=None,index=None)
    ##
    ## >>>>>>>>>>>> Molecular Function    
    ## Build of Function assotiation file
    ## Molecular Function in list
    category_F=entry_go_term[entry_go_term.Aspect.str.contains("F")==True]
    ##
    aa=category_F.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
    bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Function_association.txt',sep='\t',header=None,index=None)
    ##
    ## >>>>>>>>>>>> Cellular Component    
    ## Build of Component assotiation file
    ## Cellular Component in list
    category_C=entry_go_term[entry_go_term.Aspect.str.contains("C")==True]
    ##
    aa=category_C.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
    bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Component_association.txt',sep='\t',header=None,index=None)
    ##
    ## °°°°°°°°>>>>>>>>>>>> Enrichment for Biological Process ▬▬▬▬▬▬▬▬▬▬▬▬▬▬   
    ##
    orden_columnas=[0,4,5,9,10,11]
    ##
    output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #
    ## Open file for count terms enriched
    pro=pd.read_csv('./data/Process_Enrichment_Uniprot.csv',usecols=orden_columnas)
    min_val_pro_P=pro['P'].min()
    max_val_pro_P=pro['P'].max()
    ##
    ## °°°°°°°°>>>>>>>>>>>> Enrichment for Molecular Function ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
    output3 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    output4 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #
    ## Open file for count terms enriched
    fun=pd.read_csv('./data/Function_Enrichment_Uniprot.csv',usecols=orden_columnas)
    min_val_fun_F=fun['P'].min()
    max_val_fun_F=fun['P'].max()
    ##
    ## °°°°°°°°>>>>>>>>>>>> Enrichment for Cellular Component ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
    output5 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    output6 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #
    ## Open file for count terms enriched
    com=pd.read_csv('./data/Component_Enrichment_Uniprot.csv',usecols=orden_columnas)
    min_val_com_C=com['P'].min()
    max_val_com_C=com['P'].max()
    ###

## File without background column
else:
    ## Background complete (proteome)
    Entry_GO_id[['Entry']].drop_duplicates().to_csv('./data/Background_Proteins.txt',index=None,header=None)
    ## Mapping bwteen gene input and background
    list_input_match=pd.DataFrame.merge(list_input[['Entry']].dropna(),Entry_GO_id, how="left", on='Entry').dropna()
    ## Save Proteins list found in background database
    list_input_match[['Entry']].drop_duplicates().to_csv('./data/Protein_list.txt',header=None,index=False, float_format='%.0f')    ## Build of kegg assotiation file
    ##
    ## >>>>>>>>>>>> Biological Process    
    ## Build of Process assotiation file
    ## Biological process in list
    entry_go_term=pd.merge(Entry_GO_id,go_obo_columns,on='GO',how='left').dropna()
    category_P=entry_go_term[entry_go_term.Aspect.str.contains("P")==True]
    #
    aa=category_P.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
    bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Process_association.txt',sep='\t',header=None,index=None)
    ##
    ## >>>>>>>>>>>> Molecular Function    
    ## Build of Function assotiation file
    ## Molecular Function in list
    ##
    entry_go_term=pd.merge(Entry_GO_id,go_obo_columns,on='GO',how='left').dropna()
    category_F=entry_go_term[entry_go_term.Aspect.str.contains("F")==True]
    #
    aa=category_F.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
    bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Function_association.txt',sep='\t',header=None,index=None)
    ##
    ## >>>>>>>>>>>> Cellular Component    
    ## Build of Component assotiation file
    ## Cellular Component in list
    ##
    entry_go_term=pd.merge(Entry_GO_id,go_obo_columns,on='GO',how='left').dropna()
    category_C=entry_go_term[entry_go_term.Aspect.str.contains("C")==True]
    #
    aa=category_C.pivot_table(values='GO',index=['Entry'],aggfunc=sum).reset_index()
    bb=aa[['GO']].replace({'GO':';GO'},regex=True).replace({'$':';'},regex=True).replace({'^;':''},regex=True)
    pd.concat([aa[['Entry']],bb[['GO']]],axis=1,join='outer').to_csv('./data/Component_association.txt',sep='\t',header=None,index=None)
    ##
    ## °°°°°°°°>>>>>>>>>>>> Enrichment for Biological Process ▬▬▬▬▬▬▬▬▬▬▬▬▬▬   
    ##
    orden_columnas=[0,4,5,9,10,11]
    ##
    output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #
    ## Open file for count terms enriched
    pro=pd.read_csv('./data/Process_Enrichment_Uniprot.csv',usecols=orden_columnas)
    min_val_pro_P=pro['P'].min()
    max_val_pro_P=pro['P'].max()
    ##
    ## °°°°°°°°>>>>>>>>>>>> Enrichment for Molecular Function ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
    output3 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    output4 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #
    ## Open file for count terms enriched
    fun=pd.read_csv('./data/Function_Enrichment_Uniprot.csv',usecols=orden_columnas)
    min_val_fun_F=fun['P'].min()
    max_val_fun_F=fun['P'].max()
    ##
    ## °°°°°°°°>>>>>>>>>>>> Enrichment for Cellular Component ▬▬▬▬▬▬▬▬▬▬▬▬▬▬
    output5 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    output6 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
    #
    ## Open file for count terms enriched
    com=pd.read_csv('./data/Component_Enrichment_Uniprot.csv',usecols=orden_columnas)
    min_val_com_C=com['P'].min()
    max_val_com_C=com['P'].max()
    ###

################################# Process

if pro[(pro.P < 0.05)]['P'].count() >= 1:
    print('\n▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n                 Uniprot Annotation\n▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬\n\n               ■  Biological Process\n')
    print('======== ',pro[(pro.P < 0.05)]['P'].count(),'Biological Processes were found with P-value < 0.05')
    print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
    ## Choice of method
    Bonferroni='Bonferroni'
    FDR='FDR'
    #
    User_method=input('\n[ Step 2: Choose a Method (e.g., FDR / Bonferroni) ]\n=====> : ')
    if User_method == Bonferroni:
        print('')
    else:
        if User_method == FDR:
            print('')
        else:       
            #print('\n!!!!! Try again !!!!!')
            User_method=input('\n[ Step 2: Choose a method (e.g., FDR / Bonferroni) ]\n=====> : ')
            if User_method == Bonferroni:
                print('')
            else:
                if User_method == FDR:
                    print('')
                else:
                    print('\n!!!!! Incorrect Method !!!!!')
                    sys.exit()

    ##  cut-off for Bonferroni                
    if User_method == Bonferroni:
        print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
        Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
        match = re.search(r'[A-Za-z]{1,10}',Bon_value)
        if match:
            print('\nEnter a numeric value\n')
            print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
            Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',Bon_value)
            if match:
                print('\nIt is not a numerical value')
                sys.exit()
            else:
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                    null=''
                else:
                    print('\nIncorrect value')
                    print('\n===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                    Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                        null=''
                    else:
                        print('\nIncorrect value')
                        sys.exit()
        else:
            Bon_cut_off=float(Bon_value)
            file_name_value=Bon_value
            if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                null=''
            else:
                print('\nIncorrect value')
                print('\n===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                Bon_value=input('\n[ Step 3: Choose a Value (e.g., 0.1) ]\n=====> : ')
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_pro_P <= Bon_cut_off <= max_val_pro_P:
                    null=''
                else:
                    print('\nIncorrect value')
                    sys.exit()
            
    ## cut-off for FDR            
    else:
        if User_method == FDR:
            print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
            FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',FDR_value)
            if match:
                    print('\nEnter a numeric value\n')
                    print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                    FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                    if match:
                        print('\nIt is not a numerical value')
                        sys.exit()
                    else:
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        else:
                            print('\nIncorrect value')
                            print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                            FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                            FDR_cut_off=float(FDR_value)*100
                            file_name_value=FDR_value
                            if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                                output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                                output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                            else:
                                print('\nIncorrect value')
                                sys.exit()
            else:
                FDR_cut_off=float(FDR_value)*100
                file_name_value=FDR_value
                if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                    output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                else:
                    print('\nIncorrect value')
                    print('===== Must be in this range: ',min_val_pro_P,' - ',max_val_pro_P)
                    FDR_value=input('\n[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_pro_P <= FDR_cut_off/100 <= max_val_pro_P:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Process_association.txt ./data/GO_BP.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Process_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Process_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                    else:
                        print('\nIncorrect value')
                        sys.exit()
    ## Open raw data of enrichment GOA
    uniprotP=pd.read_csv('./data/Process_Enrichment_Uniprot.csv',usecols=orden_columnas)
    if User_method == Bonferroni:
        uniP_user_cut_off=uniprotP[(uniprotP.adj_pval <= Bon_cut_off) == True]
        uniP_method=User_method
        uniP_cutoff=Bon_cut_off 
    else:
        if User_method == FDR:
            filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
            uniP_user_cut_off=uniprotP[uniprotP[filter_fdr_T].str.contains("T") == True]
            uniP_method=User_method
            uniP_cutoff=FDR_cut_off/100
        else: 
            null=''
    ## Enrichment pathways after cut-off
    uniP_user_cut_off['Freq']=uniP_user_cut_off['Entry'].str.split().str.len()
    
else:
    print('\n        ■■■■■  Biological Process\n')
    print('No significant terms were found for Biological Processes')
    uniP_user_cut_off={'GO':[]}
    uniP_user_cut_off=pd.DataFrame(data=uniP_user_cut_off)
    uniP_method='NA'
    uniP_cutoff='NA'

#################################### Function

if fun[(fun.P < 0.05)]['P'].count() >= 1:
    print('\n              ■  Molecular Function\n')
    print('======== ',fun[(fun.P < 0.05)]['P'].count(),'Molecular Functions were found with P-value < 0.05')
    print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
    ## Choice of method
    Bonferroni='Bonferroni'
    FDR='FDR'
    #
    User_method=input('\n[ Step 4: Choose a Method (e.g., FDR / Bonferroni) ]\n=====> : ')
    if User_method == Bonferroni:
        print('')
    else:
        if User_method == FDR:
            print('')
        else:       
            #print('\n!!!!! Try again !!!!!')
            User_method=input('\n[ Step 4: Choose a method (e.g., FDR / Bonferroni) ]\n=====> : ')
            if User_method == Bonferroni:
                print('')
            else:
                if User_method == FDR:
                    print('')
                else:
                    print('\n!!!!! Incorrect Method !!!!!')
                    sys.exit()

    ##  cut-off for Bonferroni                
    if User_method == Bonferroni:
        print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
        Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
        match = re.search(r'[A-Za-z]{1,10}',Bon_value)
        if match:
            print('\nEnter a numeric value\n')
            print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
            Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',Bon_value)
            if match:
                print('\nIt is not a numerical value')
                sys.exit()
            else:
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                    null=''
                else:
                    print('\nIncorrect value')
                    print('\n===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                    Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                        null=''
                    else:
                        print('\nIncorrect value')
                        sys.exit()
        else:
            Bon_cut_off=float(Bon_value)
            file_name_value=Bon_value
            if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                null=''
            else:
                print('\nIncorrect value')
                print('\n===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                Bon_value=input('\n[ Step 5: Choose a Value (e.g., 0.1) ]\n=====> : ')
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_fun_F <= Bon_cut_off <= max_val_fun_F:
                    null=''
                else:
                    print('\nIncorrect value')
                    sys.exit()
            
    ## cut-off for FDR            
    else:
        if User_method == FDR:
            print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
            FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',FDR_value)
            if match:
                print('\nEnter a numeric value\n')
                print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                if match:
                    print('\nIt is not a numerical value')
                    sys.exit()
                else:
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    else:
                        print('\nIncorrect value')
                        print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                        FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                        else:
                            print('\nIncorrect value')
                            sys.exit()
            else:
                FDR_cut_off=float(FDR_value)*100
                file_name_value=FDR_value
                if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                    output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                else:
                    print('\nIncorrect value')
                    print('===== Must be in this range: ',min_val_fun_F,' - ',max_val_fun_F)
                    FDR_value=input('\n[ Step 5: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_fun_F <= FDR_cut_off/100 <= max_val_fun_F:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Function_association.txt ./data/GO_MF.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Function_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Function_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                    else:
                        print('\nIncorrect value')
                        sys.exit()
    ## Open raw data of enrichment GOA
    uniprotF=pd.read_csv('./data/Function_Enrichment_Uniprot.csv',usecols=orden_columnas)
    if User_method == Bonferroni:
        uniF_user_cut_off=uniprotF[(uniprotF.adj_pval <= Bon_cut_off) == True]
        uniF_method=User_method
        uniF_cutoff=Bon_cut_off 
    else:
        if User_method == FDR:
            filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
            uniF_user_cut_off=uniprotF[uniprotF[filter_fdr_T].str.contains("T") == True]
            uniF_method=User_method
            uniF_cutoff=FDR_cut_off/100
        else: 
            null=''
    ## Enrichment pathways after cut-off
    uniF_user_cut_off['Freq']=uniF_user_cut_off['Entry'].str.split().str.len()
    
else:
    print('\n        ■■■■■  Molecular Function\n')
    print('No significant terms were found for Molecular Functions')
    uniF_user_cut_off={'GO':[]}
    uniF_user_cut_off=pd.DataFrame(data=uniF_user_cut_off)
    uniF_method='NA'
    uniF_cutoff='NA'
    
#################################### Component

if com[(com.P < 0.05)]['P'].count() >= 1:
    print('\n              ■  Cellular Component\n')
    print('======== ',com[(pro.P < 0.05)]['P'].count(),'Cellular Components were found with P-value < 0.05')
    print('\nChoice a correction method for P-value (e.g., FDR .05 / Bonferroni 0.1)')
    ## Choice of method
    Bonferroni='Bonferroni'
    FDR='FDR'
    #
    User_method=input('\n[ Step 6: Choose a Method (e.g., FDR / Bonferroni) ]\n=====> : ')
    if User_method == Bonferroni:
        print('')
    else:
        if User_method == FDR:
            print('')
        else:       
            #print('\n!!!!! Try again !!!!!')
            User_method=input('\n[ Step 6: Choose a method (e.g., FDR / Bonferroni) ]\n=====> : ')
            if User_method == Bonferroni:
                print('')
            else:
                if User_method == FDR:
                    print('')
                else:
                    print('\n!!!!! Incorrect Method !!!!!')
                    sys.exit()

    ##  cut-off for Bonferroni                
    if User_method == Bonferroni:
        print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
        Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
        match = re.search(r'[A-Za-z]{1,10}',Bon_value)
        if match:
            print('\nEnter a numeric value\n')
            print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
            Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',Bon_value)
            if match:
                print('\nIt is not a numerical value')
                sys.exit()
            else:
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                    null=''
                else:
                    print('\nIncorrect value')
                    print('\n===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                    Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
                    Bon_cut_off=float(Bon_value)
                    file_name_value=Bon_value
                    if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                        null=''
                    else:
                        print('\nIncorrect value')
                        sys.exit()
        else:
            Bon_cut_off=float(Bon_value)
            file_name_value=Bon_value
            if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                null=''
            else:
                print('\nIncorrect value')
                print('\n===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                Bon_value=input('\n[ Step 7: Choose a Value (e.g., 0.1) ]\n=====> : ')
                Bon_cut_off=float(Bon_value)
                file_name_value=Bon_value
                if min_val_com_C <= Bon_cut_off <= max_val_com_C:
                    null=''
                else:
                    print('\nIncorrect value')
                    sys.exit()
            
    ## cut-off for FDR            
    else:
        if User_method == FDR:
            print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
            FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
            match = re.search(r'[A-Za-z]{1,10}',FDR_value)
            if match:
                print('\nEnter a numeric value\n')
                print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                match = re.search(r'[A-Za-z]{1,10}',FDR_value)
                if match:
                    print('\nIt is not a numerical value')
                    sys.exit()
                else:
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    else:
                        print('\nIncorrect value')
                        print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                        FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                        FDR_cut_off=float(FDR_value)*100
                        file_name_value=FDR_value
                        if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                            output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                            output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                        else:
                            print('\nIncorrect value')
                            sys.exit()
            else:
                FDR_cut_off=float(FDR_value)*100
                file_name_value=FDR_value
                if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                    output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                    output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                else:
                    print('\nIncorrect value')
                    print('===== Must be in this range: ',min_val_com_C,' - ',max_val_com_C)
                    FDR_value=input('\n[ Step 7: Choose a Value (e.g., 0.05) ]\n=====> : ')
                    FDR_cut_off=float(FDR_value)*100
                    file_name_value=FDR_value
                    if min_val_com_C <= FDR_cut_off/100 <= max_val_com_C:
                        output1 = Popen("perl ./data/GeneMerge1.4.pl ./data/Component_association.txt ./data/GO_CC.txt ./data/Background_Proteins.txt ./data/Protein_list.txt ./data/Component_Enrichment_Uniprot.csv "+str(FDR_cut_off)+"%", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read()
                        output2 = Popen("sed -i 's/\t/"'"'","'"'"/g; s/, "'"'"/"'"'"/g' ./data/Component_Enrichment_Uniprot.csv", shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).stdout.read() 
                    else:
                        print('\nIncorrect value')
                        sys.exit()
    ## Open raw data of enrichment GOA
    uniprotC=pd.read_csv('./data/Component_Enrichment_Uniprot.csv',usecols=orden_columnas)
    if User_method == Bonferroni:
        uniC_user_cut_off=uniprotC[(uniprotC.adj_pval <= Bon_cut_off) == True]
        uniC_method=User_method
        uniC_cutoff=Bon_cut_off 
    else:
        if User_method == FDR:
            filter_fdr_T= 'FDR_'+str(FDR_cut_off)+'_perc'
            uniC_user_cut_off=uniprotC[uniprotC[filter_fdr_T].str.contains("T") == True]
            uniC_method=User_method
            uniC_cutoff=FDR_cut_off/100
        else: 
            null=''
    ## Enrichment pathways after cut-off
    uniC_user_cut_off['Freq']=uniC_user_cut_off['Entry'].str.split().str.len()
    
else:
    print('\n        ■■■■■  Cellular Component\n')
    print('No significant terms were found for Cellular Components')
    uniC_user_cut_off={'GO':[]}
    uniC_user_cut_off=pd.DataFrame(data=uniC_user_cut_off)
    uniC_method='NA'
    uniC_cutoff='NA'
    
###################
print('\n▬ Summary')
results = {'Annotaion':['Uniprot-GOA','Uniprot-GOA','Uniprot-GOA','Uniprot','Uniprot','Uniprot'],
           'Aspect':['Process','Function','Component','Process','Function','Component'],
           'GO': [goaP_user_cut_off['GO'].count(),goaF_user_cut_off['GO'].count(),goaC_user_cut_off['GO'].count(),
                 uniP_user_cut_off['GO'].count(),uniF_user_cut_off['GO'].count(),uniC_user_cut_off['GO'].count()],
           'Method':[goaP_method,goaF_method,goaC_method,uniP_method,uniF_method,uniC_method],
          'cut-off':[goaP_cutoff,goaF_cutoff,goaC_cutoff,uniP_cutoff,uniF_cutoff,uniC_cutoff]}
table= pd.DataFrame(data=results)
print(table)
if os.path.exists("GeneMerge1.4.pl"): os.remove("GeneMerge1.4.pl")
if os.path.exists('data'): shutil.rmtree('data')
#print('\n',table,'\n')
##
print('\n[ Step 8: Generating graphics ]\n..........')
##############################################################
#####################       GOA         ######################
##############################################################
## Data editing goa process
if goaP_user_cut_off['GO'].count() >= 1:
    enriched_genes=DataFrame(goaP_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
    equalize=str(goaP_user_cut_off.GO*goaP_user_cut_off.Freq)
    id_pathways=re.findall('GO........',equalize)
    enriched_pathways=DataFrame(id_pathways)
    frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
    gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 
    Gene_Entry_Pathway_Enriched=gene_kegg.merge(go_obo_columns,on='GO',how='left')
    ## File name for edges goa process
    edges_file_name='edges_GOA_Enrichment_Analysis_Process_'+User_method+'_'+file_name_value+'.csv'
    Gene_Entry_Pathway_Enriched[['GO','Entry']].to_csv(edges_file_name,index=None)
    xlsx='GOA_Enrichment_Analysis_Process_'+User_method+'_'+file_name_value+'.xlsx'
    Gene_Entry_Pathway_Enriched.to_csv(xlsx,index=None)
    ## Preparing data frame for nodes file
    frame=Gene_Entry_Pathway_Enriched[['GO']].rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    frame['Freq']=goaP_user_cut_off['Freq']
    frame['num']=frame.reset_index().index + 1
    df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Entry']],index=None)
    df2=pd.read_csv(StringIO(df1))
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
    #
    entry_exp['Freq']=float(goaP_user_cut_off['Freq'].min()*.5)
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
    entry_pathways=Gene_Entry_Pathway_Enriched[['GO','Term']].drop_duplicates().rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    nodes=nodes.merge(entry_pathways,on='Entry',how='left')
    nodes['P']=goaP_user_cut_off['P']
    nodes_file_name='nodes_GOA_Enrichment_Analysis_Process_'+User_method+'_'+file_name_value+'.csv'
    nodes.drop_duplicates().to_csv(nodes_file_name,index=None)
    #
    ## Open R script from github and run
    r_script=requests.get("https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Enrichment_Plots.R").content.decode()
    os.makedirs(dir_goa+'/Process/')
    R_script_enrich = re.sub("./plots/",dir_goa+'/Process/',r_script)
    R_script_enrich = re.sub('legend = "Pathways"','legend = "Process"',R_script_enrich)
    R_script_enrich = re.sub('"Pathway"','"Process"',R_script_enrich)
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
    #print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics for Process (GOA)\n▬▬▬▬▬▬▬▬\n..........\n')
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from subprocess import call
    import shlex, subprocess
    process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
    process.wait()
    if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, oa)
    if os.path.exists(edges_file_name): shutil.move(edges_file_name, oa)
    if os.path.exists(xlsx): shutil.move(xlsx, oa)
    if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
    if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
    if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")
else:
    print('There are not enrichment terms for Biological Process (GOA)')

## Data editing goa function
if goaF_user_cut_off['GO'].count() >= 1:
    enriched_genes=DataFrame(goaF_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
    equalize=str(goaF_user_cut_off.GO*goaF_user_cut_off.Freq)
    id_pathways=re.findall('GO........',equalize)
    enriched_pathways=DataFrame(id_pathways)
    frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
    gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 
    Gene_Entry_Pathway_Enriched=gene_kegg.merge(go_obo_columns,on='GO',how='left')
    ## File name for edges goa function
    edges_file_name='edges_GOA_Enrichment_Analysis_Function_'+User_method+'_'+file_name_value+'.csv'
    Gene_Entry_Pathway_Enriched[['GO','Entry']].to_csv(edges_file_name,index=None)
    xlsx='GOA_Enrichment_Analysis_Function_'+User_method+'_'+file_name_value+'.xlsx'
    Gene_Entry_Pathway_Enriched.to_csv(xlsx,index=None)
    ## Preparing data frame for nodes file
    frame=Gene_Entry_Pathway_Enriched[['GO']].rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    frame['Freq']=goaF_user_cut_off['Freq']
    frame['num']=frame.reset_index().index + 1
    df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Entry']],index=None)
    df2=pd.read_csv(StringIO(df1))
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
    #
    entry_exp['Freq']=float(goaF_user_cut_off['Freq'].min()*.5)
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
    entry_pathways=Gene_Entry_Pathway_Enriched[['GO','Term']].drop_duplicates().rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    nodes=nodes.merge(entry_pathways,on='Entry',how='left')
    nodes['P']=goaF_user_cut_off['P']
    nodes_file_name='nodes_GOA_Enrichment_Analysis_Function_'+User_method+'_'+file_name_value+'.csv'
    nodes.drop_duplicates().to_csv(nodes_file_name,index=None)
    #
    ## Open R script from github and run
    r_script=requests.get("https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Enrichment_Plots.R").content.decode()
    os.makedirs(dir_goa+'/Function/')
    R_script_enrich = re.sub("./plots/",dir_goa+'/Function/',r_script)
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
    #print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics for Function (GOA)\n▬▬▬▬▬▬▬▬\n..........\n')
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from subprocess import call
    import shlex, subprocess
    process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
    process.wait()
    if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, oa)
    if os.path.exists(edges_file_name): shutil.move(edges_file_name, oa)
    if os.path.exists(xlsx): shutil.move(xlsx, oa)
    if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
    if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
    if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")
else:
    print('There are not enrichment terms for Molecular Function (GOA)')
    
## Data editing goa component
if goaC_user_cut_off['GO'].count() >= 1:
    enriched_genes=DataFrame(goaC_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
    equalize=str(goaC_user_cut_off.GO*goaC_user_cut_off.Freq)
    id_pathways=re.findall('GO........',equalize)
    enriched_pathways=DataFrame(id_pathways)
    frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
    gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 
    Gene_Entry_Pathway_Enriched=gene_kegg.merge(go_obo_columns,on='GO',how='left')
    ## File name for edges goa component
    edges_file_name='edges_GOA_Enrichment_Analysis_Component_'+User_method+'_'+file_name_value+'.csv'
    Gene_Entry_Pathway_Enriched[['GO','Entry']].to_csv(edges_file_name,index=None)
    xlsx='GOA_Enrichment_Analysis_Component_'+User_method+'_'+file_name_value+'.xlsx'
    Gene_Entry_Pathway_Enriched.to_csv(xlsx,index=None)
    ## Preparing data frame for nodes file
    frame=Gene_Entry_Pathway_Enriched[['GO']].rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    frame['Freq']=goaC_user_cut_off['Freq']
    frame['num']=frame.reset_index().index + 1
    df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Entry']],index=None)
    df2=pd.read_csv(StringIO(df1))
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
    #
    entry_exp['Freq']=float(goaC_user_cut_off['Freq'].min()*.5)
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
    entry_pathways=Gene_Entry_Pathway_Enriched[['GO','Term']].drop_duplicates().rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    nodes=nodes.merge(entry_pathways,on='Entry',how='left')
    nodes['P']=goaC_user_cut_off['P']
    nodes_file_name='nodes_GOA_Enrichment_Analysis_Component_'+User_method+'_'+file_name_value+'.csv'
    nodes.drop_duplicates().to_csv(nodes_file_name,index=None)
    #
    ## Open R script from github and run
    r_script=requests.get("https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Enrichment_Plots.R").content.decode()
    os.makedirs(dir_goa+'/Component/')
    R_script_enrich = re.sub("./plots/",dir_goa+'/Component/',r_script)
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
    #print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics fof Component (GOA)\n▬▬▬▬▬▬▬▬\n..........\n')
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from subprocess import call
    import shlex, subprocess
    process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
    process.wait()
    if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, oa)
    if os.path.exists(edges_file_name): shutil.move(edges_file_name, oa)
    if os.path.exists(xlsx): shutil.move(xlsx, oa)
    if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
    if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
    if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")
else:
    print('There are not enrichment terms for Cellular Component (GOA)')

##############################################################
###################       Uniprot          ###################
##############################################################

## Data editing uniprot process
if uniP_user_cut_off['GO'].count() >= 1:
    enriched_genes=DataFrame(uniP_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
    equalize=str(uniP_user_cut_off.GO*uniP_user_cut_off.Freq)
    id_pathways=re.findall('GO........',equalize)
    enriched_pathways=DataFrame(id_pathways)
    frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
    gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 
    Gene_Entry_Pathway_Enriched=gene_kegg.merge(go_obo_columns,on='GO',how='left')
    ## File name for edges uniprot process
    edges_file_name='edges_Uniprot_Enrichment_Analysis_Process_'+User_method+'_'+file_name_value+'.csv'
    Gene_Entry_Pathway_Enriched[['GO','Entry']].to_csv(edges_file_name,index=None)
    xlsx='Uniprot_Enrichment_Analysis_Process_'+User_method+'_'+file_name_value+'.xlsx'
    Gene_Entry_Pathway_Enriched.to_csv(xlsx,index=None)
    ## Preparing data frame for nodes file
    frame=Gene_Entry_Pathway_Enriched[['GO']].rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    frame['Freq']=uniP_user_cut_off['Freq']
    frame['num']=frame.reset_index().index + 1
    df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Entry']],index=None)
    df2=pd.read_csv(StringIO(df1))
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
    #
    entry_exp['Freq']=float(uniP_user_cut_off['Freq'].min()*.5)
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
    entry_pathways=Gene_Entry_Pathway_Enriched[['GO','Term']].drop_duplicates().rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    nodes=nodes.merge(entry_pathways,on='Entry',how='left')
    nodes['P']=uniP_user_cut_off['P']
    nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Process_'+User_method+'_'+file_name_value+'.csv'
    nodes.drop_duplicates().to_csv(nodes_file_name,index=None)
    #
    ## Open R script from github and run
    r_script=requests.get("https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Enrichment_Plots.R").content.decode()
    os.makedirs(dir_uniprot+'/Process/')
    R_script_enrich = re.sub("./plots/",dir_uniprot+'/Process/',r_script)
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
    #print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics for Process (Uniprot)\n▬▬▬▬▬▬▬▬\n..........\n')
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from subprocess import call
    import shlex, subprocess
    process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
    process.wait()
    if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, uu)
    if os.path.exists(edges_file_name): shutil.move(edges_file_name, uu)
    if os.path.exists(xlsx): shutil.move(xlsx, uu)
    if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
    if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
    if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")
else:
    print('There are not enrichment terms for Biological Process (Uniprot)')

## Data editing uniprot function
if uniF_user_cut_off['GO'].count() >= 1:
    enriched_genes=DataFrame(uniF_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
    equalize=str(uniF_user_cut_off.GO*uniF_user_cut_off.Freq)
    id_pathways=re.findall('GO........',equalize)
    enriched_pathways=DataFrame(id_pathways)
    frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
    gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 
    Gene_Entry_Pathway_Enriched=gene_kegg.merge(go_obo_columns,on='GO',how='left')
    ## File name for edges uniprot function
    edges_file_name='edges_Uniprot_Enrichment_Analysis_Function_'+User_method+'_'+file_name_value+'.csv'
    Gene_Entry_Pathway_Enriched[['GO','Entry']].to_csv(edges_file_name,index=None)
    xlsx='Uniprot_Enrichment_Analysis_Function_'+User_method+'_'+file_name_value+'.xlsx'
    Gene_Entry_Pathway_Enriched.to_csv(xlsx,index=None)
    ## Preparing data frame for nodes file
    frame=Gene_Entry_Pathway_Enriched[['GO']].rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    frame['Freq']=uniF_user_cut_off['Freq']
    frame['num']=frame.reset_index().index + 1
    df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Entry']],index=None)
    df2=pd.read_csv(StringIO(df1))
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
    #
    entry_exp['Freq']=float(uniF_user_cut_off['Freq'].min()*.5)
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
    entry_pathways=Gene_Entry_Pathway_Enriched[['GO','Term']].drop_duplicates().rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    nodes=nodes.merge(entry_pathways,on='Entry',how='left')
    nodes['P']=uniF_user_cut_off['P']
    nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Function_'+User_method+'_'+file_name_value+'.csv'
    nodes.drop_duplicates().to_csv(nodes_file_name,index=None)
    #
    ## Open R script from github and run
    r_script=requests.get("https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Enrichment_Plots.R").content.decode()
    os.makedirs(dir_uniprot+'/Function/')
    R_script_enrich = re.sub("./plots/",dir_uniprot+'/Function/',r_script)
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
    #print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics Function (Uniprot)\n▬▬▬▬▬▬▬▬\n..........\n')
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from subprocess import call
    import shlex, subprocess
    process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
    process.wait()
    if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, uu)
    if os.path.exists(edges_file_name): shutil.move(edges_file_name, uu)
    if os.path.exists(xlsx): shutil.move(xlsx, uu)
    if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
    if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
    if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")
else:
    print('There are not enrichment terms for Molecular Function (Uniprot)')

## Data editing uniprot component
if uniC_user_cut_off['GO'].count() >= 1:
    enriched_genes=DataFrame(uniC_user_cut_off['Entry'].str.extractall('(/[A-Za-z0-9-_]{0,20}/)')).replace({'/':''},regex=True).reset_index(drop=True)
    equalize=str(uniC_user_cut_off.GO*uniC_user_cut_off.Freq) #>>>
    id_pathways=re.findall('GO........',equalize)
    enriched_pathways=DataFrame(id_pathways)
    frames_genes_pathways_enriched=[enriched_genes,enriched_pathways]
    gene_kegg=pd.concat(frames_genes_pathways_enriched, axis=1, ignore_index=True).rename(columns={0:"Entry",1:'GO'}, index=str)#.reset_index(drop=True) 
    Gene_Entry_Pathway_Enriched=gene_kegg.merge(go_obo_columns,on='GO',how='left')
    ## File name for edges uniprot component
    edges_file_name='edges_Uniprot_Enrichment_Analysis_Component_'+User_method+'_'+file_name_value+'.csv'
    Gene_Entry_Pathway_Enriched[['GO','Entry']].to_csv(edges_file_name,index=None)
    xlsx='Uniprot_Enrichment_Analysis_Component_'+User_method+'_'+file_name_value+'.xlsx'
    Gene_Entry_Pathway_Enriched.to_csv(xlsx,index=None)
    ## Preparing data frame for nodes file
    frame=Gene_Entry_Pathway_Enriched[['GO']].rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    frame['Freq']=uniC_user_cut_off['Freq'] #>>>
    frame['num']=frame.reset_index().index + 1
    df1=DataFrame.to_csv(Gene_Entry_Pathway_Enriched[['Entry']],index=None)
    df2=pd.read_csv(StringIO(df1))
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
    #
    entry_exp['Freq']=float(uniC_user_cut_off['Freq'].min()*.5) #>>>
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],axis=0,ignore_index=True,verify_integrity=True)
    entry_pathways=Gene_Entry_Pathway_Enriched[['GO','Term']].drop_duplicates().rename(columns={'GO':'Entry'},index=str).drop_duplicates().reset_index(drop=True)
    nodes=nodes.merge(entry_pathways,on='Entry',how='left')
    nodes['P']=uniC_user_cut_off['P'] #>>>
    nodes_file_name='nodes_Uniprot_Enrichment_Analysis_Component_'+User_method+'_'+file_name_value+'.csv'
    nodes.drop_duplicates().to_csv(nodes_file_name,index=None)
    #
    ## Open R script from github and run
    r_script=requests.get("https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Enrichment_Plots.R").content.decode()
    os.makedirs(dir_uniprot+'/Component/')
    R_script_enrich = re.sub("./plots/",dir_uniprot+'/Component/',r_script)
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
    #print('▬▬▬▬▬▬▬▬\nStep 7: Generating graphics for Component (Uniprot)\n▬▬▬▬▬▬▬▬\n..........\n')
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from subprocess import call
    import shlex, subprocess
    process = subprocess.Popen(['R', 'CMD', 'BATCH', 'Enrichment_Plots.R'])
    process.wait()
    if os.path.exists(nodes_file_name): shutil.move(nodes_file_name, uu)
    if os.path.exists(edges_file_name): shutil.move(edges_file_name, uu)
    if os.path.exists(xlsx): shutil.move(xlsx, uu)
    if os.path.exists("Rplots.pdf"): os.remove("Rplots.pdf")
    if os.path.exists("Enrichment_Plots.R"): os.remove("Enrichment_Plots.R")
    if os.path.exists("Enrichment_Plots.Rout"): os.remove("Enrichment_Plots.Rout")
else:
    print('There are not enrichment terms for Cellular Component (Uniprot)')

if os.path.exists("enrichment.py"): os.remove("enrichment.py")
import datetime
end = datetime.datetime.now()
print('\n- Work done:   [ '+dir_name_plots+' ]'+
      '\n- Initial hour [ '+str(start.hour)+':'+str(start.minute)+' ]'+
      '\n- End hour     [ '+str(end.hour)+':'+str(end.minute)+' ]'+
      '\n- Total time   [ '+str(end.minute-start.minute)+' minutes ]\n')

