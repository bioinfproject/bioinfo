
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
warnings.filterwarnings("ignore")
from datetime import datetime 
inicio_total = datetime.now()
import os, fnmatch
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


# In[ ]:


## KEGG organisms list
kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()
kegg_organism=pd.read_csv(StringIO(kegg_orgs),names=['T_number','Prefix','Organism','Group'],sep='\t')


# In[ ]:


organism=input('\nStep 1: Enter an annotated gender in the KEGG database\n (e.g., Arabidopsis/Penicillium)\n\n=====> : ')
if organism == '':
    organism = 'xxxxxxxxxxxx'
species=[] # list 1
spec=[] # list 2
for index, row in kegg_organism.iterrows():
    match=re.findall(organism,row['Organism'])
    if match == []:
        match = match
    else:
        #print('\n',row['T_number']+'\t'+row['Prefix']+'\t'+row['Organism']+'\t\t'+row['Group'])
        species.append(row['T_number']+'\t\t'+row['Prefix']+'\t\t'+row['Organism']+'\t\t'+row['Group'])
        spec.append([row['T_number'],row['Prefix'],row['Organism'],row['Group']])
             
if species == []:
    print('\n!!!!!!! Organism not found !!!!!!!\n')
    organism=input('\nStep 1: Enter an annotated gender in the KEGG database\n (e.g., Arabidopsis/Penicillium)\n\n=====> : ')
    if organism == '':
        organism = 'xxxxxxxxxxxx'
    species=[] # list 1
    spec=[] # list 2
    for index, row in kegg_organism.iterrows():
        match=re.findall(organism,row['Organism'])
        if match == []:
            match = match
        else:
            #print('\n',row['T_number']+'\t'+row['Prefix']+'\t'+row['Organism']+'\t\t'+row['Group'])
            species.append(row['T_number']+'\t'+row['Prefix']+'\t'+row['Organism']+'\t\t'+row['Group'])
            spec.append([row['T_number'],row['Prefix'],row['Organism'],row['Group']])
        
    if species == []:
        print('\n!!!!!!! Organism not found !!!!!!!\n')
        sys.exit()
    else:
        #print('mostrar las opciones de especies')
        import tkinter as tk
        import tkinter.ttk as ttk

        root = Tk()
        root.title("NeVOmics")
        root.geometry("850x500")
        root.overrideredirect(1)
        Label(root, text="").pack()
        image_url = "https://raw.githubusercontent.com/eduardo1011/Programas/master/NeVOmics_logo.gif"
        image_byt = urlopen(image_url).read()
        image_b64 = base64.encodestring(image_byt)
        logo1 = tk.PhotoImage(data=image_b64)
        w1 = tk.Label(root, image=logo1).pack()

        Label(root, text="\nIf you use NeVOmics in your research, please cite:      ",font=("Arial", 11)).pack()
        Label(root, text="NeVOmics: an enrichment tool for gene ontology and functional\n"+
        "network analysis and visualization of data from OMICs technologies",font=("Arial", 11),bg="pale green").pack()
        Label(root, text="\n[ Select a species ]\n\n",font=("Arial", 15)).pack()
        Label(root, text="  T number                Org code                         Full name                                      Lineage                                             ",
              font=("Arial", 10, "bold"),bg="pink").place(x=75,y=350)
        
        yyy = tk.StringVar(root)
        xxx = tk.OptionMenu(root, yyy, *species).pack()
        yyy.set(species[0]) # default value

        def ok():
            yyy.get()
    
        button = Button(root,text="          Submit          ",bg="black", fg="white",font=("Arial", 11), command=root.destroy)
        button.pack(pady=15)
        mainloop()
        user_organism = yyy.get()
        print('\nSelected organism:',[user_organism.split("\t")[0],
                                   user_organism.split("\t")[1],
                                   user_organism.split("\t")[2],
                                   user_organism.split("\t")[-1]])
        Prefix = user_organism.split("\t")[2]
else:
    import tkinter as tk
    import tkinter.ttk as ttk

    root = Tk()
    root.title("NeVOmics")
    root.geometry("850x500")
    root.overrideredirect(1)
    Label(root, text="").pack()
    image_url = "https://raw.githubusercontent.com/eduardo1011/Programas/master/NeVOmics_logo.gif"
    image_byt = urlopen(image_url).read()
    image_b64 = base64.encodestring(image_byt)
    logo1 = tk.PhotoImage(data=image_b64)
    w1 = tk.Label(root, image=logo1).pack()

    Label(root, text="\nIf you use NeVOmics in your research, please cite:      ",font=("Arial", 11)).pack()
    Label(root, text="NeVOmics: an enrichment tool for gene ontology and functional\n"+
    "network analysis and visualization of data from OMICs technologies",font=("Arial", 11),bg="pale green").pack()
    Label(root, text="\n[ Select a species ]\n\n",font=("Arial", 15)).pack()
    Label(root, text="  T number            Org code                     Full name                                  Lineage                                             ",
          font=("Arial", 10, "bold"),bg="pink").place(x=75,y=350)

    yyy = tk.StringVar(root)
    xxx = tk.OptionMenu(root, yyy, *species).pack()
    yyy.set(species[0]) # default value

    def ok():
        yyy.get()

    button = Button(root,text="          Submit          ",bg="black", fg="white",font=("Arial", 11), command=root.destroy)
    button.pack(pady=15)
    mainloop()
    user_organism = yyy.get()
    print('\nSelected organism:',[user_organism.split("\t")[0],
                                   user_organism.split("\t")[2],
                                   user_organism.split("\t")[4],
                                   user_organism.split("\t")[-1]])
    Prefix = user_organism.split("\t")[2]
    
############

## Control of  directories
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
###


# In[ ]:


#
html=requests.get('https://www.genome.jp/kegg-bin/show_organism?org='+Prefix).content.decode()
ncbi_proteome_url = ''.join(re.findall('"ftp.*">',html))
ncbi_proteome_url = re.sub('^"|">.*','',ncbi_proteome_url)
ncbi_proteome_url = re.sub('ftp://ftp','https://ftp',ncbi_proteome_url) # convert to https

# url for download proteome from NCBI
url_download_proteome = ncbi_proteome_url+'/'+ncbi_proteome_url.split("/")[-1]+'_translated_cds.faa.gz'

os.makedirs('sequences',exist_ok=True)
# Download proteome
#url = "http://www.mywebsite.com/csv-1-0.csv.gz"
filename = 'sequences/'+url_download_proteome.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(url_download_proteome)
    f.write(r.content)
    
# uncompressed and save as bytes in variable
import gzip
uncompressed_file_gz = gzip.open('sequences/'+url_download_proteome.split("/")[-1], 'rb')
uncompressed = uncompressed_file_gz.read()
uncompressed_file_gz.close()
uncompressed=uncompressed.decode() # convert bytes to str

# save in file
fasta= open('sequences/'+Prefix+'.fasta','w')
fasta.write(uncompressed)
fasta.close()

# first Database (makeblastdb) with annotated proteins
subprocess.call(['makeblastdb','-in','sequences/'+Prefix+'.fasta','-dbtype','prot','-parse_seqids','-out','sequences/proteomes'])

# all kegg-id and pathway-id
dd=requests.get('http://rest.kegg.jp/link/pathway/'+Prefix+'').content.decode()
kegg_path_ID=pd.read_csv(StringIO(dd),sep='\t',header=None,names=['Entry_Kegg','GO']).replace({'^'+Prefix+':|path:':''},regex=True)
# all kegg-id and pathway-description
ee=requests.get('http://rest.kegg.jp/list/pathway/'+Prefix+'').content.decode()
kegg_pathways=pd.read_csv(StringIO(ee),sep='\t',header=None,names=['GO','Term']).replace({'path:|- '+organism[0:5]+'.*':''},regex=True)

# obtaining entry fasta and entry kegg ids
x=[]
for i in kegg_path_ID['Entry_Kegg'].drop_duplicates():
    #print(re.findall('>.*'+i,uncompressed))
    x.append(re.findall('>.*'+i,uncompressed))
xx = DataFrame(x).replace({'>....':'',' [[].*=':'\t'},regex=True).rename(columns={0:'a'})
entries_fasta_kegg=xx['a'].str.split('\t', expand=True).rename(columns = lambda x: "string"+str(x+1)).rename(columns={'string1':'Entry_fasta','string2':'Entry_Kegg'}).dropna()

# save entry_fasta in flat file
entries_fasta_kegg[['Entry_fasta']].to_csv('sequences/in_kegg_'+Prefix+'.txt',header=None,index=None)

# Extract (blastdbcmd) fasta sequences of annotated proteins
subprocess.call(['blastdbcmd','-db','sequences/proteomes','-entry_batch','sequences/in_kegg_'+Prefix+'.txt','-out','sequences/in_kegg_'+Prefix+'.fasta'])

# Second Database (makeblastdb) with annotated proteins
subprocess.call(['makeblastdb','-in','sequences/in_kegg_'+Prefix+'.fasta','-dbtype','prot','-parse_seqids','-out','sequences/proteomes'])

## header blastp
header=('qacc','Entry_fasta','qlen','slen','length','score','bitscore','evalue','pident','nident',
                  'mismatch','positive','gaps','gapopen','stitle')

#######################################################

## submit file

import tkinter as tk
from tkinter import filedialog
## Control of input file
print('\n[ Step 2: Submit file in .fasta format ]')
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()
root.destroy()
if file_path == '':
    print('\n!!!!!!! File not found !!!!!!!')
    import tkinter as tk
    from tkinter import filedialog
    ## Control of input file
    print('\n[ Step 2: Submit file in .fasta format ]')
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    root.destroy()
    if file_path == '':
        print('\n!!!!!!! File not found !!!!!!!')
        if os.path.exists(new_folder): shutil.rmtree(new_folder)
        if os.path.exists('data'): shutil.rmtree('data')
        if os.path.exists('sequences'): shutil.rmtree('sequences')
        if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
        if os.path.exists("HD.py"): os.remove("HD.py")
        sys.exit()
    else:
        print('=====> : ',file_path)
        my_file_handle=open(file_path)
        user_identifiers_fasta = DataFrame(re.findall('>[0-9A-Za-z_-]{0,40}',my_file_handle.read())).rename(columns={0:'qacc'}).replace({'>':''},regex=True)
        my_file_handle.close()
else:
    print('=====> : ',file_path)
    my_file_handle=open(file_path)
    user_identifiers_fasta = DataFrame(re.findall('>[0-9A-Za-z_-]{0,40}',my_file_handle.read())).rename(columns={0:'qacc'}).replace({'>':''},regex=True)
    my_file_handle.close()

#######################################################_______________________________________________

####
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
    root.geometry("850x460")
    root.overrideredirect(1)

    if prompt:
        Label(root, text="\nIf you use NeVOmics in your research, please cite:      ",font=("Arial", 11)).pack()
        Label(root, text="NeVOmics: an enrichment tool for gene ontology and functional\n"+
        "network analysis of data from OMICs technologies",font=("Arial", 11),bg="pale green").pack()
        Label(root, text="\n[ Select a method based on your sequences ]\n",font=("Arial", 15)).pack()
    v = IntVar()
    for i, option in enumerate(options):
        Radiobutton(root, text=option,font=("Arial", 13,"bold"),bg="gold2", variable=v, value=i).pack(anchor="n")
    Label(root, text="* !!!!!   It may take several minutes   !!!!!",font=("Arial", 12)).pack()
    Button(root,text="          Submit          ",bg="black", fg="white",font=("Arial", 11), command=root.destroy).pack(pady=15) #anchor must be n, ne, e, se, s, sw, w, nw, or center
    
    Label(root, text='\nThe Basic Local Alignment Search Tool (BLAST) finds regions of local similarity \n'+
'between sequences. The program compares nucleotide or protein sequences to sequence databases\n'+
'and calculates the statistical significance of matches. BLAST can be used to infer functional\n'+
'and evolutionary relationships between sequences as well as help identify members of gene families.',font=("Arial", 10)).pack()
    import webbrowser
    def callback(event):
        webbrowser.open_new(r"https://blast.ncbi.nlm.nih.gov/Blast.cgi")
    link = Label(root, text="Basic Local Alignment Search Tool", fg="blue", cursor="hand2")
    link.pack()
    link.bind("<Button-1>", callback)
    
    root.mainloop()						             #pady  separación entre la tercer linea y submit
    if v.get() == 0: return None
    return options[v.get()]

result = ask_multiple_choice_question(" ",
    [
        "   None     ",
        "   BLASTP : search protein databases using a protein query       ",
        "   BLASTX : search protein databases using a translated nucleotide query       *"
    ]
)
if format(repr(result)) == 'None':
    print('\nNo method was selected')
    if os.path.exists(new_folder): shutil.rmtree(new_folder)
    if os.path.exists('data'): shutil.rmtree('data')
    if os.path.exists('sequences'): shutil.rmtree('sequences')
    if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
    if os.path.exists("HD.py"): os.remove("HD.py")
    sys.exit()

if re.findall('[A-Z]{6}',format(repr(result)))[0] == 'BLASTP':
    method_blast = 'blastp'
else:
    method_blast = 'blastx'
##################################################________________________________________________
## blastp
subprocess.call([method_blast,'-db','sequences/proteomes','-query', file_path,'-evalue','1E-6','-outfmt',
                 '6 qacc sacc qlen slen length score bitscore evalue pident nident mismatch positive gaps gapopen stitle',
                 '-max_target_seqs','1','-max_hsps','1','-out','sequences/'+Prefix+'.tab'])

# open blastp results
blastp=pd.read_csv('sequences/'+Prefix+'.tab',sep='\t',names=header)

#filtro 70% de identidad
blastp_cut_off_70=blastp[(blastp.pident >= 70) & (blastp.pident <= 100)].reset_index(drop=True).drop_duplicates()

if float(blastp_cut_off_70['qacc'].count()) > 0:
    print('\n* BLAST Results:',int(float(blastp_cut_off_70['qacc'].count())),'Proteins found with >= 70% Identity\n')
    #blasp_qacc_Entry_fasta=blastp_cut_off_70[['qacc','Entry_fasta']]
else:
    print('\n!!!!!!!  No sequences with more than 70% identity were found  !!!!!!!')
    print('\n!!!!!!!  Finished process  !!!!!!!\n')
    sys.exit()
    
# orthologs found
orthologs = pd.merge(blastp_cut_off_70[['qacc','Entry_fasta']],entries_fasta_kegg,on='Entry_fasta',how='left')
orthologs=pd.merge(orthologs,kegg_path_ID,on='Entry_Kegg',how='left')    


# In[ ]:


## Create a folder
os.makedirs('data',exist_ok=True)

kegg_pathways.to_csv('data/Pathways.txt',sep='\t',index=None)

# 1.- Preparation of background
entries_fasta_kegg[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/Background.txt',index=None)

# 2.- Preparation of list with pathways
orthologs[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/List.txt',index=None)

# 3.- background with: Entry	GO, for association file
kegg_path_ID.rename(columns={'Entry_Kegg':'Entry'}).to_csv('data/Association.txt',index=None,sep='\t')

#exploratory analysis of P-value in data
file='Pathways.txt'
fdr_min=0.05
subprocess.call(["python","HD.py",file,str(fdr_min),'data/Enrichment_analysis_Path.tsv'])
enrich_P=pd.read_csv('data/Enrichment_analysis_Path.tsv',sep='\t')

#
if enrich_P[(enrich_P.P < 0.05)]['P'].count() >= 1:
    print('\n----------------------------------------------------\n'+
            '                KEGG Pathways Annotation\n----------------------------------------------------'+
            '\n\n*****  KEGG Pathways:  ',enrich_P[(enrich_P.P < 0.05)]['P'].count(),
            'Pathways were found with P-value < 0.05')
            
    ## Choice of correction corection method
    methods = ['Bonferroni','FDR']
    User_method=input('\n[ Step 3: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
    method_P=[]
    for x in methods:
        if x == User_method:
            method_P.append(x)
    if method_P == []:
        print('\n!!!!! Incorrect Method !!!!!')
        User_method=input('\n[ Step 3: Choose a correction method (e.g., FDR / Bonferroni) ]\n=====> : ')
        for x in methods:
            if x == User_method:
                method_P.append(x)
    if method_P == []:
        print('\n!!!!! Incorrect Method !!!!!')
        if os.path.exists(new_folder): shutil.rmtree(new_folder)
        if os.path.exists('sequences'): shutil.rmtree('sequences')
        if os.path.exists('data'): shutil.rmtree('data')
        if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
        if os.path.exists("HD.py"): os.remove("HD.py")
        sys.exit()
    ########################################
    ########################################
    ## loop control User value
    a=input('[ Step 4: Choose a Value (e.g., 0.05) ]\n=====> : ')
    if a == '':
        a='aaa'
    if re.search(r'[A-Za-z]{1,10}',a):
        print('Incorrect value')
        a=input('[ Step 4: Choose a Value (e.g., 0.05) ]\n=====> : ')
        if a == '':
            a='aaa'
        if re.search(r'[A-Za-z]{1,10}',a):
            print('Incorrect value')
            if os.path.exists(new_folder): shutil.rmtree(new_folder)
            if os.path.exists('sequences'): shutil.rmtree('sequences')
            if os.path.exists('data'): shutil.rmtree('data')
            if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
            if os.path.exists("HD.py"): os.remove("HD.py")
            sys.exit()
        else:
            if 0.005 <= float(a) <= 0.5:
                User_value_P=float(a)
            else:
                print('Incorrect value')
                if os.path.exists(new_folder): shutil.rmtree(new_folder)
                if os.path.exists('sequences'): shutil.rmtree('sequences')
                if os.path.exists('data'): shutil.rmtree('data')
                if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
                if os.path.exists("HD.py"): os.remove("HD.py")
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
                if os.path.exists('sequences'): shutil.rmtree('sequences')
                if os.path.exists('data'): shutil.rmtree('data')
                if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
                if os.path.exists("HD.py"): os.remove("HD.py")
                sys.exit()
            else:
                if 0.005 <= float(a) <= 0.5:
                    User_value_P=float(a)
                else:
                    print('Incorrect value')
                    if os.path.exists(new_folder): shutil.rmtree(new_folder)
                    if os.path.exists('sequences'): shutil.rmtree('sequences')
                    if os.path.exists('data'): shutil.rmtree('data')
                    if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
                    if os.path.exists("HD.py"): os.remove("HD.py")
                    sys.exit()
    ########################################

########################################
analysis='Pathways.txt'
if ''.join(method_P) == 'Bonferroni':
    subprocess.call(["python","HD.py",
                        analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
    enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
    results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
    proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
    GO_count_P=results_process_P[['GO']].count()[0]

else:
    if ''.join(method_P) == 'FDR':
        subprocess.call(["python","HD.py",
                            analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
        enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
        results_process_P=enrich_P[enrich_P.Sig == 'T']
        proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
        GO_count_P=results_process_P[['GO']].count()[0]
        
#####
# °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
if results_process_P['GO'].count() >= 1:
    ##
    ## 
    non_annoted = pd.DataFrame.merge(user_identifiers_fasta,orthologs,how="left", on='qacc').fillna('N')
    lis=non_annoted[non_annoted.Entry_Kegg == 'N'][['qacc']]
    lis['QACC']='QACC'
    lis['ent']=lis[['qacc']].replace({'$':'; '},regex=True)
    lis = lis.groupby('QACC')['ent'].sum().reset_index()
    # following the GeneMerge1.4 approach we use non-singletons as multiple testing value
    uuu = orthologs.groupby('GO')['qacc'].count().reset_index().sort_values(by ='qacc',ascending=False).reset_index(drop=True).drop_duplicates()
    
    report = ['\n\t\n'+
                'NeVOmics\t'+new_folder+
                '\n\nInput file name\t'+file_path+
                '\nAssociation file name\t'+analysis+
                '\nTotal number of background\t'+str(entries_fasta_kegg['Entry_Kegg'].drop_duplicates().count())+
                '\nTotal fasta sequences\t'+str(user_identifiers_fasta['qacc'].drop_duplicates().count())+
                '\n\nBackground with Pathways\t'+str(entries_fasta_kegg['Entry_Kegg'].drop_duplicates().count())+
                '\nSequences with orthologous Pathways\t'+str(orthologs['qacc'].drop_duplicates().count())+
                '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+
                '\nCorrection Method\t'+''.join(method_P)+
                '\nValue\t'+str(User_value_P)+
                '\n\t\n'+
                '\nSequences with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.Entry_Kegg == 'N'][['qacc']].count()[0])+
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
    process_goa=process_goa.merge(orthologs[['qacc','Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}),on='Entry',how='left')
    df2=process_goa[['Entry']]
    process_goa=pd.merge(process_goa,kegg_pathways[['GO','Term']].rename(columns={'GO':'Path','Term':'Description'}),on='Path',how='left').drop_duplicates()
    ## Excel file with all information
    writer = pd.ExcelWriter(new_folder+'/Enrichment_Pathways_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
    combine.to_excel(writer,'Significant KEGG Pathways',index=False)
    enrich_P.to_excel(writer,'Enrichment Results',index=False)
    process_goa.to_excel(writer,'Edges Pathways',index=False)
    writer.save()
    
    frame=results_process_P[['GO','go_list','Rank','Term','P','FDR']].rename(columns={'GO':'Entry','go_list':'Freq','Rank':'num'})
    #
    #
    entry_exp=df2
    entry_exp['Freq']=float(frame['Freq'].min()*.5)
    entry_exp['num']=float(frame['num'].count() + 1)
    nodes=pd.concat([frame,entry_exp],sort=False)
    nodes_file_name='nodes_KEGG_Enrichment_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
    nodes.drop_duplicates().to_csv(level_1_kegg+nodes_file_name,index=None)
    
    ## Open R script from github and run
    r_script=requests.get('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/Plots.R').content.decode()
    R_script_enrich = re.sub('qwertyuiop',level_1_kegg+nodes_file_name,r_script) # name edges file
    R_script_enrich = re.sub('asdfghjkl',level_1_kegg+edges_file_name,R_script_enrich) # name nodes file
    R_script_enrich = re.sub('zxcvbnm',level_2_kegg,R_script_enrich) # store plots
    R_script_enrich = re.sub('ASPECT','KEGG Pathways',R_script_enrich) # store plots
    R_script_enrich = re.sub('poiuytrewq','KEGG Pathways Annotation',R_script_enrich) # legend
    f= open(level_1_kegg+'/Kegg_Enrichment_Plots.R','w')
    f.write(R_script_enrich)
    f.close()
else:
    print('\n*****  KEGG Pathways: No significant terms were found')
    ##
    non_annoted = pd.DataFrame.merge(user_identifiers_fasta,orthologs,how="left", on='qacc').fillna('N')
    lis=non_annoted[non_annoted.Entry_Kegg == 'N'][['qacc']]
    lis['QACC']='QACC'
    lis['ent']=lis[['qacc']].replace({'$':'; '},regex=True)
    lis = lis.groupby('QACC')['ent'].sum().reset_index()
    # following the GeneMerge1.4 approach we use non-singletons as multiple testing value
    uuu = orthologs.groupby('GO')['qacc'].count().reset_index().sort_values(by ='qacc',ascending=False).reset_index(drop=True).drop_duplicates()
    
    report = ['\n\t\n'+
                'NeVOmics\t'+new_folder+
                '\n\nInput file name\t'+file_path+
                '\nAssociation file name\t'+analysis+
                '\nTotal number of background\t'+str(entries_fasta_kegg['Entry_Kegg'].drop_duplicates().count())+
                '\nTotal fasta sequences\t'+str(user_identifiers_fasta['qacc'].drop_duplicates().count())+
                '\n\nBackground with Pathways\t'+str(entries_fasta_kegg['Entry_Kegg'].drop_duplicates().count())+
                '\nSequences with orthologous Pathways\t'+str(orthologs['qacc'].drop_duplicates().count())+
                '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+
                '\nCorrection Method\t'+''.join(method_P)+
                '\nValue\t'+str(User_value_P)+
                '\n\t\n'+
                '\nSequences with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.Entry_Kegg == 'N'][['qacc']].count()[0])+
                '\n'+lis['ent'][0]]
    rep=''.join(report)
    information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    combine=pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
    
    ## Excel file with all information
    writer = pd.ExcelWriter(new_folder+'/Enrichment_Pathways_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.xlsx')
    combine.to_excel(writer,'Significant KEGG Pathways',index=False)
    enrich_P.to_excel(writer,'Enrichment Results',index=False)
    writer.save()

####
## summary
database='KEGG'
aspect='Pathway'

information=[[database],[aspect],[''.join(method_P)],[str(User_value_P)],[GO_count_P],[proteins_count_P]]
summary = DataFrame(information)
summary=summary.T.reset_index(drop=True).rename(columns={0:'DATABASE',1:'TYPE',2:'US_METH',3:'US_VAL',4:'PATHWAYS',5:'UNIQ_PROT'})

print('\nENRICHMENT SUMMARY')
print('\n',summary)

####
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
            plots_selection.append(i+''.join(fnmatch.filter(os.listdir((i)), '*.R'))) 
        # run R scripts
        for i in plots_selection:
            run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
            run_uni.wait()
            
###
if os.path.exists('data'): shutil.rmtree('data')
if os.path.exists('sequences'): shutil.rmtree('sequences')
if os.path.exists("HD.py"): os.remove("HD.py")
if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
if os.path.exists("./*.RData"): os.remove("./*.RData")
if os.path.exists(".RData"): os.remove(".RData")
if os.path.exists("./Rplots.pdf"): os.remove("./Rplots.pdf")
import os, fnmatch
def find(pattern,path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
file_uniprot=find('*.Rout','./')
for i in file_uniprot:
    if os.path.exists(i): os.remove(i)
#
#
import os, fnmatch
def find(pattern,path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
file_uniprot=find('*.R','./')
for i in file_uniprot:
    if os.path.exists(i): os.remove(i)
    
# print total time of analysis
lapso_total = datetime.now() - inicio_total
print('\n'+new_folder+': Analysis Time (hh:mm:ss.ms) {}'.format(lapso_total),'\n')

