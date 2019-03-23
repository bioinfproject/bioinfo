# coding: utf-8

# In[1]:


import datetime
start = datetime.datetime.now()
from pandas import Series, DataFrame 
import pandas as pd
from pandas.compat import StringIO
import csv
import gzip
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


# In[2]:


kegg_orgs = requests.get("http://rest.kegg.jp/list/organism").content.decode()
kegg_organism=pd.read_csv(StringIO(kegg_orgs),names=['T_number','Prefix','Organism','Group'],sep='\t')


# In[3]:


organism=input('\nStep 1: Enter an annotated gender in the KEGG database\n (e.g., Arabidopsis/Penicillium)\n\n=====> : ')
encontrado = []
for index, row in kegg_organism.iterrows():
    finder = row['Organism'].split(' ')[0]
    if organism == finder:
        encontrado.append(row.T_number)
    else:
        pass
if encontrado == []:
    print('\n!!!!!!! Organism not found !!!!!!!\n')
    organism=input('\nStep 1: Enter an annotated gender in the KEGG database\n (e.g., Arabidopsis/Penicillium)\n\n=====> : ')
    encontrado = []
    for index, row in kegg_organism.iterrows():
        finder = row['Organism'].split(' ')[0]
        if organism == finder:
            encontrado.append(row.T_number)
        else:
            pass
    if encontrado == []:
        print('\n!!!!!!! Organism not found !!!!!!!\n')
        sys.exit()
    else:
        species = []
        for i in encontrado:
            spe = kegg_organism[kegg_organism.T_number == i]
            for index, row in spe.iterrows():
                species.append(row['T_number']+'\t'+row['Prefix']+'\t'+row['Organism']+'\t\t'+row['Group'])
else:
    species = []
    for i in encontrado:
        spe = kegg_organism[kegg_organism.T_number == i]
        for index, row in spe.iterrows():
            species.append(row['T_number']+'\t'+row['Prefix']+'\t'+row['Organism']+'\t\t'+row['Group'])

#######
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
print('\nSelected organism:',user_organism.split("\t")[2],'\n')
Prefix = user_organism.split("\t")[1]


# In[4]:


############

## Control of  directories
ls=[]
for i in os.listdir("./"):
    if re.search('job_KEGG_BLAST_[0-9]{1,3}',str(i)):
        ls.append(i)
if ls == []:
    new_folder='job_KEGG_BLAST_1'
    xoxo='1'
    os.makedirs('job_KEGG_BLAST_1/job_KEGG_BLAST_1/job_KEGG_BLAST_plots_1',exist_ok=True)
    level_1_kegg=new_folder+'/job_KEGG_BLAST_'+xoxo+'/'
    level_2_kegg=new_folder+'/job_KEGG_BLAST_'+xoxo+'/job_KEGG_BLAST_plots_'+xoxo+'/'
else:
    n=max([int(x) for x in re.findall('[0-9]{1,3}',str(ls))])
    xoxo=str(n+1)
    n=str(n)
    old_folder=''.join(re.findall('job_KEGG_BLAST_'+n,str(ls)))
    new_folder=re.sub(n,xoxo,old_folder)
    os.makedirs(new_folder+'/job_KEGG_BLAST_'+xoxo+'/job_KEGG_BLAST_plots_'+xoxo+'',exist_ok=True)
    level_1_kegg=new_folder+'/job_KEGG_BLAST_'+xoxo+'/'
    level_2_kegg=new_folder+'/job_KEGG_BLAST_'+xoxo+'/job_KEGG_BLAST_plots_'+xoxo+'/'
###


# In[5]:
genome = requests.get('https://www.genome.jp/kegg-bin/show_organism?org='+Prefix).content.decode()
org_id = re.sub('mode=Info&id=','',''.join(re.findall('mode=Info&id=[0-9]{1,50}',genome)))
#ncbi_proteome_url = ''.join(re.findall('"ftp.*">',genome))
#ncbi_proteome_url = re.sub('^"|">.*','',ncbi_proteome_url)
#ncbi_proteome_url = re.sub('ftp://ftp','https://ftp',ncbi_proteome_url) # convert to https
#url_download_proteome = ncbi_proteome_url+'/'+ncbi_proteome_url.split("/")[-1]+'_translated_cds.faa.gz'
os.makedirs('sequences',exist_ok=True)


# In[6]:
link = 'https://www.uniprot.org/uniprot/?query=taxonomy:'+org_id+'&format=fasta'
file_name = 'sequences/'+Prefix+'.fasta'
with open(file_name, 'wb') as f:
    #print ("Downloading %s" % file_name)
    response = requests.get(link, stream=True)
    total_length = response.headers.get('X-Total-Results')
    if total_length is None: # no content length header
        f.write(response.content)
    else:
        dl = 0
        total_length = int(total_length)
        for data in response.iter_content(chunk_size=8192):
            dl += len(data)
            f.write(data)
            done = int(0.07 * dl / total_length)
            sys.stdout.write("\rLoading [%s%s ] %s MB" % ('■' * done, ' ' * (10-done), round(dl/1000000,2)), ) 
            sys.stdout.flush() ## ■●


# In[7]:
with open('sequences/'+Prefix+'.fasta', 'r') as seq:
    seq = seq.read()
seq1 = ''.join(re.findall('>.*',seq))
seq2 = '\n'.join(re.findall('([|][A-Z0-9]{1,50}[|]|GN=[A-Z0-9a-z_]{1,100})',seq1))
seq3 = re.sub('\n[|]','\n',seq2)
seq4 = re.sub('[|]\nGN=','\t',seq3)
seq5 = re.sub('^[|]','',seq4)
seq_uniprot = pd.read_csv(StringIO(seq5),sep='\t',names=['Entry_Uniprot','Name'])


# In[8]:
# all kegg-id and pathway-id
dd = requests.get('http://rest.kegg.jp/link/pathway/'+Prefix+'').content.decode()
dd = re.sub(Prefix+':|path:','',dd)
kegg_path_ID = pd.read_csv(StringIO(dd),sep='\t',header=None,names=['Entry_Kegg','GO'])
string = []
for i in kegg_path_ID.Entry_Kegg:
    string.append(str(i))
kegg_path_ID['Entry_Kegg'] = string


# In[9]:
# all kegg-id and pathway-description
ee=requests.get('http://rest.kegg.jp/list/pathway/'+Prefix+'').content.decode()
ee = re.sub('path:|- '+organism[0:5]+'.*','',ee)
kegg_pathways = pd.read_csv(StringIO(ee),sep='\t',header=None,names=['GO','Term'])


# In[10]:
ff = requests.get('http://rest.kegg.jp/conv/uniprot/'+Prefix+'').content.decode()
ff = re.sub(Prefix+':','',ff)
ff = re.sub('up:','',ff)
KeggID_UniprotID  = pd.read_csv(StringIO(ff),sep='\t',header=None,names=['Entry_Kegg','Entry_Uniprot'])
string = []
for i in KeggID_UniprotID.Entry_Kegg:
    string.append(str(i))
KeggID_UniprotID['Entry_Kegg'] = string

# info version
infokegg = requests.get('http://rest.kegg.jp/info/'+Prefix+'').content.decode()
infokegg = ''.join(re.findall('Release .*',infokegg))
infokegg = re.sub('Release ','',infokegg)
# In[50]:
kegg_uniprot = kegg_path_ID.merge(KeggID_UniprotID,on = 'Entry_Kegg', how = 'left').dropna()
kegg_uniprot_all = kegg_uniprot.merge(seq_uniprot, on = 'Entry_Uniprot', how = 'left')


# In[12]:
kegg_uniprot_all.Entry_Uniprot.drop_duplicates().to_csv('sequences/annotated.txt',header=None,index=None)



# https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
# first Database (makeblastdb) with all proteins
print('\n\n*** BLAST database: 1')
db = subprocess.check_output(['makeblastdb','-in','sequences/'+Prefix+'.fasta','-dbtype','prot','-parse_seqids',
                 '-out','sequences/proteomes'])


# In[17]:
# retrieve annotated proteins
db = subprocess.call(['blastdbcmd','-db','sequences/proteomes','-dbtype','prot',
                      '-entry_batch','sequences/annotated.txt','-out','sequences/annotated.fasta'])

# In[18]:
# second Database (makeblastdb) with annotated proteins
print('\n*** BLAST database: 2')
db = subprocess.check_output(['makeblastdb','-in','sequences/annotated.fasta','-dbtype','prot','-parse_seqids',
                 '-out','sequences/annotated'])
# In[19]:
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
        with open(file_path, 'r') as inp:
            inp = inp.read()
        id1 = '\n'.join(re.findall('>[A-Za-z0-9-_|]{1,100}',inp))
        id2 = re.sub('>','',id1)
        user_identifiers_fasta = pd.read_csv(StringIO(id2),header=None,names=['qacc'])
else:
    print('=====> : ',file_path)
    with open(file_path, 'r') as inp:
        inp = inp.read()
    id1 = '\n'.join(re.findall('>[A-Za-z0-9-_|]{1,100}',inp))
    id2 = re.sub('>','',id1)
    user_identifiers_fasta = pd.read_csv(StringIO(id2),header=None,names=['qacc'])


# In[20]:


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
        "   BLASTX : search protein databases using a translated nucleotide query*"
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
    print('')
    version = subprocess.call([method_blast,'-version'])
else:
    method_blast = 'blastx'
    print('')
    version = subprocess.call([method_blast,'-version'])

# In[21]:


## header blastp
header=('qacc','Entry_Uniprot','qlen','slen','length','score','bitscore','evalue','pident','nident',
                  'mismatch','positive','gaps','gapopen','stitle')


blas = subprocess.call([method_blast,'-db','sequences/annotated','-query', file_path,'-evalue','1E-6','-outfmt',
                 '6 qacc sacc qlen slen length score bitscore evalue pident nident mismatch positive gaps gapopen stitle',
                 '-max_target_seqs','10','-max_hsps','1','-out','sequences/'+Prefix+'.tab'])


# In[24]:


# open blastp results
blastp = pd.read_csv('sequences/'+Prefix+'.tab',sep='\t',names=header)

#filtro 70% de identidad
blastp_cut_off_70=blastp[(blastp.pident >= 70) & (blastp.pident <= 100)].reset_index(drop=True).drop_duplicates()
#

# removiendo duplicados
dfs = []
for i in blastp_cut_off_70.qacc.drop_duplicates():
    df = blastp_cut_off_70[blastp_cut_off_70.qacc == i].sort_values(by='pident', ascending=False)
    dfs.append(df[:1])
blastp_cut_off_70 = pd.concat(dfs)

writer = pd.ExcelWriter(new_folder+'/Blast_Results'+method_blast+'.xlsx')
blastp_cut_off_70.to_excel(writer,'Blast_Results',index=False)
writer.save()

if float(blastp_cut_off_70['qacc'].count()) > 0:
    print('\n* BLAST Results:',int(float(blastp_cut_off_70['qacc'].count())),'Proteins found with >= 70% Identity\n')
    #blasp_qacc_Entry_fasta=blastp_cut_off_70[['qacc','Entry_fasta']]
else:
    print('\n!!!!!!!  No sequences with more than 70% identity were found  !!!!!!!')
    print('\n!!!!!!!  Finished process  !!!!!!!\n')
    sys.exit()


# In[27]:


# orthologs found
orthologs = pd.merge(blastp_cut_off_70[['qacc','Entry_Uniprot']],kegg_uniprot_all,on='Entry_Uniprot',how='left')
string = []
for i in orthologs.Entry_Kegg:
    string.append(str(i))
orthologs['Entry_Kegg'] = string

# In[33]:


## Create a folder
os.makedirs('data',exist_ok=True)

kegg_pathways.to_csv('data/Pathways.txt',sep='\t',index=None)

# 1.- Preparation of background
kegg_path_ID[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/Background.txt',index=None)

# 2.- Preparation of list with pathways
orthologs[['Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}).drop_duplicates().to_csv('data/List.txt',index=None)

# 3.- background with: Entry	GO, for association file
kegg_path_ID.rename(columns={'Entry_Kegg':'Entry'}).to_csv('data/Association.txt',index=None,sep='\t')


# In[34]:


#exploratory analysis of P-value in data
file='Pathways.txt'
fdr_min=0.05
subprocess.call(["python3","HD.py",file,str(fdr_min),'data/Enrichment_analysis_Path.tsv'])
enrich_P=pd.read_csv('data/Enrichment_analysis_Path.tsv',sep='\t')


# In[46]:


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
    subprocess.call(["python3","HD.py",
                        analysis,str(User_value_P),'data/Enrichment_Analysis_Path.csv'])
    enrich_P=pd.read_csv('data/Enrichment_Analysis_Path.csv',sep='\t')
    results_process_P=enrich_P[enrich_P.Bonf_corr < User_value_P]
    proteins_count_P=DataFrame(results_process_P['entry'].str.extractall('([A-Za-z0-9-_]{0,20})')).dropna().drop_duplicates().reset_index(drop=True).count()[0]
    GO_count_P=results_process_P[['GO']].count()[0]

else:
    if ''.join(method_P) == 'FDR':
        subprocess.call(["python3","HD.py",
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
    
    # following the GeneMerge1.4 approach we use non-singletons as multiple testing value
    uuu = orthologs.groupby('GO')['qacc'].count().reset_index().sort_values(by ='qacc',ascending=False).reset_index(drop=True).drop_duplicates()
    
    report = ['\n\t\n'+
                'NeVOmics\t'+new_folder+
                '\nKEGG DB Last-Modified\t'+infokegg+
                '\nBlast method\t'+method_blast+
                '\n\nInput file name\t'+file_path+
                '\nAssociation file name\t'+analysis+
                '\nTotal number of background\t'+str(kegg_path_ID['Entry_Kegg'].drop_duplicates().count())+
                '\nTotal fasta sequences\t'+str(user_identifiers_fasta['qacc'].drop_duplicates().count())+
                '\n\nBackground with Pathways\t'+str(kegg_path_ID['Entry_Kegg'].drop_duplicates().count())+
                '\nSequences with orthologous Pathways\t'+str(orthologs['qacc'].drop_duplicates().count())+
                '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+
                '\nCorrection Method\t'+''.join(method_P)+
                '\nValue\t'+str(User_value_P)+
                '\n\t\n'+
                '\nSequences with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.Entry_Kegg == 'N'][['qacc']].count()[0])+
                '\n'+str(';'.join(lis.qacc))]
    rep=''.join(report)
    information=pd.read_csv(StringIO(rep),sep='\t',header=None,names=['GO','go_list'])
    combine=pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
    
    ##------raw
    raw_data = pd.read_csv('data/raw_list.txt',sep='\t')
    process_goa = pd.merge(results_process_P,raw_data,on='GO',how='left')[['GO','Entry']]
    string = []
    for i in process_goa.Entry:
        string.append(str(i))
    process_goa['Entry'] = string
    #
    ## save file with edges for graph
    edges_file_name='edges_KEGG_Enrichment_Analysis_'+''.join(method_P)+'_'+str(User_value_P)+'.csv'
    process_goa.to_csv(level_1_kegg+edges_file_name,index=None)
    
    
    enrich_P=enrich_P.rename(columns={'GO':'Path','go_list':'path_list','go_back':'path_back'})
    process_goa=process_goa.rename(columns={'GO':'Path'})
    ##
    process_goa=process_goa.merge(orthologs[['qacc','Entry_Uniprot','Entry_Kegg']].rename(columns={'Entry_Kegg':'Entry'}),on='Entry',how='left')
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
    
    # following the GeneMerge1.4 approach we use non-singletons as multiple testing value
    uuu = orthologs.groupby('GO')['qacc'].count().reset_index().sort_values(by ='qacc',ascending=False).reset_index(drop=True).drop_duplicates()
    
    report = ['\n\t\n'+
                'NeVOmics\t'+new_folder+
                '\nKEGG DB Last-Modified\t'+infokegg+
                '\nBlast method\t'+method_blast+
                '\n\nInput file name\t'+file_path+
                '\nAssociation file name\t'+analysis+
                '\nTotal number of background\t'+str(kegg_path_ID['Entry_Kegg'].drop_duplicates().count())+
                '\nTotal fasta sequences\t'+str(user_identifiers_fasta['qacc'].drop_duplicates().count())+
                '\n\nBackground with Pathways\t'+str(kegg_path_ID['Entry_Kegg'].drop_duplicates().count())+
                '\nSequences with orthologous Pathways\t'+str(orthologs['qacc'].drop_duplicates().count())+
                '\nNon-singletons value\t'+str(uuu.drop_duplicates().count()[0])+
                '\nCorrection Method\t'+''.join(method_P)+
                '\nValue\t'+str(User_value_P)+
                '\n\t\n'+
                '\nSequences with no information in KEGG Pathways\t'+str(non_annoted[non_annoted.Entry_Kegg == 'N'][['qacc']].count()[0])+
                '\n'+str(';'.join(lis.qacc))]
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


# In[47]:



if format(repr(result)) == 'None':
    print('\n!!!!! Graphics not generated !!!!!')
else:
    print('\nRunning: R')
    folders = [level_1_kegg]
    # Biological process
    if ''.join(re.findall('KEGG',format(repr(result)))) == 'KEGG':
        print('\nWaiting ...')
        # find R scripst process
        plots_selection=[]
        for i in folders:
            plots_selection.append(i+''.join(fnmatch.filter(os.listdir((i)), '*.R'))) 
        # run R scripts
        for i in plots_selection:
            run_uni=subprocess.Popen(['R', 'CMD', 'BATCH', i])
            run_uni.wait()


# In[49]:


# In[ ]:
if os.path.exists('data'): shutil.rmtree('data')
if os.path.exists('sequences'): shutil.rmtree('sequences')
if os.path.exists('GO.py'): os.remove('GO.py')
if os.path.exists("KEGG_BLAST.py"): os.remove("KEGG_BLAST.py")
if os.path.exists('HD.py'): os.remove('HD.py')
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
#file_uniprot=find('*.Rout','./')
#for i in file_uniprot:
#    if os.path.exists(i): os.remove(i)
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
