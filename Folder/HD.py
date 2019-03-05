
# coding: utf-8

# In[ ]:


from scipy.stats import hypergeom
from pandas import Series, DataFrame 
import pandas as pd
from pandas.compat import StringIO
import csv
import re
import shutil, os
import numpy as np
import sys


# In[ ]:


#sys.argv[1]='GO_BP.txt'


# In[ ]:


# open files
association=pd.read_csv('data/Association.txt',sep='\t')
description=pd.read_csv('data/'+sys.argv[1],sep='\t')
background=pd.read_csv('data/Background.txt',sep='\t')
List=pd.read_csv('data/List.txt',sep='\t')


# In[ ]:


# Total of proteins with GO terms in list (for hypergeometric dustribution)
total_protein_list=List['Entry'].drop_duplicates().count()
#total_protein_list


# In[ ]:


list_cat=pd.merge(List,association,on='Entry',how='left').reset_index(drop=True).dropna().drop_duplicates()
yyy=pd.merge(list_cat,description,on='GO',how='left').dropna()
yyy = yyy.drop_duplicates(["GO", "Entry"])
list_category=DataFrame(yyy.groupby('GO').Entry.size()).reset_index()
#print(list_cat['Entry'].drop_duplicates().count())


# In[ ]:


# Total of proteins with GO terms in background
total_proteins_bg=background['Entry'].drop_duplicates().count()#.drop_duplicates().count()
total_proteins_bg


# In[ ]:


# Get all background GO terms involved in Biological Process
category=pd.merge(background,association,on='Entry',how='left').reset_index(drop=True).dropna()
xxx=pd.merge(category,description,on='GO',how='left').dropna()
xxx = xxx.drop_duplicates(["GO", "Entry"])
cat=DataFrame(xxx.groupby('GO').Entry.size()).reset_index()
#print(category['Entry'].drop_duplicates().count())


# In[ ]:


#######################################################################
#########################  Statistical test   #########################
#######################################################################
## Entry_x = proteins in list
## Entry_y = proteins in background by category (P or F or C)
statistics=pd.merge(list_category,cat,on='GO',how='left').rename(columns={'Entry_x':'go_list','Entry_y':'go_back'})

## following the GeneMerge1.4 approach we use non-singletons as multiple testing value
multiple_testing_value=statistics[statistics.go_back > 1].count()[0]
#print(multiple_testing_value)
statistics['tot_list']=total_protein_list
statistics['tot_back']=total_proteins_bg

## Hypergeometric Distribution
#hypergeom.sf(k, M, n, N, loc=0)
# k = number of genes/proteins associated to the process "cell cycle"
# M = total number of genes/proteins with some annotation
# n = total number of genes/proteins annotated for "cell cycle" inside M
# N = number of genes associated to at least one “Biological Process” in the Gene Ontology.
# https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html
# example =  hypergeom.sf(8-1, total_proteins_bg, 199, total_protein_list, loc=0)
## Loop for calculate hypergeometric distribution
p_val=[]
for index, row in statistics.iterrows():
    a=(row['GO'],row['go_list'],row['go_back'])
    b=hypergeom.sf(row['go_list']-1, total_proteins_bg, row['go_back'], total_protein_list, loc=0)
    p_val.append(b)
statistics['P'] = p_val

## Loop for calculate Bonferroni correction
Bonf_cor=[]
for x in statistics.P:
    Bonf_cor.append(x*multiple_testing_value)
statistics['Bonf_corr'] = Bonf_cor
statistics=statistics.sort_values(by ='P',ascending=True).reset_index(drop=True)

## Loop for calculate FDR, sorting P-value before this test
statistics=statistics[statistics.go_list > 1]
statistics['Rank']=statistics.reset_index(drop=True).index + 1
statistics=statistics.reset_index(drop=True)


# In[ ]:


#
#sys.argv[2]=float(input('[ Step 3: Choose a Value (e.g., 0.05) ]\n=====> : '))
#
FDR_val=[]
for x in statistics.Rank:
    FDR_val.append((x/statistics.count()[0])*float(sys.argv[2]))
statistics['FDR'] = FDR_val


# In[ ]:


threshold=[]
for index, row in statistics.iterrows():
    if row['P'] <= row['FDR']:
        threshold.append(row['FDR'])
if threshold == []:
    threshold=0
    #print('\n! No statistical significance threshold was found !')

#print('\nFDR Threshold P-value: ',threshold[-1])
del threshold


# In[ ]:


## Loop to add boolean value to statistically significant
# T = True if P < FDR
# F = False if P > FDR
significant=[]
for index, row in statistics.iterrows():
    if row['P'] <= row['FDR']:
        x=(str(row['FDR']))
        significant.append('T')
    else:
        if row['P'] > row['FDR']:
            y=(str(row['FDR']))
            significant.append('F')
statistics['Sig'] = significant

statistics=pd.merge(statistics,description[['GO','Term']],on='GO',how='left')

df = []
for i in list_cat.GO.drop_duplicates():
    df1 = list_cat[list_cat.GO == i]
    df.append([i,';'.join(str(e) for e in list(df1.Entry))])
proteins_by_go = DataFrame(df, columns = ['GO','entry'])

statistics=pd.merge(statistics,proteins_by_go[['GO','entry']],on='GO',how='left')
statistics.to_csv(sys.argv[3],index=None,sep='\t')

ed = pd.merge(statistics,list_cat[['GO','Entry']],on='GO',how='left')[['GO','Entry']]
ed.to_csv('data/raw_list.txt',sep='\t',index=None)
