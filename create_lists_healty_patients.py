#!/usr/bin/env python
# coding: utf-8


#Crea due liste, una per i pazienti e una per i controlli sani da passare poi a Matlab
#Da lanciare dopo unfold_and_rename [aggiungere chiamata a questo script da unfold_and_rename]

import os
import pandas as pd
import glob



df = pd.read_csv('Directories.txt',sep='\t',index_col='Index')
global prefix 
global datadir
global codedir
global windatadir
prefix = df.loc['PREFIX','value']
datadir = df.loc['Data directory','value']
codedir = df.loc['Code directory','value']
windatadir=df.loc['Windows Data Directory','value'] #se ubuntu usare path di ubuntu

with open(datadir+'subjectlist.txt') as f:
    lines=f.readlines()
f.close()
adnimerge=pd.read_excel('/home/sbombieri/DATA/UNIQUE_FILE.xlsx')



listHealty=[]
listPatients=[]
for line in lines:
    tmp=line.split('/')
    tmp=tmp[len(tmp)-1].split('_')
    tmp=tmp[len(tmp)-1] #these are now the RIDS of subjectlist
    #for every rid create the appropriate list
    tmp=tmp.strip()
    # print(tmp)
    subj=adnimerge[adnimerge['RID']==int(tmp)]['ADNIMERGE'].values[0]
    #print(tmp+" -> "+subj)
    if subj in ['CN','SMC']: #controls
        listHealty.append(line)
        #for every list create appropriate subjectlist.txt for healty
    else:
        #for every list create appropriate subjectlist.txt for patients
        listPatients.append(line)
            
#save lists:

with open(datadir+'subjectlist_PATIENTS.txt', 'w') as fp:
    for item in listPatients:
        # write each item on a new line
        fp.write("%s" % item)
    print('Done PATIENTS')
fp.close() 
ctrls=datadir+'subjectlist_CONTROLS.txt'
with open(ctrls, 'w') as fp:
    for item in listHealty:
        # write each item on a new line
        fp.write("%s" % item)
    print('Done CONTROLS')
fp.close() 


print("Done all.")
