#!/usr/bin/env python
# coding: utf-8


#launched after unfold_and_rename.py

import csv
import os
import pandas as pd
import glob
import shutil
import itertools

def getDirectories():
        df = pd.read_csv('/home/sbombieri/CODE/Code/Directories.txt',sep='\t',index_col='Index')
        global prefix 
        global datadir
        global codedir
        global windatadir
        prefix = df.loc['PREFIX','value']
        datadir = df.loc['Data directory','value']
        codedir = df.loc['Code directory','value']
        windatadir=df.loc['Windows Data Directory','value']


if __name__ == '__main__':
    getDirectories()
    print(" DIRECTORIES ARE: \n")
    print(datadir+"\n"+codedir+"\n")
    os.chdir(datadir)
    o=1 #SET TO 0 ONLY WHEN TESTING
    if o==0:
        print("OOOOK")
        exit()

    data = pd.read_csv(datadir+'subjectlist.txt', sep=" ", header=None)
    oldpwd = os.getcwd()
    print(oldpwd)

    print(data)
    lstring=len(data.loc[0].values[0].split(sep='/'))
    for sub in range(0,len(data)):
        tmp=data.loc[sub].values[0].split(sep='/')[lstring-1]
        print("-> "+tmp)
        
        os.chdir(tmp+'/')
        print(os.getcwd())
        os.system('cp mov_warped_func.txt mov_warped_func_BACKUP_COPY.txt')
        shutil.copyfile('mov_warped_func.txt','mov_warped_func_BACKUP_COPY.txt')
    # print(os.listdir())
        vols=glob.glob("vol*.nii")
        #print(glob.glob("*.nii"))
        print('Prev. length:')
        print(len(vols))
        #until reaching right length pop last element, setp volumes
        while len(vols)>192:
            vol_to_erase=vols.pop(len(vols)-1)
            #print(vol_to_erase)
            os.remove(vol_to_erase)

        #until reaching right length pop last element, setp mov_warped_func
        df=pd.read_csv("mov_warped_func.txt",header=None)
        
        print("Removing from mov_warped_func.txt")
        if len(df)>192:
            while len(df)>192:
                df=df.drop([len(df)-1])
            df.to_csv(r'mov_warped_func.txt', header=False, index=False, sep=' ', mode='w',quoting=csv.QUOTE_NONE,escapechar=' ')
        print('Aft. length:')    
        print(len(vols))
        os.chdir(oldpwd)
