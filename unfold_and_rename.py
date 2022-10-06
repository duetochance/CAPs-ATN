"""
Iterate through subjects directories and extract volumes following an unzip for make 
data suitable to Matlab SPM12 functions. Creating subjectlist.txt filetext for CAP_TB 
script.

Just for Windows environment for now.


Da lanciare 2 volte sui sani e poi sui malati
"""



from asyncio import ALL_COMPLETED
import code
from multiprocessing.connection import wait
import os
import multiprocessing as mp
import concurrent.futures
import pandas as pd
import time

###### Directories
#datadir='/mnt/c/Users/sambo/Desktop/UNIVERSITAS/OneDrive/UNIVERSITAS/MEDICAL BIOINFORMATICS/TIROCINIO_E_TESI/'
#codedir=os.getcwd()

######

def remove_unnecessary():
    print(">> Removing unnecessary files ... ")
    #datadir= '/mnt/f/Tirocinio_HDD/Data/ADNI/' #remove in .py execution
    os.chdir(datadir)
    f = open(datadir+'subjectlist.txt')
    lines=f.readlines()
    for line in lines:
        l=line.split('/')
        l=l[len(l)-1]
        l=l.rstrip()
        os.chdir(datadir+l)
        print("->  "+datadir+l)
        if len(os.listdir()) == 0: 
            print("empty")
        else:
            if os.path.exists("mov_warped_func_BACKUP_COPY.txt"):
                os.remove('mov_warped_func_BACKUP_COPY.txt')
            if os.path.exists("mov_warped_func.nii.gz"):
                os.remove('mov_warped_func.nii.gz')
            if os.path.isdir('Backup'):
                os.chdir(datadir+l+'/Backup')
                os.remove('warped_func.nii.gz')
                os.chdir(datadir+l)
                os.rmdir('Backup')
        os.chdir(datadir)
    print('> Removing complete.')




def getMovementFile(subj,dir_tmp):
    #print("############")
    outname='mov_'+subj.split('.')[0]
    mcflirt='mcflirt -in '+subj+' -out '+outname+' -plots -spline_final'
    print(mcflirt)
    os.system(mcflirt)

    #renaming to .txt
    rename='mv '+outname+'.par'+' '+outname+'.txt'
    os.system(rename)
    #print("############")
    return "OK"
    



def unzip_and_extractvol(subj):
    #print("###########################")
    tmpdir=datadir+subj+"/"
    os.chdir(tmpdir)
    l=os.listdir()[0]
    print("------- "+str(os.getcwd())+" ------- ")
    ##create copy of original
    #newname="Copy_"+l
    #copy='mkdir Backup && cp '+l+' Backup'
    #print(copy)
    #os.system(copy)

   #Test for integrity:
    os.system("if gunzip -t warped_func.nii.gz ;then echo 'file is ok'; else      echo 'file is corrupt';pwd >> /home/sbombieri/DATA/LOGFILE.txt; fi ")
    print("################# UNZIPPING ")
    print(subj)
    print(tmpdir)
    print("#################")
    try:
        print(getMovementFile(l,tmpdir))
    except:
        print("-> ERROR get movement file")
        log_output(subj)

    ##separate from 3D to 4D
    output_basename=" volsubj"
    fslsplit='fslsplit '+l+output_basename
    #print(fslsplit)
    try:
        os.system(fslsplit)

        ##unzip for CAP_TB
        gunzip='for f in ./vol*; do gunzip  $f -f; done'
        #print(gunzip)
        os.system(gunzip)
        print('\n------- Unzipped file in '+subj+' named '+l+' --------\n')
    except:
        log_output(subj)
    os.chdir(datadir)
    

# from bash to windows format
def refinelist(subjlist):
     print("Creating *subjectlist* for inputting correct data in Matlab..\n")
    # dir="F:\Tirocinio_HDD\Data\ADNI\\"
     dir=datadir
     filePath=datadir+"subjectlist.txt"
     if os.path.exists(filePath):
        os.remove(filePath)
        print("Removed previous file")
     
     print("Writing ...")
     with open("subjectlist.txt", 'a') as the_file:
         with open (subjlist,'r') as file:
             for line in file:
                 #print("# "+line)
                 the_file.write(windatadir+line)
    
     print("Converted correctly, created file subjectlist in "+filePath)
     the_file.close()
     file.close()
     os.remove(subjlist)

def log_output(subj):
    print("------------------>Writing log to :"+datadir+"LOG_ERRORS_UNFOLD_N_RENAME.txt")

    f=open(datadir+"LOG_ERRORS_UNFOLD_N_RENAME.txt",'w')
    f.write(str("SUB ID -> "+subj+"\n"))
    f.close()
                
def getDirectories():
        df = pd.read_csv('Directories.txt',sep='\t',index_col='Index')
        global prefix 
        global datadir
        global codedir
        global windatadir
        prefix = df.loc['PREFIX','value']
        datadir = df.loc['Data directory','value']
        codedir = df.loc['Code directory','value']
        windatadir=df.loc['Windows Data Directory','value'] #se ubuntu usare path di ubuntu
        

def setTerrain():
    prev=os.getcwd()
    os.chdir(datadir)
    #obj: remove all previous volumes and files, reset to only warped.nii.gz
    l=os.listdir()
    for f in l:
        if 'sub-' in f:
            vol="rm "+datadir+f+"/vol*"
            mov="rm "+datadir+f+"/mov*"
            os.system(vol)
            os.system(mov)
    print("Default settings are OK.")
    os.chdir(prev)
    time.sleep(2.5)


if __name__ == '__main__':
    print("\n*             VERSION 2.2.1             *\n")
    #-MAIN-#
    getDirectories()
    print("## Prefix = "+prefix) #prefix of subject folder
    print("## Datadir= "+datadir) #where  to save and manage data
    print("## Codedir = "+codedir) #where this code is being launched
    print("Setting basic things...")
    os.system("rm "+datadir+"LOGFILE.txt")
    setTerrain()
    #os.system("touch "+datadir+"LOGFILE.txt")
    #create subjectlist in datadir
    os.chdir(datadir)
    command=" ls | grep "+prefix+" > subjectlist_tmp.txt"
    os.system(command)
    

    #read subjectlist 
    f=open("subjectlist_tmp.txt","r")
    subjs=f.read().split()
    f.close()
    print("LIST:")
    for s in subjs:
        print("#"+s)

    #another refinement step:
    data = pd.read_csv('EMPTY.txt', sep=" ", header=None)
    sublist=pd.read_csv("subjectlist_tmp.txt",header=None)
    #delete momentarely void folders
    for i in range(0,len(data)):
        
        sublist.drop(sublist[sublist[0]==data.loc[i].values[0]].index,inplace=True)
    sublist.to_csv(r'subjectlist_tmp.txt', header=None, index=None, sep=' ', mode='w')
    #i=0
    with concurrent.futures.ProcessPoolExecutor() as executor:
     for (subj) in subjs:
            #print("Getting movement for "+subj)
            #fs = [executor.submit(getMovementFile,subj)]
            #concurrent.futures.wait(fs,return_when=ALL_COMPLETED)
            #print("## "+str(fs[0].result))
            #i+=1
            print("Extracting volumes from "+subj)
            executor.submit(unzip_and_extractvol,subj)
            #unzip_and_extractvol(subj)
            print("##-----		DONE for subj: "+subj+"				-----##")
    #with concurrent.futures.ProcessPoolExecutor() as executor:
       # for (subj) in subjs:
         #   print("Extracting volumes from "+subj)
         #  executor.submit(unzip_and_extractvol,subj)
    

        
       
        

    refinelist("subjectlist_tmp.txt")
   # os.system("python3 "+codedir+"correctSubjlist.py")
    print('Done unfolding and renamig, now refinig the length of extracted volumes...')
    #togli volumi in più
    command="python3 "+codedir+"Togli_volumi.py"
    os.system(command)
    #crea le due liste per lo script matlab
    command="python3 "+codedir+"create_lists_healty_patients.py"
    os.system(command) 
     
    remove_unnecessary()
    print('Done!\nNow start Matlab for CAPs generation ;) ')
