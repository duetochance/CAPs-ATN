import os
import pandas as pd
import pandasql as ps

#-- Data names
# directory
dir="/mnt/f/Tirocinio_HDD/Data/ADNI/"
#files
adni_data_summary=dir+"ADNI_data_summary_02-05-2022.xlsx"
UPENNBIOMK1=dir+"UPENNBIOMK9_04_19_17.csv"
UPENNBIOMK2=dir+"UPENNBIOMK10_07_29_19.csv"
UPENNBIOMK3=dir+"UPENNBIOMK12_01_04_21.csv"

#columns
""" hippo
enthorinal
amigdala
FDG-PET (metabolismo): stesse regioni
Tau-PET: stesse regioni
Amy-PET: stesse regioni + globale
info dal CSF:
pTau
Tau
Ab42 """


#--


adni=pd.read_excel(adni_data_summary)
upenn1=pd.read_csv(UPENNBIOMK1)
upenn2=pd.read_csv(UPENNBIOMK2)
upenn3=pd.read_csv(UPENNBIOMK3)
#df[df.columns[df.columns.isin(['alcohol','hue','NON-EXISTANT COLUMN'])]]

""" print(adni.head(3))
print("----")
print(upenn1.head(3))
print("----")
print(upenn2.head(3))
print("----")
print(upenn3.head(3))
 """

outlist=[] #empty, to be filled and passed to final pandas dataframe

adni=adni[adni.columns[adni.columns.isin(['ID','RID','ADNIMERGE','ACQUISITION DATE','AGE','SEX'])]]
upenn1=upenn1[upenn1.columns[upenn1.columns.isin(['RID','EXAMDATE','RUNDATE','PTAU','TAU'])]]
upenn2=upenn2[upenn2.columns[upenn2.columns.isin(['RID','DRAWDATE','RUNDATE','ABETA42','PTAU','TAU'])]]
upenn3=upenn3[upenn3.columns[upenn3.columns.isin(['RID','EXAMDATE','RUNDATE','PTAU','TAU','AB4240'])]]

adni.rename(columns={'ACQUISITION DATE':'ACQUISITIONDATE'},inplace=True)
#select only rows , no the commets under the last row
adni=adni.dropna(subset=['RID','ACQUISITIONDATE'])
adni=adni.astype({'RID':'int'})
q1="""SELECT * FROM adni,upenn1 WHERE adni.RID=upenn1.RID """



