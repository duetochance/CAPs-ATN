import os
#Add prefix for TbCAPs
path = r"/mnt/f/Tirocinio_HDD/Data/ADNI/DATA_TO_BE_RENAMED"

directory_list = os.listdir(path)



listempty=[]
for filename in directory_list:

   #print(filename)
   if len(os.listdir(path+"/"+filename))==0:
    listempty.append(filename)
    

with open('EMPTY.txt', 'w') as f:
    for line in listempty:
        f.write(line)
        f.write('\n')

print("__")