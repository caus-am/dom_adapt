import os
import sys
import random
import numpy as np
import pandas as pd

#number of output files
#n=100

data_folder = str(sys.argv[1])
n = int(sys.argv[2])

if not os.path.exists(data_folder):
    os.makedirs(data_folder)

data_file_mouse = "../mouse.csv"
data = pd.read_csv(data_file_mouse)

# Select rows
#data[data["geno"].str.match('1796_1')]

gene_list = ['1797_1', '1796_1', '1798_1', '727_1', '1550_1', '1799_1',
       '3157_1', '3621_1', '3803_1', '3805_1', '3887_1', '4045_1', '4047_1']

target_list = ['IMPC_HEM_001_001','IMPC_HEM_002_001', 'IMPC_HEM_003_001', 'IMPC_HEM_004_001',
       'IMPC_HEM_005_001', 'IMPC_HEM_006_001', 'IMPC_HEM_007_001',
       'IMPC_HEM_008_001', 'IMPC_HEM_019_001', 'IMPC_HEM_027_001',
       'IMPC_HEM_029_001', 'IMPC_HEM_030_001', 'IMPC_HEM_031_001',
       'IMPC_HEM_032_001', 'IMPC_HEM_033_001', 'IMPC_HEM_034_001',
       'IMPC_HEM_035_001', 'IMPC_HEM_036_001', 'IMPC_HEM_037_001',

       'IMPC_HEM_038_001', 'IMPC_HEM_039_001', 'IMPC_HEM_040_001']

#Remove deterministic relations 
det_list = ['IMPC_HEM_001_001', 'IMPC_HEM_029_001','IMPC_HEM_031_001','IMPC_HEM_033_001','IMPC_HEM_035_001','IMPC_HEM_038_001','IMPC_HEM_040_001']

column_list = [item for item in target_list if item not in det_list]

def CreateDataSet(data_folder, data_file_number, v1, v2, target, I1, I2):
    dat0 = data[data['geno'].str.match('0')][[v1,v2,target,'geno']]
    dat0.geno[dat0.geno == '0'] = 1
    dat1 = data[data['geno'].str.match(I1)][[v1,v2,target,'geno']]
    dat1.geno[dat1.geno == I1] = 2
    dat1nan = data[data['geno'].str.match(I1)][[v1,v2,target,'geno']]
    dat1nan.geno[dat1nan.geno == I1] = 2
    dat1nan[target] = np.nan
    dat2 = data[data['geno'].str.match(I2)][[v1,v2,target,'geno']]
    dat2.geno[dat2.geno == I2] = 3

    data_final = pd.concat([dat0, dat1, dat2])
    data_final_nan = pd.concat([dat0, dat1nan, dat2])

    data_final.columns = ['v1', 'v2', 'v3', 'experiment']
    data_final_nan.columns = ['v1', 'v2', 'v3', 'experiment']

    data_final.to_csv(data_folder + "/real-" + str(data_file_number) + ".csv", index = False)
    data_final_nan.to_csv(data_folder + "/nan-" + str(data_file_number) + ".csv", index = False)

    text_string = data_folder + "/real-" + str(data_file_number) + ": v1=" + v1 + ", v2=" + v2 + ", v3=" + target + " (target), I_1=" + I1 + ", I_2=" + I2 + "\n" 
    with open(data_folder + "/real_data.txt", "a") as myfile:
        myfile.write(text_string)

for i in range(1,n+1):
    v1, v2, target = random.sample(column_list, 3)
    I1, I2 = random.sample(gene_list,2)
    CreateDataSet(data_folder, i, v1, v2, target, I1, I2)

