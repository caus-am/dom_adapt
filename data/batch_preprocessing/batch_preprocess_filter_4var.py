import os
import sys
import random
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


#number of output files
#data_folder = "./"
#n=100

data_folder = str(sys.argv[1])
n = int(sys.argv[2])

#ttest_threshold = 1e-16 
ttest_threshold = float(sys.argv[3])

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

def CreateDataSet(data_folder, n, column_list, gene_list, ttest_threshold):
    i = 0
    while (i <= n):
        v1, v2, v3, target = random.sample(column_list, 4)
        I1, I2, I3 = random.sample(gene_list,3)

        dat0 = data[data['geno'].str.match('0')][[v1,v2,v3,target,'geno']]
        dat0.geno[dat0.geno == '0'] = 1
        dat1 = data[data['geno'].str.match(I1)][[v1,v2,v3,target,'geno']]
        dat1.geno[dat1.geno == I1] = 2
        dat1nan = data[data['geno'].str.match(I1)][[v1,v2,v3,target,'geno']]
        dat1nan.geno[dat1nan.geno == I1] = 2
        dat1nan[target] = np.nan
        dat2 = data[data['geno'].str.match(I2)][[v1,v2,v3,target,'geno']]
        dat2.geno[dat2.geno == I2] = 3
        dat3 = data[data['geno'].str.match(I3)][[v1,v2,v3,target,'geno']]
        dat3.geno[dat3.geno == I3] = 4

        t = False

        for j in range(0,4):
            t1, p1 = ttest_ind(dat0.ix[:,j],dat1.ix[:,j], equal_var=False) 
            t2, p2 = ttest_ind(dat0.ix[:,j],dat2.ix[:,j], equal_var=False)
            t3, p3 = ttest_ind(dat0.ix[:,j],dat3.ix[:,j], equal_var=False)
            if p1<ttest_threshold or p2<ttest_threshold or p3<ttest_threshold:
                t = True
                break
#            else:
#                print str(p1) + " and " + str(p2) + " i=" + str(j) + " v1=" + v1 + " v2=" + v2 + " target=" + target + " I1=" + I1 + " I2=" + I2 + "\n"
               
        if t == True:
            data_final = pd.concat([dat0, dat1, dat2, dat3])
            data_final_nan = pd.concat([dat0, dat1nan, dat2, dat3])

            data_final.columns = ['v1', 'v2', 'v3', 'v4', 'experiment']
            data_final_nan.columns = ['v1', 'v2', 'v3', 'v4', 'experiment']

            data_final.to_csv(data_folder + "/real-" + str(i) + ".csv", index = False)
            data_final_nan.to_csv(data_folder + "/nan-" + str(i) + ".csv", index = False)

            text_string = data_folder + "/real-" + str(i) + ": v1=" + v1 + ", v2=" + v2 + ", v3=" + v3 + ", v4=" + target + " (target), I_1=" + I1 + ", I_2=" + I2 + ", I_3=" + I3 + "\n" 
            with open(data_folder + "/real_data.txt", "a") as myfile:
                myfile.write(text_string)
            i+=1

CreateDataSet(data_folder,n,column_list,gene_list,ttest_threshold)
