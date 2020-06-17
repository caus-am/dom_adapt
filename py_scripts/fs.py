import sys
import csv
import itertools
import math
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from helper import *

# Find a good subset of features with stepwise forward regression, i.e. add gredily features based on the smallest oob error of the random forest trained with them.

# Input parameters in sys.argv (in order):
# - input_file, filename has to be in form "anytext-number.csv", e.g. "/zfs/datastore/smaglia1/NIPS2017/1/sim-1.csv"
# -- in the input file, we don't have intervention variables and the regime variable is supposed to be the last column.
# - Y: the variable we are trying to predict, e.g. 2, which in terms of columns is actually Y-1 = 1 because we start from 0
# - I: the transfer learning dataset (note R here starts by 0 instead than 1 as in R, so +1), e.g. I=1 means that we will remove dataset with R=2
# - S: integer representation of the subset we want to consider (although here all variables are indexed from 0), e.g. 17 = {X1, I1}
# - currentDir: output directory

# Output: set of integers from best subset in "fs-" + number_from_input_file + "-" + str(Y) + "-" + str(I) + "-" + str(S) + ".csv"

def ForwardStepwiseRegression(X,y,subset,n_estimators):

#    print "Start forward stepwise regression using random forests"

    OoBErrorList = []

    chosenSetList = []

    indexSet = set(subset)

    for i in indexSet:
        chosenSet = set(chosenSetList)
        setTrySet = indexSet - chosenSet
        setTryList =  [j for j in setTrySet]
        errorList = []
        #Useless? setTemp = []
        for iSetTry in setTryList:
            SetTemp = [] + chosenSetList 
            SetTemp.append(iSetTry)
            #model = RandomForestRegressor(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=531,max_features=0.2, min_samples_leaf=1)
            model = RandomForestRegressor(n_estimators=n_estimators, oob_score=True, n_jobs=-1)
            model.fit(X[X.columns[SetTemp]],y)
            errorList.append(1-model.oob_score_)
            #Useless? SetTemp = []
        iBest = np.argmin(errorList)
        OoBErrorList.append(errorList[iBest])
        chosenSetList.append(setTryList[iBest])
#        print "\r{0}%".format(math.ceil(float(len(chosenSet)+1)/len(indexSet)*100))
  

    SortedOoBError = sorted(OoBErrorList)
    SortedOoBErrorIndex = [i[0] for i in sorted(enumerate(OoBErrorList), key= lambda x:x[1])]

    outputOoBErrorList = []
    outputsetList = [] + chosenSetList

    for n in SortedOoBErrorIndex:
        outputOoBErrorList.append([OoBErrorList[n],n+1])

    columns = ["OoBErr","length","Set"]

    df = pd.DataFrame(columns=columns)

    for i in outputOoBErrorList:
        row_i = [] + i
        row_i.append(outputsetList[:row_i[1]])
        df=df.append(dict(zip(columns,row_i)), ignore_index=True)

    return df


X_test_data, y_test_data, X_train_data, y_train_data = readData()
subset = subsetArg()
n_estimators = numBootstraps(1000)

# Here y_train column should not be in subset
D = ForwardStepwiseRegression(X_train_data,y_train_data,subset,n_estimators)

set_list = D["Set"].tolist()

int_list = map(SubsetToInt, set_list)

output_file = outputFilename("fs")
with open(output_file, 'wb') as output:
    wr = csv.writer(output)
    wr.writerow(int_list)


