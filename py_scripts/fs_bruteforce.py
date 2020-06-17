import sys
import csv
import itertools
import math
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from helper import *

# Find the best possible subset of features (w.r.t oob error of a trained random forest) by trying all possible subsets.

# Input parameters in sys.argv (in order):
# - input_file, filename has to be in form "anytext-number.csv", e.g. "/zfs/datastore/smaglia1/NIPS2017/1/sim-1.csv"
# -- in the input file, we don't have intervention variables and the regime variable is supposed to be the last column.
# - Y: the variable we are trying to predict, e.g. 2, which in terms of columns is actually Y-1 = 1 because we start from 0
# - I: the transfer learning dataset (note R here starts by 0 instead than 1 as in R, so +1), e.g. I=1 means that we will remove dataset with R=2
# - S: integer representation of the subset we want to consider (although here all variables are indexed from 0), e.g. 17 = {X1, I1}
# - currentDir: output directory

# Output: set of integers from best subset in "fs-" + number_from_input_file + "-" + str(Y) + "-" + str(I) + "-" + str(S) + ".csv"

def Subsets(subset):
    combs = []

    for i in xrange(1, len(subset)+1):
        els = [list(x) for x in itertools.combinations(subset, i)]
        combs.append(els)
    return combs

def RegressionAllSubsets(X,y,subset,n_estimators):

    allSubsets = Subsets(subset)
   
    length=len(subset)
 
    OoBErrorList = []
    setList = []

    for i in xrange(0, length):
        for j, iSet in enumerate(allSubsets[i]):
            model = RandomForestRegressor(n_estimators=n_estimators, oob_score=True, n_jobs=-1)
            model.fit(X[X.columns[iSet]],y)
            OoBErrorList.append(1-model.oob_score_)
            setList.append(allSubsets[i][j])

    SortedOoBError = sorted(OoBErrorList)
    SortedOoBErrorIndex = [i[0] for i in sorted(enumerate(OoBErrorList), key= lambda x:x[1])]
    
    outputOoBErrorList = []
    outputsetList = [] + setList
    
    for n in SortedOoBErrorIndex:
        outputOoBErrorList.append([OoBErrorList[n],n+1])
    
    columns = ["OoBErr","length","Set"]

    df = pd.DataFrame(columns=columns)

    for i in outputOoBErrorList:
        row_i = [] + i
        row_i.append(outputsetList[row_i[1]-1])
        df=df.append(dict(zip(columns,row_i)), ignore_index=True)

    return df

X_test_data, y_test_data, X_train_data, y_train_data = readData()
subset = subsetArg()
n_estimators = numBootstraps(1000)

#Here y_train column should not be in subset
D = RegressionAllSubsets(X_train_data,y_train_data,subset,n_estimators)

set_list = D["Set"].tolist()

int_list = map(SubsetToInt, set_list)

output_file = output_file = outputFilename("fs")
with open(output_file, 'wb') as output:
    wr = csv.writer(output)
    wr.writerow(int_list)


