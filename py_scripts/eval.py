import sys
import csv
import itertools
import math
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from helper import *

# Fit a random forest on a subset of columns of the training data, then write the mean squared error in a file.

# Note: this piece of code creates intervention variables with the diagonal design, starting from the regime. 
# Note2: here the variables are indexed from 0, while in ASP and R they start at 1.

# Input parameters in sys.argv (in order):
# - input_file, filename has to be in form "anytext-number.csv", e.g. "/zfs/datastore/smaglia1/NIPS2017/1/sim-1.csv"
# -- in the input file, we don't have intervention variables and the regime variable is supposed to be the last column.
# - Y: the variable we are trying to predict, e.g. 2, which in terms of columns is actually Y-1 = 1 because we start from 0
# - I: the transfer learning dataset (note R here starts by 0 instead than 1, so +1), e.g. I=1 means that we will remove dataset with R=2
# - S: integer representation of the subset we want to consider (although here all variables are indexed from 0), e.g. 17 = {X1, I1}
# - currentDir: output directory

# Output: mserror in "eval-" + number_from_input_file + "-" + str(Y) + "-" + str(I) + "-" + str(S) + ".csv"

X_test_data, y_test_data, X_train_data, y_train_data = readData()
subset = subsetArg()
n_estimators = numBootstraps(1000)

# Here y_train column should not be in subset

# Train the random forest model:
#model = RandomForestRegressor(n_estimators=1000, oob_score=True, n_jobs=-1, random_state=531,max_features=0.2, min_samples_leaf=1)
model = RandomForestRegressor(n_estimators=n_estimators, oob_score=True, n_jobs=-1)
model.fit(X_train_data[X_train_data.columns[subset]], y_train_data)

msError = [str(mean_squared_error(y_test_data, model.predict(X_test_data[X_test_data.columns[subset]])))]

output_file = outputFilename("eval")
with open(output_file, 'wb') as output:
    wr = csv.writer(output)
    wr.writerow(msError)


