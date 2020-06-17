import sys
import csv
import itertools
import math
import os
import numpy as np
import pandas as pd

# utilities used for encoding sets in ASP (binary encoding)
def BinToInt(bin_str):
    return int(bin_str, 2)

def IntToBin(int_str):
    return str(bin(int_str))[2:]

def IntToSubset(int_str):
    bin_string = IntToBin(int_str)[::-1]
    bin_list = []
    n=0
    for c in bin_string:
        if int(c) == 1:
            bin_list.append(n)
        n += 1
    return bin_list

def SubsetToInt(subset):
    string = ''
    for n in range(max(subset)+1):
        if n in subset:
            string += '1'
        else:
            string += '0'
    bin_str=string[::-1]
    return BinToInt(bin_str)

def readData():
    input_file = sys.argv[1]
    Y = int(sys.argv[2])
    I = int(sys.argv[3])
    #S = int(sys.argv[4])
    currentDir = sys.argv[5]

    input_data = pd.read_csv(input_file)

    # Read metadata file (for numbers of variables and blinding)
    is_split = False
    li = input_file.rsplit('real-', 1)
    if len(li) > 1:
        is_split = True
    meta_file = 'meta-'.join(li)
    li = meta_file.rsplit('nan-', 1)
    if len(li) > 1:
        is_split = True
    meta_file = 'meta-'.join(li)
    if is_split and os.path.isfile(meta_file):
        metadata_present = True
        metadata = pd.read_csv(meta_file)
        num_sys_vars = metadata.loc[0,'num_sys_vars']
        num_ints_vars = metadata.loc[0,'num_ints_vars']
        IBlind = metadata.loc[0,'IBlind']
        IBlind_visible_value = metadata.loc[0,'IBlind_visible_value']
    else:
        metadata_present = False

    if not metadata_present or input_data.columns.size == num_sys_vars + 1:
        # Context variables not provided in data files: construct a diagonal design
        regime_column = input_data.columns[input_data.columns.size-1]

        # Construct train and test data
        unknown_rows = input_data[regime_column].isin([I+1])

        # Create the experimental design with diagonal matrix and then ignore the regime:
        #dummy_data = pd.get_dummies(input_data[regime_column],drop_first=True) # drop_first not on cluster
        dummy_data_withfirst = pd.get_dummies(input_data[regime_column])
        data = pd.concat([input_data, dummy_data_withfirst.iloc[:,1:]], axis=1)
    elif input_data.columns.size == num_sys_vars + 1 + num_ints_vars:
        # Context variables are provided in data files

        # The test set consists of those rows where context variable IBlind's value is *not* equal to
        # IBlind_visible_value
        IBlind_column = input_data.columns[num_sys_vars+IBlind] # not +1 because Python uses 0-based indexing
        unknown_rows = ~(input_data[IBlind_column].isin([IBlind_visible_value]))

        # No need to add any columns
        data = input_data
    else:
        print "Error: expected either "+str(num_sys_vars + 1)+" (sys+R) or "+str(num_sys_vars + 1 + num_ints_vars)+" (sys+R+context) columns; got "+str(input_data.columns.size)
        sys.exit()

    #print "training points: ", y_train_data.shape[0]
    #print "test points: ", y_test_data.shape[0]

    data_unknown_rows = data.loc[unknown_rows]
    data_known_rows = data.loc[~unknown_rows]

    X_test_data = data_unknown_rows.copy()
    y_test_data = data_unknown_rows[data_known_rows.columns[Y-1]].copy()
    X_train_data = data_known_rows.copy()
    y_train_data = data_known_rows[data_known_rows.columns[Y-1]].copy()

    return X_test_data, y_test_data, X_train_data, y_train_data

def subsetArg():
    Y = int(sys.argv[2])
    S = int(sys.argv[4])

    subset = IntToSubset(S)

    # Double check that Y is not in the subset S.
    if Y-1 in subset:
        print "Error: Y=" + str(Y) + " is in the subset S=" + str(S) + " (" + str(subset) + ")."
        sys.exit()

    return subset

def outputFilename(tag):
    input_file = sys.argv[1]
    Y = int(sys.argv[2])
    I = int(sys.argv[3])
    S = int(sys.argv[4])
    currentDir = sys.argv[5]

    i = input_file.split("-")[1].split(".")[0]
    return os.path.join(currentDir, tag + "-" + i + "-" + str(Y) + "-" + str(I) + "-" + str(S) + ".csv")

def numBootstraps(defaultNum):
    num = defaultNum
    if len(sys.argv) > 6:
        num = int(sys.argv[6])
    return num
