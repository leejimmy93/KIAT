#!/bin/python 

# pip install pandas 
import pandas
import matplotlib
# import numpy as np
import matplotlib.pyplot as plt

# read file into a Pandas dataframe
from pandas import DataFrame, read_csv
f = 'https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data'
df = read_csv(f)
df=df[0:10]
print(df)
print(df.describe())

# calculate correlation
cor = df.corr()
print(cor)

# function to calculate correlation in python and remove highly correlated ones
def correlation(dataset, threshold):
    col_corr = set() # Set of all the names of deleted columns
    corr_matrix = dataset.corr()
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if corr_matrix.iloc[i, j] >= threshold:
                colname = corr_matrix.columns[i] # getting the name of column
                col_corr.add(colname)
                if colname in dataset.columns:
                    del dataset[colname] # deleting the column from the dataset

    print(dataset)

correlation(df, 0.9)

# correlation matric plot function 

def correlation_matrix(df):
    from matplotlib import pyplot as plt
    from matplotlib import cm as cm

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cmap = cm.get_cmap('jet', 30)
    cax = ax1.imshow(df.corr(), interpolation="nearest", cmap=cmap)
    ax1.grid(True)
    plt.title('Abalone Feature Correlation')
    labels=['Sex','Length','Diam','Height','Whole','Shucked','Viscera','Shell','Rings',]
    ax1.set_xticklabels(labels,fontsize=6)
    ax1.set_yticklabels(labels,fontsize=6)
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    fig.colorbar(cax, ticks=[.75,.8,.85,.90,.95,1])
    plt.show()

correlation_matrix(df)
