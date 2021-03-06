import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt # bokeh for interactive plotting
import glob
import re
import sys, os


file1='/Users/desho/Desktop/yolanda/pp_ppt_coord'
phage_all = pd.read_csv(file1, sep='\t')

file2='/Users/desho/Desktop/yolanda/RAST_phage.csv'
rast = pd.read_csv(file2, sep=',')

phage_count = phage_all[['orgId', 'scaffoldId', 'phage', 'start', 'stop']].copy()
phage_count['count'] = 0
phage_count['percentage'] = 0

#iterate through regions
#c = count; t = total
for i in range(phage_count.shape[0]):
    c = 0
    t = 0
    scaffold = phage_count.loc[i][1]
    start = (int) (phage_count.loc[i][3])
    end = (int) (phage_count.loc[i][4])

    #iterate through RAST annotations
    for j in range(rast.shape[0]):
        if rast.loc[j][1] == scaffold:
            t += 1
            anno_start = (int) (rast.loc[j][4].replace(",",""))
            anno_end = (int) (rast.loc[j][5].replace(",",""))

            #annotation in this region
            if (anno_end > start and anno_start < start) or (anno_end > end and anno_start < end) or (anno_end < end and anno_start > start):
                c += 1

    phage_count.loc[i, 'count'] = c
    if t != 0:
        phage_count.loc[i, 'percentage'] = c/t

print(phage_count)
phage_count.to_csv('phage_count.csv')