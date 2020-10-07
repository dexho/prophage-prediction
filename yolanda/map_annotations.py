import pandas as pd
import sys, os


file1='/Users/desho/Desktop/yolanda/pp_ppt_coord.txt'
pp_ppt_coord = pd.read_csv(file1, sep='\t')

file2='/Users/desho/Desktop/yolanda/RAST_phage.csv'
rast = pd.read_csv(file2, sep=',')
rast['found'] = False

#iterate through regions
for i in range(pp_ppt_coord.shape[0]):
    scaffold = pp_ppt_coord.loc[i][1]
    start = (int) (pp_ppt_coord.loc[i][3])
    end = (int) (pp_ppt_coord.loc[i][4])

    #iterate through RAST annotations
    for j in range(rast.shape[0]):
        if rast.loc[j][1] == scaffold:
            
            anno_start = (int) (rast.loc[j][4].replace(",",""))
            anno_end = (int) (rast.loc[j][5].replace(",",""))

            #annotation in this region
            if (anno_end > start and anno_start < start) or (anno_end > end and anno_start < end) or (anno_end < end and anno_start > start):
                rast.loc[j, 'found'] = True

rast.to_csv('map_annotations.csv')