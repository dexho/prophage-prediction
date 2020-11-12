"""
INPUT: CSV WITH RAST ANNOTATIONS
       CSV WITH PREDICTED PHAGE
OUTPUT: ANNOTATIONS THAT COINCIDE  
"""

import pandas as pd
import sys, os


file1='/Users/desho/Desktop/yolanda/data/pp_pptk_coord.txt'
pp_ppt_coord = pd.read_csv(file1, sep='\t')

file2='/Users/desho/Desktop/yolanda/data/rast_phage_caps.csv'
rast = pd.read_csv(file2, sep=',')
rast['found'] = False

#iterate through regions
for i in range(pp_ppt_coord.shape[0]):
    scaffold = pp_ppt_coord.loc[i][1].split(',')[0]
    start = (int) (pp_ppt_coord.loc[i][3])
    end = (int) (pp_ppt_coord.loc[i][4])

    #iterate through RAST annotations
    for j in range(rast.shape[0]):
        if rast.loc[j][2] == scaffold:
            anno_start = rast.loc[j][5]
            anno_end = rast.loc[j][6]

            #annotation in this region
            if (anno_end > start and anno_start < start) or (anno_end > end and anno_start < end) or (anno_end < end and anno_start > start):
                rast.loc[j, 'found'] = True

rast.to_csv('map_annotations_caps.csv')