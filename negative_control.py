import pandas as pd
import math

file1='/Users/desho/Desktop/yolanda/pp_ppt_coord.txt'
pp_ppt_coord = pd.read_csv(file1, sep='\t')

file2='/Users/desho/Desktop/yolanda/RAST_phage.csv'
rast = pd.read_csv(file2, sep=',')

ls_control = []

scaffoldID = ''
#iterate through RAST orgID/scaffold pairs
#cast everything to type str so no comparing str to int or floats

for i in range(rast.shape[0]):
    #skip if scaffoldID is identical to the next one
    if (str)(rast.loc[i][1]) == scaffoldID:
        continue
    if type(rast.loc[i][0]) is str:
        orgID = (str)(rast.loc[i][0])
    scaffoldID = (str)(rast.loc[i][1])

    found = False
    #iterate through pp_ppt_coord for orgID/scaffoldId
    for j in range(pp_ppt_coord.shape[0]):
        if str(pp_ppt_coord.loc[j][1]) == scaffoldID and str(pp_ppt_coord.loc[j][0]) == orgID:
            found = True
            continue

    if not found:
        ls_control.extend([{'orgId': orgID, 'scaffoldId': scaffoldID}])

control = pd.DataFrame(ls_control)
control.to_csv('control.csv')
    
