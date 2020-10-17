'''
INPUT: RAST-ANNOTATED GENOMES
INPUT: ORGID/SCAFFOLD REFERENCE
OUTPUT: A CSV CONTAINING ORGID, SCAFFOLDID, SOURCE(SERVER_RAST), PHAGE GENE NAME, START, END, LENGTH
IF THE GENE NAME CONTAINS A PHAGE-RELATED SEARCH TERM
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import sys, os

genomes = '/Users/desho/Desktop/yolanda/genomes/RAST-annotated-fit33'
os.chdir(genomes)

#create dataframe with column names
df = pd.DataFrame(columns=['orID', 'scaffoldId', 'source', 'phage_gene', 'start', 'end', 'length'])

#phage-related search terms
search = ['phage', 'terminase', 'capsid', 'repressor', 'tail', 'fiber', 'sheath', 'connector', 'portal', 'tube', 'baseplate', 'plate', 'coat']

#iterate through genomes
for gb in os.listdir(genomes):
    seq_record = SeqIO.parse(gb, "genbank")

    orgID = gb.split('.')[0]
    
    # add entry in dataframe if it contains one of the search terms
    for record in seq_record:
        scaffoldID = record.id[len(orgID)+1:]
        print(scaffoldID)
        for feature in record.features:
            phage = False
            if feature.type == 'CDS':
                try:
                    desc = feature.qualifiers['product'][0]
                    for s in search:
                        if s in desc:
                            phage = True
                    if phage:
                        df.iloc[0] = feature
                        
                except:
                    continue





# df.to_csv('map_annotations.csv')