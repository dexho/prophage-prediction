'''
INPUT: RAST-ANNOTATED GENOMES
INPUT: ORGID/SCAFFOLD REFERENCE
OUTPUT: A CSV CONTAINING ORGID, SCAFFOLDID, SOURCE(SERVER_RAST), PHAGE GENE NAME, START, END, LENGTH
IF THE GENE NAME CONTAINS A PHAGE-RELATED SEARCH TERM
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import pandas as pd
import sys, os

genomes = '/Users/desho/Desktop/yolanda/genomes/RAST-annotated-fit33'
os.chdir(genomes)

#create dataframe with column names
df = pd.DataFrame(columns=['orgID', 'scaffoldId', 'source', 'phage_gene', 'start', 'end', 'length'])

#phage-related search terms
search = ['phage', 'terminase', 'capsid', 'tail', 'fiber', 'sheath', 'connector', 'portal', 'tube', 'baseplate', 'plate', 'coat']
exclude = ['Macrophage', 'phage shock']

#list of dictionaries
ls_dict = []

#iterate through genomes
for gb in os.listdir(genomes):

    #new dicitonary => 1 row of dataframe for phage annotation
    phage_dict = {}

    seq_record = SeqIO.parse(gb, "genbank")
    orgID = gb.split('.')[0]
    
    # add entry in dataframe if it contains one of the search terms
    for record in seq_record:
        scaffoldID = record.id[len(orgID)+1:]
        scaffoldID = scaffoldID.split('.')[0]

        for feature in record.features:
            if feature.type == 'CDS':
                try:
                    desc = feature.qualifiers['product'][0]
                    for s in search:
                        if s in desc and 'Macrophage' not in desc and 'phage shock' not in desc:
                            phage_dict['orgID'] = orgID
                            phage_dict['scaffoldID'] = scaffoldID
                            phage_dict['source'] = 'server_RAST'
                            phage_dict['phage_gene'] = desc
                            phage_dict['start'] = feature.location.start
                            phage_dict['end'] = feature.location.end
                            phage_dict['len'] = feature.location.end - feature.location.start
                            continue
                except:
                    continue

            ls_dict.extend([phage_dict])

df = pd.DataFrame(ls_dict)
df.to_csv('rast_phage.csv')