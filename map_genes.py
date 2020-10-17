print("""
    map_genes: retrieves genes from df_search_genes based on regions in df_query_regions. 
    Default is to first group by columns orgId and scaffoldId. Groupby is optional.
    Default is to expand query region by 0 kb. Inclusion of flanking region is optional.
    Arguments: (df_query_regions_regions, df_search_genes_genes, query_cols=['start', 'stop'],
        search_cols = ['begin', 'end'], col_groupby=['orgId','scaffoldId'],
        flanking_len=0)
    """)

import pandas as pd

def map_genes(df_query_regions, df_search_genes, query_cols=['start', 'stop'], search_cols = ['begin', 'end'], 
               col_groupby=['orgId','scaffoldId'], flanking_len = 0):
    df_out_genes = pd.DataFrame()      
    # With groupby
    if len(col_groupby)!=0:
        df_query_regions = df_query_regions.groupby(col_groupby)
        df_search_genes = df_search_genes.groupby(col_groupby)
        
        # convert groupby df into dictionaries
        df_regions_grouped = dict(list(df_query_regions))
        df_genes_grouped = dict(list(df_search_genes))
        
        # groupby keys that are in both dfs
        intersect_keys = set(df_regions_grouped.keys()).intersection(set(df_genes_grouped.keys()))

        for key in df_regions_grouped:
            # if key is present in both dfs. {key} to convert tuple key into a set of single tuple for comparison.
            if {key}.issubset(intersect_keys):
                # assign positions for query and search
                region_start = df_regions_grouped[key][query_cols[0]] - flanking_len
                region_stop = df_regions_grouped[key][query_cols[1]] + flanking_len
                gene_start = df_genes_grouped[key][search_cols[0]]
                gene_stop = df_genes_grouped[key][search_cols[1]]

                # check that start and stop are sequential in phage regions (genes aren't)
                if (region_start >= region_stop).any():
                    return print('Error: {} query start position is greater than stop position!') \
                                 .format(key)
                
                # compares genes against region positions individually:
                for i in range(len(region_start)):
                    # if start or stop of gene is within query region
                    compare_overlap = (gene_start >= region_start.iloc[i]) & (gene_start <= region_stop.iloc[i]) \
                                   | (gene_stop >= region_start.iloc[i]) & (gene_stop <= region_stop.iloc[i])
                    # genes part of current phage region
                    df_genes_single_region = df_genes_grouped[key][compare_overlap]
                    # single row of phage currently being compared
                    df_phage = pd.DataFrame(df_regions_grouped[key].iloc[i])\
                        .transpose().add_prefix('pp_')
                    df_phage = df_phage.rename(columns={'pp_'+col_groupby[0]:col_groupby[0],
                                                        'pp_'+col_groupby[1]:col_groupby[1]})
                    # merge phage data with genes from the same region
                    df_genes_phage = df_genes_single_region.merge(df_phage, on = col_groupby, how='outer')
                    # append genes and phage df to df for output
                    df_out_genes = df_out_genes.append(df_genes_phage)    
            else:
                print('Group not shared between the two dataset: {}'.format([x for x in key]))
#     # no groupby, iterate through each row
#     else:
#         # series of positions for df_search_genes
#         gene_start = df_search_genes[search_cols[0]]
#         gene_stop = df_search_genes[search_cols[1]]
#         for index, row_query in df_query_regions.iterrows():
#             region_start = row_query[query_cols[0]] - flanking_len
#             region_stop = row_query[query_cols[1]] + flanking_len
            
#             # check that start and stop are sequential
#             if (region_start >= region_stop).any():
#                 return print('Error: at index {}, query start position is greater than stop position!') \
#                              .format(index)
#             elif (gene_start >= gene_stop).any():
#                 return print('Error: at index {}, search start position is greater than stop position!') \
#                              .format(index)
                            
#             # compares df_search_genes against query positions individually:
#             for i in range(len(region_start)):
#                 # if start or stop of search is within query region
#                 compare_overlap = (gene_start >= region_start.iloc[i]) & (gene_start <= region_stop.iloc[i]) \
#                                | (gene_stop >= region_start.iloc[i]) & (gene_stop <= region_stop.iloc[i])
#                 # if query is within the entire search region
#                 compare_within = (gene_start <= region_start.iloc[i]) & (gene_stop >= region_stop.iloc[i])
#                 # append matched row to df for output
#                 df_out_genes = df_out_genes.append(df_genes_grouped[key][compare_overlap|compare_within])    
    df_out_genes = df_out_genes.drop_duplicates().reset_index(drop=True)
    return df_out_genes