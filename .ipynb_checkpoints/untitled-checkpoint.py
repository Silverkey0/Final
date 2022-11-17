import pandas as pd
 
tsv1 = pd.read_csv("/Users/henry/Documents/GitHub/Final/Gene_Counts_Modified/Gene_Counts1.tsv", sep='\t')
tsv2 = pd.read_csv("/Users/henry/Documents/GitHub/Final/Gene_Counts_Modified/Gene_Counts2.tsv", sep='\t')    
Output_df = pd.merge(tsv1, tsv2, on='gene_id',
                     how='outer')
Output_df.to_csv("/Users/henry/Documents/GitHub/Final/pancrease_sample/output.tsv",
                 sep="\t", header=True,
                 index=False)