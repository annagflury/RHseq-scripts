import pandas as pd
import numpy as np
import csv
from functools import reduce

###  ED3077 data  ######
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_1-3.txt', sep='\t')
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_4-5.txt', sep='\t')
df3 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_6-7_2.txt', sep='\t')
df4 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_10-17.txt', sep='\t')
df5 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_15-16.txt', sep='\t')

df6 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_1-2.txt', sep='\t')
df7 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_10-5.txt', sep='\t')
df8 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_11-13.txt', sep='\t')
df9 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_15-16_2.txt', sep='\t')
df10 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_21-22.txt', sep='\t')
df11 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_23-24.txt', sep='\t')
df = df1.append([df2, df3, df4, df5, df6, df7, df8, df9, df10, df11])
 
####  N2 data  ######
# df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_3-4.txt', sep='\t')
# df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_6-7.txt', sep='\t')
# df3 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_8-9_2.txt', sep='\t')
# df4 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_12-14.txt', sep='\t')
# df5 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_17-18.txt', sep='\t')
# df6 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_19-20.txt', sep='\t')
# df7 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_8-9.txt', sep='\t')
# df8 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_11-13_2.txt', sep='\t')
# df9 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_12-14_2.txt', sep='\t')
# df10 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_18-2.txt', sep='\t')
# df11 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_19-20_2.txt', sep='\t')
# df12 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_21-22_2.txt', sep='\t')
# df13 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_23-24_2.txt', sep='\t')

# df = df1.append([df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13])

df = df.groupby(['gene'])['log_normalized_ratio'].apply(lambda floats: ",".join(list(map(str, floats)))).reset_index()
df['gene'] = df['gene'].map(lambda x: x.lstrip('n2'))

df=df.set_index('gene')
print(df)


pd.DataFrame.to_csv(df, '/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/ed_combined_reps.txt', sep=',', na_rep='.', index=True)




