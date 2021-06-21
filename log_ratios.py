import pandas as pd
import numpy as np

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC1.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC2.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_1-2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC10.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC5.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_10-5.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC11.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC13.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_11-13.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC15.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC16.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_15-16.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC21.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC22.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_21-22.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC23.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC24.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_23-24.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC1.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC3.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_1-3.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC4.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC5.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_4-5.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC6.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC7.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_6-7_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC10.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC17.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_10-17.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC15.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC16.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_15-16_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC4.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC5.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_4-5.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC3.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC4.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_3-4.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC6.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC7.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_6-7.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC8.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC9.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_8-9.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC12.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC14.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_12-14.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC17.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC18.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_17-18.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC19.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF1BC20.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_19-20.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC8.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC9.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_8-9_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC11.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC13.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_11-13_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC12.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC14.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_12-14_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC18.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC2.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_18-2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC19.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC20.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_19-20_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC21.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC22.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_21-22_2.txt', sep='\t')

##############################################################################################################################

#import control data
df1 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC23.fq_pooled_reads', sep='\t')

#import treatment data
df2 = pd.read_csv('/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/RBAF2_BC24.fq_pooled_reads', sep='\t')
df = pd.merge(df1, df2, on="location", how='outer').fillna(value=1)

#filter on conditions: 
#greater than n in control and less than n in treatment, or less than n in control and greater than n in treatment, or greater than n in control and treatment
df = df.loc[(df['n_x']>=100) & (df['n_y']<=100) | (df['n_x']<=100) & (df['n_y']>=100) | (df['n_x']>=100) & (df['n_y']>=100)]
df = df[df.annotation_x != 'NC']
df = df[df.annotation_y != 'NC']

# #normalize counts by dividing by total pooled reads & take log of normalized ratio (treatment/control)
df['n_x_normalized'] = df['n_x']/(df1['n'].sum())
df['n_y_normalized'] = df['n_y']/(df2['n'].sum())
df['log_normalized_ratio'] = np.log2((df['n_y_normalized'])/(df['n_x_normalized']))

df['annotation_x'] = df.apply(lambda x: x['annotation_x'] if x['annotation_x'] != 1 else x['annotation_y'], axis=1)

df = df.rename(columns={"annotation_x": "gene"})
df = df[['gene','log_normalized_ratio','n_x','n_y','n_x_normalized','n_y_normalized']]

df = df.set_index('gene')

df.to_csv(r'/Users/anna/Documents/buck/rhseq/mapping_90_percent_identity/cutoff_100/log_ratios_23-24_2.txt', sep='\t')

##############################################################################################################################
