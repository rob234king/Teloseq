import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

input_file = sys.argv[1]
input_file2 = sys.argv[2]
input_file3 = sys.argv[3]

def parse_seqkit(fname):
    cols = {
        'Read':str, 'Ref':str, 'Acc':float, 'ReadLen':int, 'RefCov':float, 'EndPos':int,'RightClip':int,'LeftClip':int, 'Strand':str,
        'ReadAln':int, 'ReadCov':float, 'IsSec':bool, 'IsSup':bool}
    df = pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys())
    #add column
    df['Clipped'] = df['ReadLen'] - df['ReadAln']
    #new column primary, secondary or supplementary
    df['Type'] = 'Primary'
    df.loc[df['IsSec'], 'Type'] = 'Secondary'
    df.loc[df['IsSup'], 'Type'] = 'Supplementary'
    return df

def parse_motif(fname):
    cols = { 'seqID':str, 'start':int}
    df = pd.read_csv(fname, sep="\t", dtype=cols, usecols=cols.keys())
    df.columns =  ["Read","Telomere_length"]
    return df

#parse both files
dfreads = parse_seqkit(input_file)
dfmotif = parse_motif(input_file2)

#keep only primary alignments
rslt_df = dfreads[dfreads['Type'] == 'Primary'] 
rslt_df2 = rslt_df[rslt_df['Strand'] == '1'] 
dfmerged=rslt_df2.merge(dfmotif,on='Read')
#print(rslt_df.head())
#coverage
Summary=dfmerged.groupby('Ref').agg(**{
        'Coverage': ('Read', 'count'),
        'avg_accuracy': ('Acc', 'mean'),
        'SD_accuracy': ('Acc', 'std'),
        'Telomere_mean': ('Telomere_length', 'mean'),
        'Telomere_sd': ('Telomere_length', 'std'),
        'Telomere_max': ('Telomere_length', 'max')})

#import reference fai as dataframe for list of chr that should be in list
cols = { 'Ref':str, 'length':int, 'length2':int, 'length3':int, 'length4':int}
referenceset = pd.read_csv(input_file3, sep="\t", header=None, usecols=cols.keys(), names=cols.keys(), dtype=cols)

Summary2=referenceset['Ref'].to_frame().merge(Summary,on='Ref',how='left').fillna(0)

Summary2.to_csv('Coverage.csv',index=False)


#summarise read coverage for each chr arm with boxplot
##FIGURES
fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(20)
 # Divide the figure into a 2x1 grid, and give me the first section
ax1 = fig.add_subplot(111)
dfmerged.boxplot(column='Telomere_length',by='Ref',fontsize=12,ax=ax1,rot=90)  
ax1.xaxis.set_label_text("Chr arm")
fig.subplots_adjust(hspace=0.5)

fig.savefig('Boxplot_of_Telomere_length.pdf')