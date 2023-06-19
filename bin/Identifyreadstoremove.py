import numpy as np
import pandas as pd
import sys

input_file = sys.argv[1]

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

def getIDtoRemove(df):
    rslt_df = df[df['Type'] == 'Primary'] 
    modeEND=rslt_df.groupby('Ref',as_index=False)['EndPos'].apply(lambda x: x.median())
    rslt_df=rslt_df.merge(modeEND,on='Ref')
    rslt_df.columns =  ["Read","Ref","EndPos","Acc","ReadLen","RefCov","ReadAln","ReadCov","Strand","LeftClip", "RightClip", "IsSec","IsSup","Clipped","Type","med2"]
    strand=rslt_df[(rslt_df.Strand=="-1") | (rslt_df.Acc<94) | (rslt_df.RightClip>49) | (rslt_df.EndPos<rslt_df.med2-100) | (rslt_df.EndPos>rslt_df.med2+800) ]
    return strand


df = parse_seqkit(input_file)
#get ID list of read IDs to remove from bam and dataframe
df6=getIDtoRemove(df)

df7 = df6.iloc[:, 0]
df7.to_csv('ID.txt',index=False,header=False)

