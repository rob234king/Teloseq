import pandas as pd

# Read the data from your file into a pandas DataFrame
df = pd.read_csv('RK_motif.txt', delimiter='\t', header=None,
                 names=['seqID', 'patternName', 'pattern', 'strand', 'start', 'end', 'matched'])

# Convert the 'start' column to numeric
df['start'] = pd.to_numeric(df['start'], errors='coerce')

# Group the DataFrame by 'seqID'
grouped = df.groupby('seqID')

# Iterate through groups and count rows within 500 of the start values
for name, group in grouped:
    count = group[(group['start'] - group['start'].iloc[0]).abs() <= 500].shape[0]
    if count >= 5:
        print(name)