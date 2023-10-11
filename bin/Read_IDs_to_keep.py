import pandas as pd
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description='Filter and extract data from input files.')
parser.add_argument('text_file', help='Input text file name')
parser.add_argument('read_ids_to_remove', help='Input file containing read IDs to remove')
parser.add_argument('bed_file_cutsites', help='Input bed file name')
parser.add_argument('bed_file_telomereend', help='Input bed file name')
parser.add_argument('output_file', help='Output file name')

# Parse the command line arguments
args = parser.parse_args()

# Read the text file into a DataFrame
text_file = pd.read_csv(args.text_file, sep='\t')

# Read the bed file into a DataFrame
bed_filecut = pd.read_csv(args.bed_file_cutsites, sep='\t', header=None, names=['Ref', 'StartPos','junk'],dtype={'Ref': str})
#replace flicked reads
# replacements = [
#     {'Ref': 'chr21_PATERNAL_P', 'StartPos': 14500, 'junk': 14500},
#     {'Ref': 'chr21_MATERNAL_P', 'StartPos': 9000, 'junk': 9000}
# ]
# # Iterate through the list of replacements and update the DataFrame
# for replacement in replacements:
#     condition = (bed_filecut['Ref'] == replacement['Ref'])
#     bed_filecut.loc[condition, ['Ref', 'StartPos', 'junk']] = [replacement['Ref'], replacement['StartPos'], replacement['junk']]


bed_filetel = pd.read_csv(args.bed_file_telomereend, sep='\t', header=None, names=['Ref', 'telPos1','junk2'],dtype={'Ref': str})

# Merge the two DataFrames based on the 'Ref' column
merged_data = pd.merge(text_file, bed_filecut, on='Ref', how='inner')
merged_data2 = pd.merge(merged_data, bed_filetel, on='Ref', how='inner')

# Read the file containing read IDs to remove into a DataFrame
read_ids_to_remove = pd.read_csv(args.read_ids_to_remove, sep='\t', header=None, names=['Read'], dtype={'Read': str})


# Save the result to the specified output file
#merged_data2.to_csv("raw.txt", sep='\t', index=False, header=True)
# Filter the rows where 'EndPos' is within 100 of 'StartPos' AND 'IsSec' and 'IsSup' both have a value of 0
filtered_data = merged_data2[((merged_data2['StartPos'] - 15 <= merged_data2['EndPos']) &
                             (merged_data2['EndPos'] <= merged_data2['StartPos'] + 15)) &
                            ((merged_data2['EndPos'] - merged_data2['ReadAln']) < merged_data2['telPos1']) &
                            (merged_data2['IsSec'] == 0) &
                            (merged_data2['IsSup'] == 0)]

# Apply an independent filter for 'ReadAln' of 7000 for the chromosome arms 'chr21_PATERNAL_P' and 'chr21_MATERNAL_P'
chr21_filter = (merged_data2['Ref'].isin(['chr21_PATERNAL_P', 'chr21_MATERNAL_P'])) & (merged_data2['ReadAln'] > 7000) & (merged_data2['MapQual'] > 58) 


# Filter out "Read" matches based on the read_ids_to_remove DataFrame
filtered_data = filtered_data[~filtered_data['Read'].isin(read_ids_to_remove['Read'])]

# Add the filtered chromosome arms to the filtered_data set
filtered_data = pd.concat([filtered_data, merged_data2[chr21_filter]])


# Extract the 'Read' column
output_data = filtered_data[['Read']]

# Save the result to the specified output file
output_data.to_csv(args.output_file, sep='\t', index=False, header=False)

