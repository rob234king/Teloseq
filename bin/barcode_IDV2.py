import sys
import pandas as pd
from itertools import combinations

def main():
    if len(sys.argv) < 1:
        print("Usage: python script.py input_csv1.csv input_csv2.csv ...")
        return

    input_csv_files = sys.argv[1:]
    all_unique_names = {}
    shared_counts = {}

    for input_csv in input_csv_files:
        barcode_name = input_csv.replace(".csv", "")
        df = pd.read_csv(input_csv, usecols=['Name'], delimiter=',')
        unique_names = set(df['Name'].unique())
        all_unique_names[barcode_name] = unique_names

    for name1, name2 in combinations(all_unique_names.keys(), 2):
        common_names = all_unique_names[name1].intersection(all_unique_names[name2])
        shared_counts[f"Shared {name1} and {name2}"] = len(common_names)

    for barcode_name, name_list in all_unique_names.items():
        other_name_lists = [names for name, names in all_unique_names.items() if name != barcode_name]
        common_in_others = set.union(*other_name_lists)
        unique_names_exclusive = name_list - common_in_others
        shared_counts[f"Unique {barcode_name} (exclusive)"] = len(unique_names_exclusive)
        unique_names_output_file = f'unique_names_{barcode_name}_exclusive.txt'
        with open(unique_names_output_file, 'w') as f:
            for name in unique_names_exclusive:
                f.write(name + '\n')

    # Write shared counts to a file
    with open('Barcode_statistics.txt', 'w') as f:
        for label, count in shared_counts.items():
            f.write(f"{label}: {count}\n")

    # Print statistics to the command line
    for label, count in shared_counts.items():
        print(f"{label}: {count}")

if __name__ == "__main__":
    main()

