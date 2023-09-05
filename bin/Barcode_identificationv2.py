import sys
import edlib
from Bio import SeqIO
import pandas as pd

# Define barcode patterns
barcode_patterns = {
    "Barcode05": "AAGGTTACACAAACCCTGGACAAG",
    "Barcode02": "ACAGACGACTACAAACGGAATCGA",
    "Barcode01": "CACAAAGACACCGACAACTTTCTT"
}

def identify_start_adapter(input_file, pattern, kvalue):
    results = []
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            sequence_75 = seq[:120]
            result = edlib.align(pattern, sequence_75, task="path", mode="HW", k=kvalue)



        # Remove overlapping positions
            cleaned_locations = []
            for alignment_result in result["locations"]:
                start_position = alignment_result[0]
                end_position = alignment_result[1]
                overlaps = False

                for existing_result in cleaned_locations:
                    added_start, added_end = existing_result
                    if start_position <= added_end and end_position >= added_start:
                        overlaps = True
                        break

                if not overlaps:
                    cleaned_locations.append((start_position, end_position))

        # Process non-overlapping positions
            for start_position, end_position in cleaned_locations:
                adapter = sequence_75[start_position+1:end_position+1]
                results.append((record.name, adapter, start_position+1, end_position+1))

    return results




if __name__ == "__main__":
    input_file = sys.argv[1]
    barcode_and_kvalues = sys.argv[2].split(',')

    all_results = []
    
    for barcode_and_k in barcode_and_kvalues:
        barcode, kvalue = barcode_and_k.split(':')
        kvalue = int(kvalue)
        
        if barcode in barcode_patterns:
            pattern = barcode_patterns[barcode]
            start_results = identify_start_adapter(input_file, pattern, kvalue)
            df = pd.DataFrame(start_results, columns=['Name', 'Adapter', 'Start', 'End'])
            df.to_csv(f"{barcode}.csv", index=False)
        else:
            print(f"Unknown barcode: {barcode}")






        