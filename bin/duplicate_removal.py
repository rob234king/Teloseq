import sys
import pysam

# Check for the correct number of command-line arguments
if len(sys.argv) != 3:
    print("Usage: python deduplicate_bam.py input.bam output.bam")
    sys.exit(1)

# Input and output BAM file paths
input_bam = sys.argv[1]
output_bam = sys.argv[2]

# Open the input BAM file
with pysam.AlignmentFile(input_bam, "rb") as bamfile:

    # Create a dictionary to keep track of seen primary alignments by read name
    primary_seen = {}

    # Create a set to store seen alignment lines
    seen_lines = set()

    # Create an output BAM file for writing
    with pysam.AlignmentFile(output_bam, "wb", header=bamfile.header) as output:

        # Iterate over the alignments in the input BAM file
        for alignment in bamfile:

            # Get the read name
            read_name = alignment.query_name

            # Check if the alignment is a primary alignment
            is_primary = not alignment.is_secondary and not alignment.is_supplementary

            if is_primary:
                # If it's a primary alignment, check if it's the first primary for this read
                if read_name not in primary_seen:
                    primary_seen[read_name] = True  # Mark it as seen
                    output.write(alignment)  # Write the primary alignment to the output BAM file

            # Check if the alignment line is identical to any previous alignment line
            alignment_line = alignment.to_string()
            if alignment_line not in seen_lines:
                seen_lines.add(alignment_line)  # Mark it as seen

# Index the deduplicated BAM file
pysam.index(output_bam)
