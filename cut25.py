#run  by running 'python cut25.py input.fastq.gz output.fastq.gz'
#meant to subsample barcodes for rarefaction by reads
#libraries
import gzip
import random
import sys


def select_random_reads(input_file, output_file, percent=25):
    # Open the input and output files
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        # Store all the reads in a list
        reads = []
        while True:
            # Read 4 lines (one read) at a time (fastq format)
            read = [infile.readline().strip() for _ in range(4)]
            if not read[0]:
                break
            reads.append(read)

        # Determine the number of reads to select (10%), does not account for read length
        num_reads_to_select = len(reads) * percent // 100

        # Randomly shuffle and select the reads
        selected_reads = random.sample(reads, num_reads_to_select)

        # Write selected reads to the output file
        for read in selected_reads:
            for line in read:
                outfile.write(f"{line}\n")

    print(f"Selected {num_reads_to_select} out of {len(reads)} reads and saved to {output_file}")
