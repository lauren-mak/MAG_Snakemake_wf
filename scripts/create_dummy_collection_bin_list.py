import os
import sys


file_with_bins = sys.argv[1]
with open(file_with_bins, 'r') as bin_list:
    with open(sys.argv[2], 'w') as table:
        for r in bin_list.readlines():
            filename = r.strip()

            if filename.endswith(".fa"):
                print(os.path.basename(filename)[:-3])
                with open(filename, 'r') as contigs_file:
                    for r in contigs_file.readlines():
                        if r.startswith(">"):
                            r = r[1:].strip()
                            table.write(r + "\t" + "ALL" + "\n")

            if filename.endswith(".fasta"):
                print(os.path.basename(filename)[:-5])
                with open(filename, 'r') as contigs_file:
                    for r in contigs_file.readlines():
                        if r.startswith(">"):
                            r = r[1:].strip()
                            table.write(r + "\t" + "ALL" + "\n")
