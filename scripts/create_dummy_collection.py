import os
import sys


directory = os.fsencode(sys.argv[1])
with open(sys.argv[2], 'w') as table:
    for file in os.scandir(directory):
        filename = os.fsdecode(file)

        if filename.endswith(".fa"):
            print(os.path.basename(filename)[:-3])
            with open(file, 'r') as contigs_file:
                for r in contigs_file.readlines():
                    if r.startswith(">"):
                        r = r[1:].strip()
                        table.write(r + "\t" + "ALL" + "\n")

        if filename.endswith(".fasta"):
            print(os.path.basename(filename)[:-5])
            with open(file, 'r') as contigs_file:
                for r in contigs_file.readlines():
                    if r.startswith(">"):
                        r = r[1:].strip()
                        table.write(r + "\t" + "ALL" + "\n")
