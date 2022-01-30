from Bio import SeqIO
import sys

def usage():
    print("Usage: python filter_assembly_using_bin_list.py <file_with_list_of_bins> <assembly_file> <output_assembly_file>")

if len(sys.argv) != 4:
    usage()
    exit(0)

used = set()

with open(sys.argv[1], 'r') as file_with_paths:
    for filename in file_with_paths.readlines():
        filename = filename.strip()
        for record in SeqIO.parse(filename, "fasta"):
            used.add(record.id)



with open(sys.argv[2]) as assembly_file:
    with open(sys.argv[3], 'w') as assembly_file_new:
        for record in SeqIO.parse(assembly_file, "fasta"):
            if record.id is used:
                SeqIO.write(assembly_file_new, record)
