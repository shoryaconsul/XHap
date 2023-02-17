"""
This script reads in the fragment file generated by extractHAIRS from
HapCUT2 and generates the corresponding file for SDHap, by determining
the read count and column count.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="File generated by extractHAIRS")
parser.add_argument("-o", "--output", help="File for SDHap")

args = parser.parse_args()

# frag_file = "SDHap/Fosmid_all_chromosomes/chr15.matrix.SORTED"
# out_file = "SDHap/Fosmid_all_chromosomes/chr15.matrix_sdhap"

frag_file = args.input
out_file = args.output

read_start = -1
read_count = 0
col_num = 0
col = []
line_list = []

with open(frag_file, "r") as f:
    while True:
        line = f.readline().split()
        if not line:
            break
        else:
            line_list.append(line)

        if len(line) >= 5:
            if read_start == -1:
                read_start = len(line_list) - 1  # line with first read
            read_count = read_count + 1
            for i in range(2, len(line) - 1, 2):
                col.append(int(line[i]))
                if int(line[i]) > col_num:
                    col_num = int(line[i])

with open(out_file, "w") as f:
    f.write(str(read_count + 1))
    f.write('\n')
    f.write(str(col_num))
    f.write('\n')
    for line in line_list[read_start:]:
        f.write(' '.join(line))
        f.write('\n')
