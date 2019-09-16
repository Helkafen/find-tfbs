import gzip
import sys
import re

# This script selects the peaks with large variation in number of TFBS

input_filename = sys.argv[1]

f = gzip.open(input_filename, 'rt')
for line in f:
    if "CHROM" in line:
        print(line.strip())
    else:
        counts = [int(x) for x in re.findall("COUNT=[\d,]+;",line)[0][:-1].replace("COUNT=","").split(",")]
        a = min(counts)
        b = max(counts)
        if b - a >= 3:
            print(line.strip())