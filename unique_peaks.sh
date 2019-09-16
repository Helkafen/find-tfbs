#!/bin/bash

echo "Removing the regions that are present in any other bed file. Creating a 4 columns output (chr,start,end,empty)"

for f in $*; do
    cat `ls *.bed | grep -v ^${f}$` | sort -k1,1 -k2,2n | cut -f 1,2,3 | bedtools merge > $f.others # Merge all the other BED files
    bedtools subtract -a $f -b $f.others | awk -v OFS='\t' '{print $1, $2, $3, ""}' > $f.unique     # Create a new BED file that only contains the regions that are unique to f
    echo Created $f.unique
    rm $f.others
done

echo "Done"