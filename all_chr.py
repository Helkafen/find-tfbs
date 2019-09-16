#!/usr/bin/env python
import sys

command = " ".join(sys.argv[1:])

for c in ["chr%s" % str(i) for i in range(1,23)] + ["chrX"]:
    print(command.replace("chr1", c))