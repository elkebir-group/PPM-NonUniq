#!/usr/bin/python
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: %s <INPUT_TREE>\n" % sys.argv[0])
    sys.exit(1)

with open(sys.argv[1]) as f:
    f.readline()
    edges = 0
    lines = []
    count = True
    for line in f:
        if line.find("GL") == -1 and line.find("_") == -1:
            if line.find("#leaves") != -1:
                count = False
            elif count:
                edges += 1
                lines.append(line)
    #print edges, "#edges"
    for line in lines:
        print line,
