#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:40:25 2021

@author: balthasar
"""

import sys

with open(sys.argv[1]) as f:
    filterlines = {l[:-1] for l in f}

newfastq = []
save = False
with sys.stdin as f:
    with open(sys.argv[2], "w") as fw:
        for line in f:
            ID = line[1:].split(" ")[0]
            if ID in filterlines:
                fw.write(line)
                fw.write(f.readline())
                fw.write(f.readline())
                fw.write(f.readline())
