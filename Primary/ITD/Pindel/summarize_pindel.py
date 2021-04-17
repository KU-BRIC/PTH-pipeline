#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 16:37:51 2021

@author: balthasar
"""

csv_content = ["SAMPLE\tCHR\tPOS\tSVLEN\tDEPTH\n"]

import os
for filename in os.listdir('pindel_TD'):
    with open('pindel_TD/'+filename) as f:
        for line in f:
            if not line[0] == "#":
                sp = line[:-1].split("\t")
                attr = sp[7].split(";")
                tmpdict = {a.split("=")[0]:a.split("=")[1] for a in attr}
                
                # filename is not good
                csv_content.append("\t".join([filename.split("_")[0],sp[0],f"{sp[1]}-{tmpdict['END']}",tmpdict["SVLEN"],sp[9].split(":")[1]+"\n"]))

with open("pindel_results.csv","w") as f:
    f.writelines(csv_content)
