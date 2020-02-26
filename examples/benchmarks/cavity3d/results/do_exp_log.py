#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:57:35 2019
@author: fuyuan
"""
import re
import os, sys
import numpy as np

#some_list = ['abc-123', 'def-456', 'ghi-789', 'abc-456']
#matching = [s for s in some_list if "abc" in s]
#print(matching)

f = open(sys.argv[1], 'r')
jobid=sys.argv[1].split(".")[-2]
print(jobid)

mflups=[]

for line in f.readlines():
    temp = line.split(" ")
    for s in temp:
        if "Mega" in s:
            mflups.append(temp[3])
        
npmflups = np.asarray(mflups, dtype=np.float32)

for i in range(0, 21 * 22, 21):
  np.savetxt(sys.stdout.buffer, npmflups[i : (i + 21)].reshape(3,7).transpose(), fmt="%.4f")
  print()
  print()