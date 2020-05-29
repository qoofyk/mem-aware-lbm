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
mb=[]
setup=[]
dim=[]

for line in f.readlines():
  temp = line.split(" ")
  for s in temp:
    if "Mega" in s:
        mflups.append(temp[3]) #start from index0
    if ("MB" in s):
      mb.append(temp[10]) # extra space in this line is splited
      dim.append(temp[3])
    if "cavitySetup" in s:
    	setup.append(temp[4].strip('\n'))
        
np_mflups = np.asarray(mflups, dtype=np.float32)
np_mb = np.asarray(mb, dtype=np.int64)
np_setup = np.asarray(setup, dtype=np.float32)

for i in range(0, 21 * 50, 21):
  np.savetxt(sys.stdout.buffer, np_mflups[i : (i + 21)].reshape(3,7).transpose(), fmt="%.4f")
  print()
  print()

# for i in range(0, 21 * 31, 21):
#   np.savetxt(sys.stdout.buffer, np_mb[i : (i + 21)].reshape(3,7).transpose()[:, 0], fmt="%d")
#   print()
#   print()

# for i in range(0, 21 * 50, 21):
#   np.savetxt(sys.stdout.buffer, np_setup[i : (i + 21)].reshape(3,7).transpose()[:, 0], fmt="%.6f")
#   print()
#   print()

# for i in range(0, 21 * 50, 21):
#   for j in range(i, i+7, 1):
#     print(dim[j])
#   print()
#   print()