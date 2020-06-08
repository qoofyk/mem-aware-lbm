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


mflups=[]
mb=[]
setup=[]
dim=[]
threads=[]

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
    if "MPI" in s:
      if (len(temp) > 9): 
        threads.append(temp[9])
        
np_mflups = np.asarray(mflups, dtype=np.float32)
np_mb = np.asarray(mb, dtype=np.int64)
np_setup = np.asarray(setup, dtype=np.float32)
np_threads = np.asarray(threads, dtype=np.int64)
np_dim = np.asarray(dim)

repeat=3
num_exp = np_mflups.shape[0] // repeat
print(np_mflups.shape, np_mflups.shape[0], num_exp, repeat)
# print(np_mflups[0 : repeat * num_exp])
print(jobid)
# np.savetxt(sys.stdout.buffer, np.flipud(np_mflups), fmt="%.4f")

#cube
# np.savetxt(sys.stdout.buffer, np.flipud(np_mflups[0 : repeat * num_exp].reshape(repeat, num_exp).transpose()), fmt="%.4f")
# if np_threads.size != 0:
#   np.savetxt(sys.stdout.buffer, np.flipud(np_threads[0 : repeat * num_exp].reshape(repeat, num_exp).transpose())[:,0], fmt="%d")
# if np_dim.size != 0:
#   np.savetxt(sys.stdout.buffer, np.flipud(np_dim[0 : repeat * num_exp].reshape(repeat, num_exp).transpose())[:,0], fmt="%s")

#rectangular cs
group=7
# for i in range(0, repeat * num_exp, repeat * group):
#   end = np.where(i + repeat * group < repeat * num_exp, i + repeat * group, repeat * num_exp)
#   num_exp_per_group = (end - i) // repeat
#   if np_dim.size != 0:
#     np.savetxt(sys.stdout.buffer, np.flipud(np_dim[i : end].reshape(num_exp_per_group, repeat))[:,0], fmt="%s")
#   print()

# for i in range(0, repeat * num_exp, repeat * group):
#   end = np.where(i + repeat * group < repeat * num_exp, i + repeat * group, repeat * num_exp)
#   num_exp_per_group = (end - i) // repeat
#   if np_threads.size != 0:
#     np.savetxt(sys.stdout.buffer, np.flipud(np_threads[i : end].reshape(num_exp_per_group, repeat))[:,0], fmt="%d")
#   print()

for i in range(0, repeat * num_exp, repeat * group):
  end = np.where(i + repeat * group < repeat * num_exp, i + repeat * group, repeat * num_exp)
  num_exp_per_group = (end - i) // repeat
  np.savetxt(sys.stdout.buffer, np.flipud(np_mflups[i : end].reshape(num_exp_per_group, repeat)), fmt="%.4f")
  print()

# for i in range(0, repeat * num_exp, repeat * group):
#   end = np.where(i + repeat * group < repeat * num_exp, i + repeat * group, repeat * num_exp)
#   num_exp_per_group = (end - i) // repeat
#   np.savetxt(sys.stdout.buffer, np.flipud(np_mflups[i : end].reshape(num_exp_per_group, repeat)), fmt="%.4f")
#   if np_threads.size != 0:
#     np.savetxt(sys.stdout.buffer, np.flipud(np_threads[i : end].reshape(num_exp_per_group, repeat))[:,0], fmt="%d")
#   if np_dim.size != 0:
#     np.savetxt(sys.stdout.buffer, np.flipud(np_dim[i : end].reshape(num_exp_per_group, repeat))[:,0], fmt="%s")

# for i in range(0, 21 * 50, 21):
#   np.savetxt(sys.stdout.buffer, np_mflups[i : (i + 21)].reshape(3,7).transpose(), fmt="%.4f")
#   print()
#   print()

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