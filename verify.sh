#!/usr/bin/env bash

# First created: 
# Last modified: 2018 Aug 30

# Author: Yuankun Fu
# email: qoofyk@gmail.com

icc -DORG pillar_test.cpp -o org
icc -DPILLAR pillar_test.cpp -o pillar

./org 14 6 14 8 > tmp1
./pillar 14 6 14 8 > tmp2
diff tmp1 tmp2

./org 62 30 62 32 > tmp1
./pillar 62 30 62 32 > tmp2
diff tmp1 tmp2

./org 126 126 126 32 > tmp1
./pillar 126 126 126 32 > tmp2
diff tmp1 tmp2

./org 28 126 126 32 > tmp1
./pillar 28 126 126 32 > tmp2
diff tmp1 tmp2