#!/usr/bin/env bash

# First created: 
# Last modified: 2020 Sep 6

# Author: Yuankun Fu
# email: qoofyk@gmail.com

############### sequential
cd fuse
sed '/compileFlags =*/a\step2_whole_Flags = false\nstep2_3parts_Flags = false\nstep2_omp_Flags = false\nstep2_unroll_Flags = false\nstep2_pyramid_Flags = false' ../bridges.seq.Makefile > Makefile

cd ../fuse_prism
sed '/compileFlags =*/a\step2_whole_Flags = false\nstep2_3parts_Flags = false\nstep2_omp_Flags = false\nstep2_unroll_Flags = false\nstep2_pyramid_Flags = false' ../bridges.seq.Makefile > Makefile

cd ../2step
sed '/compileFlags =*/a\step2_whole_Flags = false\nstep2_3parts_Flags = true\nstep2_omp_Flags = false\nstep2_unroll_Flags = false\nstep2_pyramid_Flags = false' ../bridges.seq.Makefile > Makefile

cd ../2step_3parts_prism
sed '/compileFlags =*/a\step2_whole_Flags = false\nstep2_3parts_Flags = true\nstep2_omp_Flags = false\nstep2_unroll_Flags = false\nstep2_pyramid_Flags = true' ../bridges.seq.Makefile > Makefile

cd ../2step_whole_prism
sed '/compileFlags =*/a\step2_whole_Flags = true\nstep2_3parts_Flags = false\nstep2_omp_Flags = false\nstep2_unroll_Flags = false\nstep2_pyramid_Flags = true' ../bridges.seq.Makefile > Makefile

cd ../2step_whole_prism_unroll
sed '/compileFlags =*/a\step2_whole_Flags = true\nstep2_3parts_Flags = false\nstep2_omp_Flags = false\nstep2_unroll_Flags = true\nstep2_pyramid_Flags = true' ../bridges.seq.Makefile > Makefile

# ############ # omp, choose one of them to present result
cd ../2step_whole_prism_unroll_omp
sed '/compileFlags =*/a\step2_whole_Flags = true\nstep2_3parts_Flags = false\nstep2_omp_Flags = true\nstep2_unroll_Flags = true\nstep2_pyramid_Flags = true' ../bridges.omp.Makefile > Makefile
