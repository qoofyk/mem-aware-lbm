# sequential
cd fuse
make clean;make
cd ../fuse_prism
make
# cd ../2step
# make clean;make
# cd ../2step_3parts_prism
# make clean;make
# cd ../2step_whole_prism
# make clean;make
# cd ../2step_whole_prism_unroll
# make clean;make

# omp, choose one of them to present result
cd ../2step_whole_prism_omp
make clean;make
