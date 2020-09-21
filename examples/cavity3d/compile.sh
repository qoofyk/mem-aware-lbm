# sequential
cd fuse
make clean;make
cd ../fuse_prism
make
cd ../2step
make
cd ../2step_3parts_prism
make
cd ../2step_whole_prism
make
cd ../2step_whole_prism_unroll
make

# omp, choose one of them to present result
cd ../2step_whole_prism_unroll_omp
make
