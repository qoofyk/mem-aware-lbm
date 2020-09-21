# sequential
rm ../../lib/*.a
cd fuse
make clean
cd ../fuse_prism
make clean
cd ../2step
make clean
cd ../2step_3parts_prism
make clean
cd ../2step_whole_prism
make clean
cd ../2step_whole_prism_unroll
make clean

# omp, choose one of them to present result
cd ../2step_whole_prism
make clean