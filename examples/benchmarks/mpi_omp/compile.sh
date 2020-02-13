rm -f fs_*
mpiicpc hello.c -qopenmp -o hello
icpc false_sharing_omp_nopad.cpp -g -std=c++11 -qopenmp -o fs_nopad_vec
icpc false_sharing_omp_pad.cpp -g -std=c++11 -qopenmp -o fs_pad_vec

icpc false_sharing_omp_nopad_array.cpp -g -std=c++11 -qopenmp -o fs_nopad_arr
icpc false_sharing_omp_pad_array.cpp -g -std=c++11 -qopenmp -o fs_pad_arr
