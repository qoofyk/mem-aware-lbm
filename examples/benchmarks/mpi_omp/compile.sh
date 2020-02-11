mpiicpc hello.c -qopenmp -o hello
icpc false_sharing_omp_nopad.c -qopenmp -o fs_nopad
icpc false_sharing_omp_pad.c -qopenmp -o fs_pad
