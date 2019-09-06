# mem-aware-lbm

## Prerequisite
download scons
```bash
tar xf scons-3.1.1.tar.gz
mkdir build
python3 setup.py install --prefix=./build/
```


## Benchmark
```bash
cd ./examples/benchmarks/cavity3d
make
mpirun -n 3 ./cavity3d 100
```
./cavity3d N
where N is the resolution. The benchmark cases published 
on the Palabos Wiki use N=100, N=400, N=1000, or N=4000.

### Sample output
```
Starting benchmark with 101x101x101 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 3
After 874 iterations: 18.6564 Mega site updates per second.
```