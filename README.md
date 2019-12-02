# mem-aware-lbm

## Prerequisite
download scons
```bash
tar xf scons-3.1.1.tar.gz
cd scons-3.1.1
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

#### 2D cavity
On Bridges interactive job 1 node
```
$ mpirun -n 1 ./cavity2d 127
Starting benchmark with 128x128 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 1
numCells: 16384 After 18311 iterations: 25.3742 Mega site updates per second.

Starting benchmark with 1024x1024 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 1
numCells: 1048576 After 286 iterations: 29.4854 Mega site updates per second.

Starting benchmark with 8192x8192 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 1
numCells: 67108864 After 4 iterations: 30.9548 Mega site updates per second.

```

#### 3D cavity
On my desktop
```
Starting benchmark with 101x101x101 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 3
After 874 iterations: 18.6564 Mega site updates per second.
```

On Bridges interactive job 1 node
```
Starting benchmark with 101x101x101 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 1
After 291 iterations: 9.28136 Mega site updates per second.

Starting benchmark with 101x101x101 grid points (approx. 2 minutes on modern processors).
Number of MPI threads: 3
After 874 iterations: 23.436 Mega site updates per second.
```

### 3D MPI openMP
domain-(1, 1) (1, 11) (1, 11), Test
```
mpirun -n 1 ./cavity3d_step2 10 0 0 2 4
```

### Post-operationi on script
`less 1node.6473699.out | grep 'Mega' | cut -d ":" -f 2 | cut -d " " -f 2`
