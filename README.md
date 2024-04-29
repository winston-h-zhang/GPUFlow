### Build guide
A C++17 compliant compiler is needed to build. 
```bash
mkdir build && cd build
cmake .. && make
cd ..
```

We also require some version of OpenMP.

### Example usage
The `maxflow` binary reads a max flow problem in the [DIMACS][Dimacs_descr] format and solves it using the specified algorithm. For example, the following uses the sequential push-relabel solver:
```bash
./maxflow prf -f example.inp
``` 
For the list of possible algorithms and other options, use:
```bash
> ./maxflow --help
usage: ./maxflow <solver> [-f <path>] [-p <number>]

list of possible solvers: [prf, ppr]

prf:    push-relabel algorithm with FIFO vertex selection
ppr:    parallel push-relabel algorithm

[-f <path>]:    specify path to a maxflow problem instance in DIMACS format. Reads from stdin if omitted.
[-p <number>]:  specify max number of threads for parallel solvers. Ignored for sequential solvers. Default = number of hw threads.
```

### Dataset

Look at the report for more details. You'll have to actually download and unzip the data by hand to run the tests. 

- DIMACS generated dataset
    - http://archive.dimacs.rutgers.edu/pub/netflow/
    - Washington-RLG-Wide: easy family of random grid graphs
    - RMF-Wide: hardest DIMACS dataset for modern push–relabel algorithms [Gol09]

- Image segmentation dataset: https://vision.cs.uwaterloo.ca/data/maxflow
    - liver – 128x128x119, short paths
    - adhead – 256x256x192, long paths
    - LB07-bunny – 3D volumetric surface of the Stanford Bunny

We've also included a benchmarking script `run.sh` for ease of testing. However, you will have to manually copy to benchmarks into the following file structure:

```bash
.
├── generate
│   ├── adhead.n6c10
│   │   ├── adhead.n6c10.max
│   │   └── adhead.n6c10.sol
│   ├── gengraph
│   ├── LB07-bunny-med
│   │   ├── LB07-bunny-med.max
│   │   └── LB07-bunny-med.sol
│   ├── liver.n6c10
│   │   ├── liver.n6c10.max
│   │   └── liver.n6c10.sol
│   └── washington.n6c7
│       └── washington.n6c7.max
...
```

To run the benchmarks, setup the dataset as shown above, then run `./run.sh`.

[Dimacs_descr]: http://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm