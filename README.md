# New small regular graphs of given girth: the cage problem and beyond

This repository contains code and data related to the manuscript "New small regular graphs of given girth: the cage problem and beyond". 
Code can be found in the `src` directory, and data in the `data` directory.

The following notation will be used:
- $(k,g)$-graph: A $k$-regular graph with girth $g$
- $(k,g,\lambda_e)$-graph or $egr$-graph: A $k$-regular graph with girth $g$ where every edge is contained in $\lambda_e$ cycles of length $g$.
- $(k,g,\lambda_v)$-graph or $vgr$-graph: A $k$-regular graph with girth $g$ where every vertex is contained in $\lambda_v$ cycles of length $g$.
- $(k,g,\underline{g+1})$-graph: A $(k,g)$-graph without $(g+1)$-cycles.

The next sections detail how the code can be used and what data is available.

## Data
The `data` directory is structured as follows:
```
data/
├── groups
│   ├── groups_1-50.txt
├── graphs
│   ├── cage
│   ├── egr
│   ├── kgnogpo
│   ├── spectra.zip
└── └── vgr.zip
```

The groups subdirectory contains all non-isomorphic groups of order 1 to 50 with up to 2000 automorphisms for each group.

The graphs subdirectory contains all graphs stored in graph6 format found in the manuscript

They are organized by type and use the following naming convention:
- `cage`: Small $(k,g)$-graphs where `nW_kX_gY.g6` denotes a file containing $(X,Y)$-graphs with $W$ vertices.
- `egr`: Small $(k,g,\lambda_e)$-graphs where `nW_kX_gY_lZ.g6` denotes a file containing $(X,Y,Z)$-graphs with $W$ vertices.
- `kgnogpo`: Small $(k,g,\underline{g+1})$-graphs where `nW_kX_gY.g6` denotes a file containing $(X,Y,\underline{Y+1})$-graphs with $W$ vertices.
- `spectra.zip`: Zipped to preserve space. Contains $(k,g)$-graphs where `nW_kX_gY.g6` denotes a file containing $(X,Y)$-graphs with $W$ vertices.
- `vgr.zip`: Zipped to preserve space. Contains small $(k,g,\lambda_v)$-graphs where `nW_kX_gY_lZ.g6` denotes a file containing $(X,Y,Z)$-graphs with $W$ vertices.

## Code
The `src` directory contains multiple files and subdirectories. 
The `src/GAP` subdirectory contains GAP code to generate the groups including their automorphisms.
The `src/nautyAndMultigraph` subdirectory contains slightly modified versions of the Nauty library and the generator multigraph.

The remaining code revolves around generating lifts to obtain $k$-regular graphs achieving some minimum girth. 
It can be compiled from the top-level directory using the following commands:
```bash
cmake
make
```

The executables are available under the `bin` directory, but should be executed using:
```bash
. ./run.sh <partition> <modulo> <n> <k> <groupMin> <groupMax> <minGirth>
```
This will construct all non-isomorphic lifts achieving the given `minGirth` using $k$-regular multigraphs on $n$ vertices with loops and semi-edges 
and groups with order between `groupMin` and `groupMax` (inclusive). The parameters `partition` and `modulo` can be used to evenly split the base graphs over multiple invocations.

For example, to construct all lifts with girth at least 5 using 3-regular base graphs on 2 vertices 
and groups with up to 5 elements, you can run the following command:
```bash 
. ./run.sh 0 1 2 3 1 5 5
```

This will output the following information:
```
Partition: 0
Modulo: 1
n: 2
k: 3
groupMin: 1
groupMax: 5
minGirth: 5

Constructing multigraphs...
Found 1 base multigraphs for degree distribution: 0 0 2 
Resulting in 1 new multigraphs for this partition.
Found 0 base multigraphs for degree distribution: 0 1 1 
Resulting in 0 new multigraphs for this partition.
Found 1 base multigraphs for degree distribution: 0 2 0 
Resulting in 1 new multigraphs for this partition.
Found 0 base multigraphs for degree distribution: 1 0 1 
Resulting in 0 new multigraphs for this partition.
Found 0 base multigraphs for degree distribution: 1 1 0 
Resulting in 0 new multigraphs for this partition.
Found 1 base multigraphs for degree distribution: 2 0 0 
Resulting in 3 new multigraphs for this partition.
Found 5 3-regular multigraphs with 2 vertices (2 without semi edges) in 0ms

Graph 1 with at most 0 semi edges at any vertex.
Adjacency:
0:  1, 1, 1, 
1:  0, 0, 0, 
Graph 2 with at most 1 semi edges at any vertex.
Adjacency:
0:  1, 1, 0, 
1:  0, 0, 1, 
Graph 3 with at most 2 semi edges at any vertex.
Adjacency:
0:  1, 0, 0, 
1:  0, 1, 1, 
Graph 4 with at most 2 semi edges at any vertex.
Adjacency:
0:  1, 0, 0, 
1:  0, 1, 1, 
Graph 5 with at most 0 semi edges at any vertex.
Adjacency:
0:  1, 0, 0, 
1:  0, 1, 1, 

Reading groups...
Found 6 groups in 0ms.

Calculating usable edge automorphisms...
Graph 1 has 5 generators, 1 automorphisms, 1 orbits.
Graph 2 has 4 generators, 1 automorphisms, 2 orbits.
Graph 3 has 4 generators, 7 automorphisms, 2 orbits.
Graph 4 has 4 generators, 1 automorphisms, 2 orbits.
Graph 4 has 4 generators, 1 automorphisms, 2 orbits.
Done calculating edge automorphisms in 1ms.

Starting main algorithm...
Starting graph 1/5
Starting graph 2/5
Progress: 10/30 combinations done in 0ms.
Starting graph 3/5
Starting graph 4/5
Progress: 20/30 combinations done in 0ms.
Starting graph 5/5
(k,g)-graph - 3 5 10 - IsP@OkWHG
egr-graph - 3 5 4 10 - IsP@OkWHG
vgr-graph - 3 5 6 10 - IsP@OkWHG
Progress: 30/30 combinations done in 0ms.
Done in 0ms.

Total runtime: 1ms.
```

The following types of lines are written to standard output:

- `(k,g)-graph - <k> <g> <n> - graph6`: A $(k,g)$-graph with $n$ vertices and the given graph6 encoding.
- `egr-graph - <k> <g> <lambda_e> <n> - graph6`: An $(k,g,\lambda_e)$-graph with $n$ vertices and the given graph6 encoding.
- `vgr-graph - <k> <g> <lambda_v> <n> - graph6`: A $(k,g,\lambda_v)$-graph with $n$ vertices and the given graph6 encoding.
- `(k,g,g+1)-graph - <k> <g> <n> - graph6`: A $(k,g,\underline{g+1})$-graph with $n$ vertices and the given graph6 encoding.

The remaining text is written to standard error. 
The output is filtered up to isomorphism, so the same graph will be written at most 4 times (once for each type of graph).
Note that the invocation above correctly finds the Petersen graph as a lift of the dumbbell graph.