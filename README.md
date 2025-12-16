# New small regular graphs of given girth: the cage problem and beyond

This repository contains the code and data related to the Master's dissertation "Algoritmes voor reguliere lifts van grafen: het cage probleem en varianten" available at [Limo Libis](https://kuleuven.limo.libis.be/discovery/fulldisplay?docid=alma9995771928401488&context=L&vid=32KUL_KUL:KULeuven&search_scope=All_Content&tab=all_content_tab&lang=en). For the repository about the manuscript "New small regular graphs of given girth: the cage problem and beyond", please visit [The following GitHub Repository](https://github.com/AGT-Kulak/smallRegGirthGraphs). 
Code can be found in the `src` directory, while data is available in the `data` directory.

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

The groups subdirectory contains all non-isomorphic groups of order 1 to 50 with up to 1000 automorphisms for each group.

The graphs subdirectory contains all graphs stored in graph6 format mentioned in the manuscript.

They are organized by type and use the following naming conventions:
- `cage`: Small $(k,g)$-graphs where `nW_kX_gY.g6` denotes a file containing $(X,Y)$-graphs with $W$ vertices.
- `egr`: Small $(k,g,\lambda_e)$-graphs where `nW_kX_gY_lZ.g6` denotes a file containing $(X,Y,Z)$-graphs with $W$ vertices.
- `kgnogpo`: Small $(k,g,\underline{g+1})$-graphs where `nW_kX_gY.g6` denotes a file containing $(X,Y,\underline{Y+1})$-graphs with $W$ vertices.
- `spectra.zip`: Zipped to preserve space. Contains $(k,g)$-graphs where `nW_kX_gY.g6` denotes a file containing $(X,Y)$-graphs with $W$ vertices.
- `vgr.zip`: Zipped to preserve space. Contains small $(k,g,\lambda_v)$-graphs where `nW_kX_gY_lZ.g6` denotes a file containing $(X,Y,Z)$-graphs with $W$ vertices.

## Code
The `src` directory contains multiple files and subdirectories. 
The `src/GAP` subdirectory contains GAP code that can generate the files containing groups and automorphisms used to construct lifts.
The `src/nautyAndMultigraph` subdirectory contains slightly modified versions of the Nauty library and the generator multigraph.

### Compilation

The remaining code revolves around constructing lifts to obtain $k$-regular graphs achieving a given minimum girth. 
It can be compiled from the top-level directory using the following commands (tested with GCC 11.4):
```bash
cmake CMakeLists.txt
make
```

### Usage
The executables are available under the `bin` directory, but should be executed using:
```bash
. ./run.sh <partition> <modulo> <n> <k> <groupMin> <groupMax> <minGirth> <pureGraph6> <Optional:groupfile>
```
This will construct all non-isomorphic connected lifts achieving the given `minGirth` using all possible $k$-regular multigraphs on $n$ vertices with loops and semi-edges as base graph
and groups with order between `groupMin` and `groupMax` (inclusive). The parameters `partition` and `modulo` can be used to evenly split the base graphs over multiple invocations.
The `pureGraph6` parameter can be set to 1 to only output the graph6 encoding of the resulting graphs, or 0 to also output additional information for each found graph.
The `groupfile` parameter can be used to specify a file containing groups to use instead of the default groups.

For example, to construct all lifts with girth at least 5 using 3-regular base graphs on 2 vertices 
and groups with up to 5 elements, you can run the following command:
```bash 
. ./run.sh 0 1 2 3 1 5 5 0
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
pureGraph6: 1
groupFile: ./data/groups/groups_1-50.txt

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

The following types of lines are written to cout:

- `(k,g)-graph - <k> <g> <n> - graph6`: A $(k,g)$-graph with $n$ vertices and the given graph6 encoding.
- `egr-graph - <k> <g> <lambda_e> <n> - graph6`: An $(k,g,\lambda_e)$-graph with $n$ vertices and the given graph6 encoding.
- `vgr-graph - <k> <g> <lambda_v> <n> - graph6`: A $(k,g,\lambda_v)$-graph with $n$ vertices and the given graph6 encoding.
- `(k,g,g+1)-graph - <k> <g> <n> - graph6`: A $(k,g,\underline{g+1})$-graph with $n$ vertices and the given graph6 encoding.

The remaining text is written to clog. 
The output is filtered up to isomorphism, so the same graph will be written at most 4 times (once for each type).

If you are only interested in the actual graphs in graph6 format, you can obtain them as follows:
```bash 
. ./run.sh 0 1 2 3 1 5 5 1 > graphs.g6
```

This would write each graph in a separate line. The command above would result in a file containing:
```
IsP@OkWHG
```

### Using different groups
This repository contains a file with all small groups of order 1 to 50 with up to 1000 automorphisms for each group.
If you want to use a different set of groups, you can generate this using GAP.

Assuming GAP is installed at `/path/to/gap`, you can run the following command:
```bash
/path/to/gap -q ./src/GAP/SmallGroups.g > groups_1-50.txt
```

to generate the `data/groups/groups_1-50.txt` file.

If you are e.g. only interested in lifts from cyclic groups, you can create a file under `./src/GAP/CyclicGroups.g` with the following content:
```gap
Read("./src/GAP/Automorphisms.g");
Read("./src/GAP/Orbits.g");

lowerBound := 1;
upperBound := 50;
automorphismCount := 1000;

for order in [lowerBound..upperBound] do
  G := CyclicGroup(order);
  M := MultiplicationTable(G);

  Print("Order: ", order, ", Number: 1\n");
  Print(M);
  Print("\nA\n");
  PrintAutomorphismMappings(G, automorphismCount);
  Print("\nO\n");
  PrintElementOrbits(G);
  Print("####");
  Print("\n");
od;

quit;
```

The following command will then create a file with these groups:
```bash
/path/to/gap -q ./src/GAP/CyclicGroups.g > cyclic_groups_1-50.txt
```

To now construct all lifts with girth at least 5 using 3-regular base graphs on 2 vertices
and cyclic groups with up to 10 elements, you can run the following command:
```bash 
. ./run.sh 0 1 2 3 1 10 5 0 ./cyclic_groups_1-50.txt
```

Note that group files are assumed to contain groups sorted by their order.
