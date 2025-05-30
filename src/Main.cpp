#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <variant>
#include <regex>
#include <unordered_set>
#include <thread>

#include "group/Group.h"
#include "graph/Graph.h"
#include "graph/GraphParser.h"
#include "algorithms/BTA.h"
#include "algorithms/TabuSearch.h"
#include "algorithms/Helpers/Other.h"
#include "group/GroupParser.h"

void run(
    const vector<Graph>& graphs,
    const vector<Group>& groups,
    runSetup &setup
    ) {

    clog << "Starting main algorithm..." << endl;
    const auto start = chrono::high_resolution_clock::now();

    int totalCombinationCount = graphs.size() * groups.size();
    int currentCombination = 0;
    int currentGraph = 0;
    for (auto graph : graphs) {
        clog << "Starting graph " << currentGraph+1<< "/" << graphs.size() << endl;
        currentGraph++;

        for (const auto& group : groups) {
            switch (setup.algorithm) {
                case runSetup::AlgorithmType::BTA: {
                    bool finished = filteredBTA(graph, group, setup);
                    if (!finished) {
                        clog << "BTA did not finish at group " << group.name << endl;
                    }
                    break;
                }
                case runSetup::AlgorithmType::TabuSearch: {
                    tabuSearch(graph, group, setup);
                    break;
                }
            }

            ++currentCombination;
            if (currentCombination % setup.ms.logEvery == 0) {
                auto end = chrono::high_resolution_clock::now();
                clog << "Progress: " << currentCombination << "/" << totalCombinationCount
                    << " combinations done in " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms."
                    << endl;
            }
        }
    }
    auto end = chrono::high_resolution_clock::now();
    clog << "Done in " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms." << endl << endl;
}

vector<Group> getGroups(int groupMin, int groupMax) {
    clog << "\nReading groups..." << endl;
    auto start = chrono::high_resolution_clock::now();

    string path;
    if (groupMax <= 50) {
        path = "./data/groups/groups_1-50.txt";
    } else {
        cerr << "Group size too large." << endl;
        exit(1);
    }

    ifstream groupPipe(path);
    vector<Group> result;
    string line;
    string groupLine;

    while (getline(groupPipe, line)) {
        if (line == "####") {
            Group group = fromString(groupLine);
            if (group.multiplicationTable.size() < groupMin) {
                groupLine = "";
                continue;
            }
            if (group.multiplicationTable.size() > groupMax) {
                break;
            }
            result.push_back(group);

            groupLine = "";
        } else {
            groupLine += line;
        }
    }
    groupPipe.close();

    auto end = chrono::high_resolution_clock::now();
    clog << "Found " << result.size() << " groups in " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms." << endl;

    return result;
}

void parseDreadnaut(vector<Graph>& graphs, const int partition) {
    clog << "\nCalculating usable edge automorphisms..." << endl;
    auto start = chrono::high_resolution_clock::now();

    ifstream orbitPipe("/tmp/pipe_" + to_string(partition) + "_out");
    string line;
    string orbitLine;
    string genLine;
    vector<vector<int>> orbits;
    vector<string> generatorLines;
    vector<string> currentGenerators;
    int orbitIndex = 0;

    while (getline(orbitPipe, line)) {
        // Skip
        if (line.starts_with("level")) {
            continue;
        }

        // Continuation of generator
        if (line.starts_with("   (")) {
            genLine += line.substr(3);
            continue;
        }

        // New generator
        if (line.starts_with("(")) {
            if (!genLine.empty()) {
                generatorLines.push_back(genLine);
                genLine.clear();
            }
            genLine = line;
            continue;
        }

        // Orbits will start
        if (line.starts_with("tctotal=")) {
            // Parse the last generator line
            if (!genLine.empty()) {
                generatorLines.push_back(genLine);
                genLine.clear();
            }

            while (getline(orbitPipe, line)) {
                if (line.starts_with("(") || line.contains("orbits")) {
                    orbits = parseOrbit(orbitLine);
                    graphs[orbitIndex].setOrbits(orbits);

                    auto automorphisms = getAutomorphisms(generatorLines, graphs[orbitIndex], 2000);
                    clog << "Graph " << orbitIndex+1 << " has " << generatorLines.size()
                            << " generators, " << automorphisms.size()
                            << " automorphisms, " << graphs[orbitIndex].orbits.size() << " orbits." << endl;
                    graphs[orbitIndex].edgeAutomorphisms = automorphisms;

                    orbitIndex++;
                    orbits.clear();
                    generatorLines.clear();
                    orbitLine = "";
                    if (line.starts_with("(")) {
                        genLine = line;
                    }
                    break;
                }
                orbitLine += line;
            }
            if (!(line.starts_with("(") || line.contains("orbits"))) {
                orbits = parseOrbit(orbitLine);
                graphs[orbitIndex].setOrbits(orbits);

                auto automorphisms = getAutomorphisms(generatorLines, graphs[orbitIndex], 1000);
                clog << "Graph " << orbitIndex << " has " << generatorLines.size()
                        << " generators, " << automorphisms.size()
                        << " automorphisms, " << graphs[orbitIndex].orbits.size() << " orbits." << endl;
                graphs[orbitIndex].edgeAutomorphisms = automorphisms;
                break;
            }
        }
    }

    orbitPipe.close();

    auto end = chrono::high_resolution_clock::now();
    clog << "Done calculating edge automorphisms in " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms." << endl << endl;
}

void writeToDreadnaut(const vector<Graph>& graphs, const int partition) {
    auto start = chrono::high_resolution_clock::now();

    ofstream dreadnautPipe("/tmp/pipe_" + to_string(partition) + "_in");
    for (Graph graph : graphs) {
        vector<vector<int>> lineGraphAdjacency = graph.getLineGraph();
        string s = "+d n=" + to_string(lineGraphAdjacency.size()) + " $=0 g\n";
        for (int i = 0; i < lineGraphAdjacency.size(); i++) {
            s += to_string(i) + " : ";
            for (int j : lineGraphAdjacency[i]) {
                s += to_string(j) + " ";
            }
            s.pop_back();
            s += ";\n";
        }
        s += "$$xo\n";
        dreadnautPipe << s;
    }
    dreadnautPipe << "q" << endl;
    dreadnautPipe.close();

    auto end = chrono::high_resolution_clock::now();
}

vector<vector<int>> encodeMultigraphToSimpleGraph(const vector<vector<int>>& multigraph) {
    int n = multigraph.size();
    vector<vector<int>> simpleGraph;
    simpleGraph.resize(n);

    int nextNode = n;

    for (int u = 0; u < n; ++u) {
        unordered_map<int, int> edgeCount;

        for (int v : multigraph[u]) {
            edgeCount[v]++;
        }

        for (auto& [v, count] : edgeCount) {
            if (u == v) {
                // Each is replaced with a dummy node attachment
                for (int i = 0; i < count; ++i) {
                    simpleGraph.emplace_back();
                    int dummy = nextNode++;
                    simpleGraph[u].push_back(dummy);
                    simpleGraph[dummy].push_back(u);
                }
                continue;
            }
            if (u >= v) {
                continue;
            }

            // For each occurrence of u-v, add a dummy node u-dummy-v
            for (int i = 0; i < count; ++i) {
                simpleGraph.emplace_back();
                int dummy = nextNode++;
                simpleGraph[u].push_back(dummy);
                simpleGraph[v].push_back(dummy);
                simpleGraph[dummy].push_back(u);
                simpleGraph[dummy].push_back(v);
            }

        }
    }

    return simpleGraph;
}

string computeCanonicalLabelMultigraph(const Graph& graph) {
    // Convert the multigraph to a simple graph
    vector<vector<int>> adj = encodeMultigraphToSimpleGraph(graph.adjacency);

    Graph simpleGraph;
    simpleGraph.adjacency = adj;
    return simpleGraph.getGraph6(true);
}

void addLoops(const int k, vector<Graph>& result, const bool noSemiEdges) {
    // Add all loops in a single go
    if (noSemiEdges) {
        for (auto& graph : result) {
            for (int v = 0; v < size(graph.adjacency); v++) {
                int missing = (k - size(graph.adjacency[v]))/2;
                if (missing <= 0) continue;

                // Fill with loops
                for (int i = 0; i < missing; i++) {
                    graph.adjacency[v].push_back(v);
                }
            }
        }
        return;
    }

    // Add one loop at a time and keep each resulting graph
    bool changed = true;
    int indexInResult = 0;
    unordered_set<string> g6s;
    while (changed) {
        changed = false;

        vector<Graph> graphsToAdd;
        for (int i = indexInResult; i < size(result); i++) {
            auto& graph = result[i];
            for (int v = 0; v < size(graph.adjacency); v++) {
                int missing = k;
                for (int element : graph.adjacency[v]) {
                    missing--;
                    if (element == v) {
                        missing--;
                    }
                }
                if (missing <= 1) continue;

                Graph g;
                vector<vector<int>> adjacency = graph.adjacency;
                g.adjacency = adjacency;
                g.adjacency[v].push_back(v);

                // Skip isomorphic graphs
                string g6 = computeCanonicalLabelMultigraph(g);
                if (g6s.contains(g6)) {
                    continue;
                }
                g6s.insert(g6);

                graphsToAdd.push_back(g);
                changed = true;
            }
        }

        indexInResult = size(result);
        for (auto& g : graphsToAdd) {
            result.push_back(g);
        }
    }
}

bool willHaveSemiEdges(const int k, const Graph& graph) {
    for (auto& adj : graph.adjacency) {
        int missing = k;
        for (int i = 0; i < size(adj); i++) {
            missing--;
            if (adj[i] == i) {
                missing--;
            }
        }
    }
    return false;
}

void getMultigraphBaseGraphs(const int k, const int n, const int partition, const int modulo, const int used, const int index,
                             vector<int>& ns, vector<Graph>& result, int& count, const bool noSemiEdges) {
    if (index == k) {
        if (used == n) {

            // Odd ns should be 0 when avoiding semi edges
            if (noSemiEdges) {
                for (int i = 0; i < k; i++) {
                    if (ns[i] > 0 && i % 2 == k % 2) {
                        return;
                    }
                }
            }

            int initialSize = size(result);

            vector<Graph> baseMultigraphs = parseMultiGraphs(ns);
            clog << "Found " << baseMultigraphs.size() << " base multigraphs for degree distribution: ";
            for (int i = 0; i < k; i++) {
                clog << ns[i] << " ";
            }
            clog << "\n";

            for (auto& graph : baseMultigraphs) {
                vector<Graph> extended;
                extended.push_back(graph);

                addLoops(k, extended, noSemiEdges);

                for (auto& g : extended) {
                    count++;
                    // Partitioning
                    if (count % modulo != partition) {
                        continue;
                    }

                    // Semi edge filter
                    if (noSemiEdges && willHaveSemiEdges(k, g)) {
                        continue;
                    }

                    result.push_back(g);
                }
            }

            clog << "Resulting in " << size(result) - initialSize << " new multigraphs for this partition." << endl;
        }
        return;
    }

    for (int i = 0; i <= n - used; i++) {
        ns[index] = i;
        getMultigraphBaseGraphs(k, n, partition, modulo, used + i, index + 1, ns, result, count, noSemiEdges);
    }
    ns[index] = 0;
}

void fillWithSemiEdges(const int k, int& semiEdgeGraphs, vector<Graph>& graphs)
{
    for (auto& graph: graphs) {
        bool graphChanged = false;
        graph.createEdgesWithDefaultVoltage();
        // check if graph needs semi edges
        for (int v = 0; v < size(graph.adjacency); v++) {
            int missing = k - size(graph.adjacency[v]);
            if (missing == 0) continue;

            // Fill with semi edges
            for (int i = 0; i < missing; i++) {
                graph.adjacency[v].push_back(v);
                auto e = make_shared<Edge>(v, v);
                graph.edges.push_back(e);
                graph.neighbourToEdge[v][v].push_back(e);
                graph.edgesWithoutSpanningTree.push_back(e);
                if (!graphChanged) {
                    semiEdgeGraphs++;
                    graphChanged = true;
                }
            }
        }
        // Verify the graph is k regular
        for (int v = 0; v < size(graph.adjacency); v++) {
            if (size(graph.adjacency[v]) != k) {
                cerr << "Graph is not k regular!" << endl;
                graph.printAdjacency();
                exit(1);
            }
        }
    }

    // Set the largest amount of semi edges at any vertex for each graph
    for (auto& graph : graphs) {
        for (int v = 0; v < size(graph.adjacency); v++) {
            int semi_edges = 0;
            for (const auto& edge : graph.neighbourToEdge[v][v]) {
                if (edge->reverseEdge == nullptr) {
                    semi_edges++;
                }
            }
            graph.maxSemiEdgesPerVertex = max(graph.maxSemiEdgesPerVertex, semi_edges);
        }
    }
}

vector<Graph> getMultiGraphs(const int k, const int n, const int partition, const int modulo, const bool noSemiEdges) {
    clog << "\nConstructing multigraphs..." << endl;
    auto start = chrono::high_resolution_clock::now();
    int semiEdgeGraphs = 0;
    int currentNumber = -1;

    vector<int> ns;
    ns.reserve(k);
    for (int i = 0; i < k; i++) {
        ns.push_back(0);
    }

    vector<Graph> result;
    getMultigraphBaseGraphs(k, n, partition, modulo, 0, 0, ns, result, currentNumber, noSemiEdges);
    fillWithSemiEdges(k, semiEdgeGraphs, result);

    if (noSemiEdges) {
        vector<Graph> noSemiEdgeGraphs;
        for (auto& graph : result) {
            if (graph.maxSemiEdgesPerVertex == 0) {
                noSemiEdgeGraphs.push_back(graph);
            }
        }
        result = noSemiEdgeGraphs;
    }

    auto end = chrono::high_resolution_clock::now();
    clog << "Found " << result.size() << " " << k << "-regular multigraphs with " << n << " vertices (" << result.size() - semiEdgeGraphs <<" without semi edges) in "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms\n" << endl;

    if (result.size() <= 1000) {
        int gc = 1;
        for (auto& graph : result) {
            clog << "Graph " << gc << " with at most " << graph.maxSemiEdgesPerVertex << " semi edges at any vertex." << endl;
            graph.printAdjacency();
            gc++;
        }
    }

    return result;
}

int main(int argc, char *argv[] ) {
    if (argc != 8) {
        cerr << "Usage: <partition> <modulo> <n> <k> <groupMin> <groupMax> <minGirth>" << endl;
        return 1;
    }

    const int partition = atoi(argv[1]);
    const int modulo = atoi(argv[2]);
    const int n = atoi(argv[3]);
    const int k = atoi(argv[4]);
    const int groupMin = atoi(argv[5]);
    const int groupMax = atoi(argv[6]);
    const int minGirth = atoi(argv[7]);
    clog << "Partition: " << partition << endl;
    clog << "Modulo: " << modulo << endl;
    clog << "n: " << n << endl;
    clog << "k: " << k << endl;
    clog << "groupMin: " << groupMin << endl;
    clog << "groupMax: " << groupMax << endl;
    clog << "minGirth: " << minGirth << endl;

    mainSetup main_setup = {
        .seed = 0,
        .k = k,
        .n = n,
        .minGirth = max(minGirth, 3),
        .logEvery = 10,
        .useGraphAutomorphisms = true,
        .useGroupAutomorphisms = true
    };

    auto setup_BTA = runSetup(
        main_setup,
        runSetup::AlgorithmType::BTA,
        bruteForceSetup{
            .timeLimit = 200000
        }
    );

    auto setup_TS = runSetup(
        main_setup,
        runSetup::AlgorithmType::TabuSearch,
        tabuSearchSetup{
            .timeLimit = 20,
            .maxIterations = 100000,
            .tabuSizeMult = 3,
            .PerturbAfterNoImprovement = 100,
            .neighbourMinGirth = minGirth,
            .costFunction = createWalkSamplerCostFunction(500, main_setup.minGirth)
            // .costFunction = createWalkRegularityCostFunction()
        }
    );

    auto start = chrono::high_resolution_clock::now();

    bool noSemiEdges = false;
    vector<Graph> graphs = getMultiGraphs(main_setup.k, main_setup.n, partition, modulo, noSemiEdges);

    vector<Group> groups = getGroups(groupMin, groupMax);

    if (graphs.size() <= 1000) {
        std::thread parser(parseDreadnaut, std::ref(graphs), partition);
        writeToDreadnaut(graphs, partition);
        parser.join();
    } else {
        clog << "Skipping automorphism calculation, too many graphs." << endl;
    }

    run(graphs, groups, setup_BTA);
    auto end = chrono::high_resolution_clock::now();
    clog << "Total runtime: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms." << endl;

    return 0;
}