#include "CostFunctions.h"

#include <iostream>
#include <queue>
#include <stack>
#include <unordered_set>

#include "Other.h"
#include "RNG.h"

double averageSimplifiedShortestCycles(Graph& graph, const Group& group, const int amountOfCycles, const vector<int>& weights) {
    if (weights.size() != amountOfCycles && !weights.empty()) {
        throw invalid_argument("The amount of weights should be equal to requested amount of cycles.");
    }

    Graph lifted = lift(graph, group);
    vector<int> cycleLengths;
    int largestIndex = -1;

    // Verify that the graph is connected
    if (!lifted.isConnected()) {
        return -1;
    }

    for (int v = 0; v < size(lifted.adjacency); v++) {
        for (int u : lifted.adjacency[v]) {
            if (u <= v) {
                continue;
            }

            int dist;
            try {
                dist = lifted.distance(v, u, false).first + 1;
            } catch (const invalid_argument& e) {
                continue;
            }

            if (size(cycleLengths) < amountOfCycles) {
                cycleLengths.push_back(dist);
                continue;
            }

            if (largestIndex == -1) {
                largestIndex = 0;
                for (int i = 1; i < amountOfCycles; i++) {
                    if (cycleLengths[i] > cycleLengths[largestIndex]) {
                        largestIndex = i;
                    }
                }
            }

            // replace largest if new cycle is smaller
            if (dist < cycleLengths[largestIndex]) {
                cycleLengths[largestIndex] = dist;
                largestIndex = -1;
            }

        }
    }

    int result = 0;
    int weightSum = 0;
    for (int i = 0; i < size(cycleLengths); i++) {
        int weight = 1;
        if (!weights.empty()) {
            weight = weights[i];
        }

        result += cycleLengths[i] * weight;
        weightSum += weight;
    }

    return static_cast<double>(result) / weightSum;
}

vector<vector<int>> getFundamentalCycles(Graph& graph){
    // const auto& edgesNotInTree = graph.getEdgesWithoutDefaultVoltage();

    // Create spanning tree
    set<int> connected;
    set<int> explored;
    set<int> unexplored;
    set<pair<int, int>> notInSpanningTree;
    vector<vector<int>> treeAdjacency(size(graph.adjacency));
    for (int i = 0; i < size(graph.adjacency); i++) {
        treeAdjacency[i] = vector<int>();
    }

    unexplored.insert(0);
    while (size(explored) != size(graph.adjacency)) {
        int current = *unexplored.begin();
        unexplored.erase(current);
        explored.insert(current);
        connected.insert(current);

        for (auto v: graph.adjacency[current]) {
            if (!connected.contains(v)) {
                connected.insert(v);
                unexplored.insert(v);

                treeAdjacency[current].push_back(v);
                treeAdjacency[v].push_back(current);
            } else {
                pair edge = {current, v};
                notInSpanningTree.insert(edge);

                pair back = {v, current};
                notInSpanningTree.insert(back);
            }
        }
    }


    // Construct Fundamental Cycles
    vector<vector<int>> fundamentalCycles;

    fundamentalCycles.reserve(notInSpanningTree.size());
    for (const auto& edge: notInSpanningTree) {
        // Only consider i -> j if i < j
        if (edge.first > edge.second) {
            continue;
        }

        // Do a BFS from the vertex
        vector distances(size(graph.adjacency), -1);
        vector parent(size(graph.adjacency), -1);
        queue<int> queue;

        distances[edge.first] = 0;
        queue.push(edge.first);

        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();

            for (int neighbour: treeAdjacency[current]) {
                if (current == edge.first && neighbour == edge.second) {
                    continue;
                }

                if (distances[neighbour] == -1) {
                    distances[neighbour] = distances[current] + 1;
                    parent[neighbour] = current;
                    queue.push(neighbour);
                }

                // Stop if we reached the second vertex
                if (neighbour == edge.second) {
                    vector<int> path;
                    for (int v = edge.second; v != -1; v = parent[v]) {
                        path.push_back(v);
                    }
                    fundamentalCycles.push_back(path);
                }
            }
        }
    }

    return fundamentalCycles;
}


set<set<pair<int, int>>> getCycles(Graph& graph){
    auto fundamentalCycles = getFundamentalCycles(graph);

    if (fundamentalCycles.empty()) {
        return {};
    }

    // Convert fundamental cycles to vector<set<pair<int, int>>>
    vector<set<pair<int, int>>> fundamentalCyclesSet;
    for (const auto& cycle: fundamentalCycles) {
        set<pair<int, int>> cycleSet;
        for (int i = 0; i < size(cycle); i++) {
            int start = cycle[i];
            int end = cycle[(i + 1) % size(cycle)];
            if (start > end) {
                swap(start, end);
            }
            cycleSet.insert({start, end});
        }
        fundamentalCyclesSet.push_back(cycleSet);
    }

    set<set<pair<int, int>>> S;
    S.insert(fundamentalCyclesSet.front());
    set<set<pair<int, int>>> Q;
    Q.insert(fundamentalCyclesSet.front());
    set<set<pair<int, int>>> R;
    set<set<pair<int, int>>> RStar;
    for (int i = 1; i < size(fundamentalCyclesSet); i++) {
        auto& Bi = fundamentalCyclesSet[i];

        for (auto& T : Q){
            bool intersectionIsEmpty = true;
            for (auto& t : T) {
                if (Bi.contains(t)) {
                    intersectionIsEmpty = false;
                    break;
                }
            }
            set<pair<int, int>> xorSet = T;
            for (const auto& b : Bi) {
                if (xorSet.contains(b)) {
                    xorSet.erase(b);
                } else {
                    xorSet.insert(b);
                }
            }


            if (intersectionIsEmpty) {
                RStar.insert(xorSet);
            } else {
                R.insert(xorSet);
            }
        }

        set<set<pair<int, int>>> RNew;
        for (auto& r : R) {
            RNew.insert(r);
        }

        for (auto& U : R) {
            for (auto& V : R) {
                if (U == V) {
                    continue;
                }
                bool isSubset = true;
                for (auto& u : U) {
                    if (!V.contains(u)) {
                        isSubset = false;
                        break;
                    }
                }
                if (isSubset) {
                    RNew.erase(V);
                    RStar.insert(V);
                }
            }
        }

        // S = S U R U {Bi}
        for (auto& r : RNew) {
            S.insert(r);
        }
        S.insert(Bi);

        // Q = Q U R U R* U {Bi}
        for (auto& r : RNew) {
            Q.insert(r);
        }
        for (auto& r : RStar) {
            Q.insert(r);
        }
        Q.insert(Bi);

        // Reset R and R*
        R.clear();
        RStar.clear();
    }

    // // print S
    // for (const auto& cycle : S) {
    //     for (const auto& edge : cycle) {
    //         clog << edge.first << " -> " << edge.second << " ";
    //     }
    //     clog << endl;
    // }

    return S;
}

int getRandomCycleLength(Graph& graph) {
    // This cannot find triangles
    // Construct 1 random cycle using DFS
    set<int> visited;
    vector parents(size(graph.adjacency), -1);
    stack<int> stack;

    int current = getRandomInt(0, size(graph.adjacency) - 1);
    stack.push(current);
    parents[current] = -1;

    while (!stack.empty()) {
        current = stack.top();
        stack.pop();

        if (!visited.contains(current)) {
            visited.insert(current);

            vector neighbours = graph.adjacency[current];
            shuffleVec(neighbours);

            for (int neighbour : neighbours) {
                if (!visited.contains(neighbour)) {
                    stack.push(neighbour);
                    parents[neighbour] = current;
                } else if (parents[current] != neighbour) {
                    // Found a cycle
                    int length = 1;
                    int temp = current;
                    while (temp != neighbour) {
                        temp = parents[temp];
                        length++;
                    }
                    return length+1;
                }

            }
        }
    }

    return -1;
}

double monteCarloAverageCycleLength(Graph& graph, const Group& group, const int N) {
    Graph lifted = lift(graph, group);

    // Verify that the graph is connected
    if (!lifted.isConnected()) {
        return -1;
    }

    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += getRandomCycleLength(lifted);
    }

    return sum / N;
}

double averageCycleLength(Graph& graph, const Group& group, const int ignoreCount) {
    Graph lifted = lift(graph, group);

    // Verify that the graph is connected
    if (!lifted.isConnected()) {
        return -1;
    }

    const auto cycles = getCycles(lifted);

    map<size_t, int> cycleLengths;
    for (auto& cycle: cycles) {
        cycleLengths[cycle.size()]++;
    }

    // Average ignoring the smallest ignoreCount cycles
    vector<int> lengths;
    lengths.reserve(cycleLengths.size());
    for (const auto& cycleLength: cycleLengths) {
        lengths.push_back(cycleLength.first);
    }
    sort(lengths.begin(), lengths.end());

    int sum = 0;
    int ignored = 0;
    for (const int length : lengths) {
        if (ignored < ignoreCount) {
            ignored += cycleLengths[length];
            if (ignored > ignoreCount) {
                sum += (ignored - ignoreCount) * length;
            }
        } else {
            sum += cycleLengths[length] * length;
        }

    }

    return static_cast<double>(sum) / (cycles.size() - ignoreCount);
}

double averageFundamentalCycleLength(Graph &graph, const Group &group, int ignoreCount) {
    Graph lifted = lift(graph, group);

    // Verify that the graph is connected
    if (!lifted.isConnected()) {
        return -1;
    }

    const auto cycles = getFundamentalCycles(lifted);

    vector<int> cycleLengths;
    cycleLengths.reserve(cycles.size());
    for (const auto& cycle: cycles) {
        cycleLengths.push_back(size(cycle));
    }
    sort(cycleLengths.begin(), cycleLengths.end());

    double sum = 0;
    int count = 0;
    for (const auto& l: cycleLengths) {
        if (count >= ignoreCount) {
            sum += l;
        }
        count++;
    }

    return sum / (count - ignoreCount);
}


double voltageDiversity(Graph &graph, const Group &group) {
    vector<shared_ptr<Edge>> edges = graph.getEdgesWithoutDefaultVoltage();
    vector<shared_ptr<Edge>> legalEdges = filterInverses(edges);

    int numIdentities = 0;
    int numInverses = 0;
    int numRepeatedOrders = 0;
    int numRepeatedVoltages = 0;
    unordered_set<int> uniqueOrders;
    unordered_set<int> seenVoltages;
    for (const auto& edge : legalEdges) {
        int order = group.powerToIdentity[edge->voltage];
        if (order == 1) {
            numIdentities++;
        } else if (seenVoltages.contains(group.inverse[edge->voltage])) {
            numInverses++;
        } else if (seenVoltages.contains(edge->voltage)) {
            numRepeatedVoltages++;
        } else if (order > 1) {
            if (uniqueOrders.contains(order)) {
                numRepeatedOrders++;
            } else {
                uniqueOrders.insert(order);
            }
        }
        seenVoltages.insert(edge->voltage);
    }

    return -20 * numIdentities - 10 * numInverses - 10 * numRepeatedVoltages - numRepeatedOrders;
}

int countMaxLenPaths(int current, int target, int depth, int maxLength, const shared_ptr<Edge>& lastTaken, int currentVoltage, Graph& graph, const Group& group) {
    if (depth >= maxLength) return 0;

    // This skips paths consisting of a single loop or semi edge, but those should be caught by the static checks
    int amountFound = 0;
    for (const auto& [neighbor, edgesToNeighbor] : graph.neighbourToEdge[current]) {
        for (const auto& edge : edgesToNeighbor) {
            // Skip if voltage not set
            if (edge->voltage == -1) continue;

            // Dont take the reverse edge of the previously taken edge
            if (edge->reverseEdge == lastTaken) continue;

            // Semi edges are their own inverse
            if (edge->start == edge->end && edge->reverseEdge == nullptr && edge == lastTaken) continue;

            currentVoltage = group.multiplicationTable[currentVoltage][edge->voltage];
            if (neighbor == target && currentVoltage == 0 && depth == maxLength - 1) {
                amountFound++;
            }

            amountFound += countMaxLenPaths(neighbor, target, depth + 1, maxLength, edge, currentVoltage, graph, group);
            currentVoltage = group.multiplicationTable[currentVoltage][group.inverse[edge->voltage]];
        }
    }
    return amountFound;
}

double walkRegularity(Graph &graph, const Group &group) {
    vector vertex_regularity(graph.adjacency.size(), 0);
    vector edge_regularity(graph.edges.size(), 0);

    int girth = lift(graph, group).getGirth();
    for (int i = 0; i < size(graph.edges); i++) {
        auto& edge = graph.edges[i];
        if (edge->reverseEdge == nullptr || edge->start > edge->end) {
            continue;
        }

        int count = countMaxLenPaths(edge->start, edge->end, 1, girth, edge->reverseEdge, 0, graph, group);
        vertex_regularity[edge->start] += count;
        vertex_regularity[edge->end] += count;
        edge_regularity[i] += count;
    }

    double varVertexRegularity = 0;
    double meanVertexRegularity = 0;
    for (const int r : vertex_regularity) {
        meanVertexRegularity += r;
    }
    meanVertexRegularity /= vertex_regularity.size();
    for (const int r : vertex_regularity) {
        varVertexRegularity += (r - meanVertexRegularity) * (r - meanVertexRegularity);
    }
    varVertexRegularity /= vertex_regularity.size();

    double varEdgeRegularity = 0;
    double meanEdgeRegularity = 0;
    for (const int r : edge_regularity) {
        meanEdgeRegularity += r;
    }
    meanEdgeRegularity /= edge_regularity.size();
    for (const int r : edge_regularity) {
        varEdgeRegularity += (r - meanEdgeRegularity) * (r - meanEdgeRegularity);
    }
    varEdgeRegularity /= edge_regularity.size();
    return -min(varVertexRegularity, varEdgeRegularity);
}

pair<int, int> sampleWalk(int current, int target, int depth, int maxLength, const shared_ptr<Edge>& lastTaken, int currentVoltage, Graph& graph, const Group& group) {
    if (depth >= maxLength) return {-1, -1};

    auto neighbours = graph.adjacency[current];
    shuffleVec(neighbours);

    for (const auto & neighbor : neighbours) {
        auto& edgesToNeighbor = graph.neighbourToEdge[current][neighbor];
        shuffleVec(edgesToNeighbor);
        for (const auto& edge : edgesToNeighbor) {
            // Skip if voltage not set
            if (edge->voltage == -1) continue;

            // Dont take the reverse edge of the previously taken edge
            if (edge->reverseEdge == lastTaken) continue;

            // Semi edges are their own inverse
            if (edge->start == edge->end && edge->reverseEdge == nullptr && edge == lastTaken) continue;

            currentVoltage = group.multiplicationTable[currentVoltage][edge->voltage];
            if (neighbor == target && depth >= maxLength - 2) {
                return {depth, currentVoltage};
            }

            auto p = sampleWalk(neighbor, target, depth + 1, maxLength, edge, currentVoltage, graph, group);
            if (p.first != -1) {
                return p;
            }

            currentVoltage = group.multiplicationTable[currentVoltage][group.inverse[edge->voltage]];
        }
    }
    return {-1, -1};
}

double walkSampler(const Group& group, Graph& graph, const shared_ptr<Edge>& edge, const int minGirth, const int N) {
    auto& edges = graph.getEdgesWithoutDefaultVoltage();
    auto& e = edges[getRandomInt(0, size(edges) - 1)];
    shared_ptr<Edge> takenEdge = e->reverseEdge;
    if (takenEdge == nullptr) {
        takenEdge = e;
    }
    double totalCost = 0;
    for (int i = 0; i < N; i++) {
        auto p = sampleWalk(e->start, e->end, 1, minGirth-1, takenEdge, takenEdge->voltage, graph, group);
        if (p.first == -1) {
            clog << "Failed to sample walk" << endl;
            continue;
        }
        if (p.second == 0) {
            // Net voltage is 0
            totalCost += 1000;
        } else {
            totalCost += 1. / (p.first * group.powerToIdentity[p.second]);
        }
    }
    return -totalCost;
}
