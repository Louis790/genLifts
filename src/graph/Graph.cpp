#include "Graph.h"
#include "GraphParser.h"
#include "../nautyAndMultigraph/labelg.h"

#include <array>
#include <bitset>
#include <climits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

using std::array;
using std::bitset;
using std::find;
using std::ifstream;
using std::invalid_argument;
using std::map;
using std::ofstream;
using std::queue;
using std::set;
using std::shared_ptr;
using std::size;
using std::string;
using std::to_string;
using std::vector;
using std::min;
using std::make_shared;

Graph::Graph() = default;

vector<int> & Graph::operator[](int i) {
    return adjacency[i];
}

void Graph::printAdjacency() {
    std::clog << "Adjacency:\n";
    for (int v = 0; v < size(adjacency); v++) {
        std::string s = " ";
        for (auto c : adjacency[v]) {
            s += std::to_string(c) + ", ";
        }
        std::clog << v << ":" << s << '\n';
    }
}

string Graph::getGraph6(bool makeCanonical) {
    // Modified from https://github.com/JorikJooken/edgeGirthRegularGraphs/blob/master/Code/generateRGLambdaGraphs.c
    const int n = size(adjacency);
    string graphString;
    graphString.reserve(8 + n*(n - 1)/2);

    if (n > 68719476735) {
        std::cerr << "Error: number of vertices too large.\n";
        std::exit(1);
    }

    //  Save number of vertices in the first one, four or 8 bytes.
    if(n <= 62) {
        graphString.push_back(static_cast<char>(n + 63));
    }
    else if(n <= 258047) {
        graphString.push_back(63 + 63);
        for (int i = 2; i >= 0; i--) {
            const char c = static_cast<char>(n >> i * 6 & 0x3F); // extract 6 bits
            graphString.push_back(static_cast<char>(c + 63));
        }
    }
    else {
        graphString.push_back(63 + 63);
        graphString.push_back(63 + 63);
        for (int i = 5; i >= 0; i--) {
            const char c = static_cast<char>(n >> i * 6 & 0x3F); // Extract 6 bits
            graphString.push_back(static_cast<char>(c + 63));
        }
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < n; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = static_cast<char>(charToPrint << 1);
            if (find(adjacency[i].begin(), adjacency[i].end(), j) != adjacency[i].end()) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString.push_back(static_cast<char>(charToPrint + 63));
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    // Verify every byte is in the range 63-126.
    for (char c: graphString) {
        if (c < 63 || c > 126) {
            std::cerr << "Error: Invalid character in graph6 string.\n";
            std::exit(1);
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = static_cast<char>(charToPrint << 1);
            if(++counter == 6) {
                graphString.push_back(static_cast<char>(charToPrint + 63));
            }
        }
    }

    if (makeCanonical) {
        #pragma omp critical
        {
            char* argv[] = {static_cast<char*>("-q")};
            char* resCode = runLabelG(1, argv, const_cast<char*>(graphString.c_str()));
            graphString = string(resCode);
        }
    }

    return graphString;
}

void Graph::saveAsMatrix(const string& path) {
    vector<vector<int>> adjacencyMatrix;
    const int n = size(adjacency);
    for (int i = 0; i < n; i++) {
        adjacencyMatrix.emplace_back();
        for (int j = 0; j < n; j++) {
            adjacencyMatrix[i].push_back(0);
        }
    }

    for (const auto& edge : edges) {
        adjacencyMatrix[edge->start][edge->end] = 1;
    }

    vector<string> lines;
    for (int i = 0; i < n; i++) {
        string line;
        for (int j = 0; j < n; j++) {
            line += to_string(adjacencyMatrix[i][j]) + " ";
        }
        line += "\n";
        lines.push_back(line);
    }

    // Write to file
    ofstream file(path);
    for (const auto& line : lines) {
        file << line;
    }
    file.close();
}

int Graph::getGirth() {
    int girth = INT_MAX;
    // Check for loops
    for (const auto& edge: edges) {
        if (edge->start == edge->end) {
            return 1;
        }
    }

    // Check for double edges
    set<int> neighbours;
    for (auto & neighs : adjacency) {
        for (int neighbour: neighs) {
            if (neighbours.contains(neighbour)) {
                return 2;
            }
            neighbours.insert(neighbour);
        }
        neighbours.clear();
    }

    for (int vertex = 0; vertex<size(adjacency); vertex++) {
        // Do a BFS from the vertex
        vector<int> distances(size(adjacency), -1);
        vector<int> parents(size(adjacency), -1);
        queue<int> queue;

        distances[vertex] = 0;
        queue.push(vertex);

        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();

            for (int neighbour: adjacency[current]) {
                if (distances[neighbour] == -1) {
                    distances[neighbour] = distances[current] + 1;
                    parents[neighbour] = current;
                    queue.push(neighbour);
                }
                if (parents[current] != neighbour && parents[neighbour] != current) {
                    girth = min(girth, distances[current] + distances[neighbour] + 1);
                }
            }
        }
    }

    return girth;
}

bool Graph::isOneEdgeConnected() {
    for (int i = 0; i<size(adjacency); i++) {
        for (int j : adjacency[i]) {
            if (i == j) {
                continue;
            }
            adjacency[i].erase(find(adjacency[i].begin(), adjacency[i].end(), j));
            adjacency[j].erase(find(adjacency[j].begin(), adjacency[j].end(), i));
            if (!isConnected()) {
                // Sort adjacency lists
                for (auto& vertex: adjacency) {
                    sort(vertex.begin(), vertex.end());
                }
                return false;
            }
            adjacency[i].push_back(j);
            adjacency[j].push_back(i);
        }
    }
    // Sort adjacency lists
    for (auto& vertex: adjacency) {
        sort(vertex.begin(), vertex.end());
    }
    return true;
}

bool Graph::isConnected() {
    vector visited(size(adjacency), false);
    queue<int> queue;
    queue.push(0);
    visited[0] = true;

    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();

        for (int neighbour: adjacency[current]) {
            if (!visited[neighbour]) {
                visited[neighbour] = true;
                queue.push(neighbour);
            }
        }
    }

    return find(visited.begin(), visited.end(), false) == visited.end();
}

bool Graph::isRegular() {
    const int degree = size(adjacency[0]);
    for (const auto& vertex: adjacency) {
        if (size(vertex) != degree) {
            return false;
        }
    }
    return true;
}

bool Graph::isBipartite() {
    vector<int> color(size(adjacency), -1);
    queue<int> queue;
    queue.push(0);
    color[0] = 0;

    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();

        for (int neighbour: adjacency[current]) {
            if (color[neighbour] == -1) {
                color[neighbour] = 1 - color[current];
                queue.push(neighbour);
            } else if (color[neighbour] == color[current]) {
                return false;
            }
        }
    }
    return true;
}

void Graph::setOrbits(vector<vector<int>> orbits) {
    this->orbits = orbits;

    // The orbits correspond to the line graph nodes
    // Map nodes of line graph to edges of this graph
    for (auto & orbit: orbits) {
        vector<shared_ptr<Edge>> orbitEdges;
        for (const int index: orbit) {
            orbitEdges.push_back(edges[index]);
        }

        if (orbitEdges.front()->start == orbitEdges.front()->end) {
            // orbits contain loops -> split into loops and semi edges
            vector<shared_ptr<Edge>> loops;
            vector<shared_ptr<Edge>> semiEdges;
            for (const auto& edge: orbitEdges) {
                if (edge->reverseEdge != nullptr) {
                    loops.push_back(edge);
                } else {
                    semiEdges.push_back(edge);
                }
            }
            edgeOrbits.push_back(loops);
            edgeOrbits.push_back(semiEdges);
        } else {
            edgeOrbits.push_back(orbitEdges);
        }
    }
}


int Graph::getShortestCycle(int vertex) {
    int smallest = INT_MAX;
    for (int neighbour: adjacency[vertex]) {
        int distance = 0;
        try {
            distance = this->distance(vertex, neighbour, false).first + 1;
        } catch (invalid_argument) {
            continue;
        }

        if (distance < smallest) {
            smallest = distance;
        }
    }
    return smallest;
}

vector<int> Graph::getShortestCycles() {
    vector<int> shortestCycles(size(adjacency), INT_MAX);
    for (int vertex = 0; vertex<size(adjacency); vertex++) {
        shortestCycles[vertex] = getShortestCycle(vertex);
    }
    return shortestCycles;
}

pair<int, vector<int>> Graph::distance(const int vertex1, const int vertex2, const bool allowDirect,
                                       const bool onlyTreeEdges) {
    // Do a BFS from the vertex
    vector distances(size(adjacency), -1);
    vector parent(size(adjacency), -1);
    queue<int> queue;

    distances[vertex1] = 0;
    queue.push(vertex1);

    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();

        for (int neighbour: adjacency[current]) {
            if (current == vertex1 && neighbour == vertex2 && !allowDirect) {
                continue;
            }

            if (onlyTreeEdges) {
                bool found = false;
                for (const shared_ptr<Edge>& e : neighbourToEdge[current][neighbour]) {
                    if (e->isPartOfTree) {
                        found = true;
                        break;
                    }
                }
                if (!found) continue;
            }

            if (distances[neighbour] == -1) {
                distances[neighbour] = distances[current] + 1;
                parent[neighbour] = current;
                queue.push(neighbour);
            }

            // Stop if we reached the second vertex
            if (neighbour == vertex2) {
                vector<int> path;
                for (int v = vertex2; v != -1; v = parent[v]) {
                    path.push_back(v);
                }
                return {distances[vertex2], path};
            }
        }
    }

    throw invalid_argument("The two vertices are not connected.");
}

vector<shared_ptr<Edge>> &Graph::getEdgesWithoutDefaultVoltage() {
    if (size(edgesWithoutSpanningTree) == 0) {
        createEdgesWithDefaultVoltage();
    }

    return edgesWithoutSpanningTree;
}

void Graph::createEdgesWithDefaultVoltage() {
    set<int> connected;
    set<int> explored;
    set<int> unexplored;
    set<array<int, 2>> spanningTree;

    unexplored.insert(0);

    while (size(explored) != size(adjacency)) {
        int current = *unexplored.begin();
        unexplored.erase(current);
        explored.insert(current);
        connected.insert(current);

        for (auto v: adjacency[current]) {
            if (!connected.contains(v)) {
                connected.insert(v);
                unexplored.insert(v);

                array edge = {current, v};
                spanningTree.insert(edge);

                array back = {v, current};
                spanningTree.insert(back);
            }
        }
    }

    for (int v = 0; v<size(adjacency); v++) {
        int loops = 0;
        for (auto u: adjacency[v]) {
            array edge = {v, u};
            auto e = make_shared<Edge>(v, u);

            if (v == u) {
                // add return edge for loop
                auto e2 = make_shared<Edge>(u, v);
                e->reverseEdge = e2;
                e2->reverseEdge = e;

                edges.push_back(e2);
                neighbourToEdge[u][v].push_back(e2);
                edgesWithoutSpanningTree.push_back(e2);
                loops++;
            }

            if (!spanningTree.contains(edge)) {
                edgesWithoutSpanningTree.push_back(e);
            } else {
                e->voltage = 0;
                e->isPartOfTree = true;
                // Erase in case of double edges
                spanningTree.erase(edge);
            }
            edges.push_back(e);
            neighbourToEdge[v][u].push_back(e);
        }
        for (int i = 0; i < loops; i++) {
            adjacency[v].push_back(v);
        }
    }

    // Set reverse edges
    for (const auto& e: edges) {
        if (e->reverseEdge == nullptr) {
            vector<shared_ptr<Edge>> options = neighbourToEdge[e->end][e->start];
            for (const auto& option: options) {
                if (option->reverseEdge == nullptr && option->isPartOfTree == e->isPartOfTree) {e->reverseEdge = option;
                    option->reverseEdge = e;
                    break;
                }
            }
        }
    }
}

vector<vector<int>> Graph::getLineGraph() {
    vector lineGraph(size(edges), vector<int>());

    for (int i = 0; i<size(edges); i++) {
        for (int j = i; j<size(edges); j++) {
            if (edges[i]->start == edges[j]->start || edges[i]->start == edges[j]->end ||
                edges[i]->end == edges[j]->start || edges[i]->end == edges[j]->end) {
                lineGraph[i].push_back(j);
                lineGraph[j].push_back(i);
            }
        }
    }

    return lineGraph;
}

vector<vector<int>> Graph::getCanonicalDoubleCover() {
    const int newSize= 2*size(adjacency);
    const int n = size(adjacency);
    vector doubleCover(newSize, vector<int>());

    for (int u = 0; u<n; u++) {
        for (int v: adjacency[u]) {

            if (u > v) {
                continue;
            }

            int idx1=0*n+u;
            int idx2=1*n+v;
            doubleCover[idx1].push_back(idx2);
            doubleCover[idx2].push_back(idx1);

            idx1=0*n+v;
            idx2=1*n+u;
            doubleCover[idx1].push_back(idx2);
            doubleCover[idx2].push_back(idx1);
        }
    }

    return doubleCover;
}

bool dfsHamiltonian(int node, const vector<vector<int>>& adj, vector<bool>& visited, vector<int>& path, const int remaining) {
    visited[node] = true;
    path.push_back(node);

    if (remaining == 1) return true; // Found a Hamiltonian path

    for (const int neighbor : adj[node]) {
        if (!visited[neighbor]) {
            if (dfsHamiltonian(neighbor, adj, visited, path, remaining - 1)) {
                return true;
            }
        }
    }

    visited[node] = false;
    path.pop_back();
    return false;
}

vector<int> Graph::findHamiltonianPath() {
    const int n = adjacency.size();
    vector<int> path;
    vector visited(n, false);

    for (int start = 0; start < n; ++start) {
        path.clear();
        fill(visited.begin(), visited.end(), false);
        if (dfsHamiltonian(start, adjacency, visited, path, n)) {
            return path; // Return the first Hamiltonian path found
        }
    }

    return {}; // Return empty vector if no Hamiltonian path exists
}



int Graph::main_Graph() {
    Graph graph = parseGraph6String("IsP@OkWHG");
    graph.printAdjacency();
    graph.createEdgesWithDefaultVoltage();
    printf("Girth: %d\n", graph.getGirth());
    vector<int> shortestCycles = graph.getShortestCycles();
    printf("Shortest cycles: ");
    for (int i: shortestCycles) {
        printf("%d ", i);
    }
    return 0;
}
