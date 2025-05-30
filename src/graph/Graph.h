#ifndef GRAPH_H
#define GRAPH_H

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Edge.h"

using std::vector;
using std::string;
using std::map;
using std::shared_ptr;
using std::pair;

#define MAXVERTICES 4096
// Unsafe because no defined behaviour if character = 0. ctz and clz work with 32 bit numbers.
#define unsafePrev(character, current) (__builtin_ctz(character) - current >= 0 ? -1 : current -__builtin_clz((character) << (32 - current)) - 1)
#define custom_prev(character,current) (character ? unsafePrev(character,current) : -1)

class Graph {
public:
    Graph();

    vector<int> findHamiltonianPath();
    static int main_Graph();

    /**
     * Get the adjacency list of vertex i of this graph.
     * @param i the vertex.
     * @return the adjacency list of the vertex.
     */
    vector<int> &operator[](int i);
    void printAdjacency();
    string getGraph6(bool makeCanonical);

    void saveAsMatrix(const string& path = R"(../../Data/SavedGraph.mat)");


    /**
     * Get the girth of the graph.
     * @return the girth of the graph.
     */
    int getGirth();

    /**
     * Get whether the graph is one-edge-connected. I.e. whether there is a path between every pair of vertices if one edge is removed.
     * @return A boolean indicating whether the graph is one-edge-connected.
     */
    bool isOneEdgeConnected();

    /**
     * Get whether the graph is connected. I.e. whether there is a path between every pair of vertices.
     * @return A boolean indicating whether the graph is connected.
     */
    bool isConnected();

    /**
     * Get whether the graph is regular.
     * @return A boolean indicating whether the graph is regular.
     */
    bool isRegular();

    /**
     * Get whether the graph is bipartite.
     * @return A boolean indicating whether the graph is bipartite.
     */
    bool isBipartite();

    /**
     * Set the orbits of this graph to the given orbits.
     * Also populates edgeOrbits.
     */
    void setOrbits(vector<vector<int>> orbits);

    /**
     * Get the length of the shortest cycle containing the vertex.
     * @param vertex
     * @return the length of the shortest cycle containing the vertex.
     */
    int getShortestCycle(int vertex);

    /**
     * Get the distance between two vertices.
     * @param vertex1
     * @param vertex2
     * @param allowDirect whether to allow direct edges between the two vertices.
     * @param onlyTreeEdges whether only edges that are part of the tree should be used.
     * @return the distance between the two vertices.
     */
    pair<int, vector<int>> distance(int vertex1, int vertex2, bool allowDirect, bool onlyTreeEdges = false);

    vector<shared_ptr<Edge>> &getEdgesWithoutDefaultVoltage();

    /**
     * Get the lengths of the shortest cycles in the graph.
     * @return A vector containing the lengths of the shortest cycles in the graph for each vertex.
     */
    vector<int> getShortestCycles();

    /**
     * Constructs a random spanning tree of the graph and returns the edges that are not in the spanning tree.
     * @return A vector containing the edges that are not in the spanning tree.
     */
    void createEdgesWithDefaultVoltage();

    vector<vector<int>> adjacency;
    map<int, map<int, vector<shared_ptr<Edge>>>> neighbourToEdge;
    vector<shared_ptr<Edge>> edges;
    vector<shared_ptr<Edge>> edgesWithoutSpanningTree;
    shared_ptr<Graph> baseGraph = nullptr;
    int maxSemiEdgesPerVertex = 0;

    /**
     * The orbits the line graph of this graph.
     * {{0 1} {2}} means that vertices 0 and 1 are in the same orbit and vertex 2 is in its own orbit.
     */
    vector<vector<int>> orbits;

    /**
     * The edge automorphisms of this graph.
     * {{0 2 1}} means that swapping edges 1 and 2 and keeping edge 0 is an edge automorphism.
     * Order is the same as the edges vector.
     */
    vector<vector<shared_ptr<Edge>>> edgeAutomorphisms;

    /**
     * Get the adjacency matrix of the line graph of this graph.
     * @return the adjacency matrix of the line graph of this graph.
     */
    vector<vector<int>> getLineGraph();

    /**
     * Get the adjacency matrix of the double cover of this graph.
     * @return the adjacency matrix of the double cover of this graph.
     */
    vector<vector<int>> getCanonicalDoubleCover();

    /**
     * The edge orbits of this graph.
     * Two edges (u, v) and (x, y) are in the same orbit if u and x are in the same orbit and v and y are in the same orbit.
     * {{0 1} {2}} means that edges 0 and 1 are in the same orbit and edge 2 is in its own orbit.
     * Uses indexes of the edges in the edges vector.
     */
    vector<vector<shared_ptr<Edge>>> edgeOrbits;
};

#endif //GRAPH_H