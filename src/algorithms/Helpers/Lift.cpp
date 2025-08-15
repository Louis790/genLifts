#include "Lift.h"

#include <iostream>

#include "../../group/GroupParser.h"

Graph lift(Graph graph, const Group& group) {
    int vertexCount = size(graph.adjacency) * size(group.multiplicationTable);
    Graph lifted;
    for (int i = 0; i < vertexCount; i++) {
        lifted.adjacency.emplace_back();
    }

    const int multiplicationSize = size(group.multiplicationTable);

    // [v1g1, v1g2, ..., v2g1, ...]
    for (int v = 0; v < vertexCount; v++) {
        for (int u = v; u < vertexCount; u++) {
            int vGroup = v % multiplicationSize;
            int uGroup = u % multiplicationSize;

            int vVert = (v - vGroup) / multiplicationSize;
            int uVert = (u - uGroup) / multiplicationSize;

            bool inserted = false;

            // (v, u) in graph
            if (!graph.neighbourToEdge[vVert][uVert].empty()) {
                // uGroup == vGroup * voltage(edge)
                for (const shared_ptr<Edge>& e : graph.neighbourToEdge[vVert][uVert]) {
                    if (uGroup == group.multiplicationTable[vGroup][e->voltage]) {
                        lifted.adjacency[v].push_back(u);
                        lifted.adjacency[u].push_back(v);
                        // inserted = true;
                        // break;
                    }
                }
            }
            // if (inserted) {
            //     continue;
            // }

            // // (u, v) in graph
            // if (!graph.neighbourToEdge[uVert][vVert].empty()) {
            //     // vGroup == uGroup * voltage(edge)
            //     for (const shared_ptr<Edge>& e : graph.neighbourToEdge[uVert][vVert]) {
            //         if (vGroup == group.multiplicationTable[uGroup][e->voltage]) {
            //             lifted.adjacency[v].push_back(u);
            //             lifted.adjacency[u].push_back(v);
            //             break;
            //         }
            //     }
            // }
        }
    }

    // lifted.baseGraph = make_shared<Graph>(graph);

    return lifted;
}

int girthOfLift(const Graph& graph, const Group& group) {
    return lift(graph, group).getGirth();
}

int main_Indexing() {
    vector<int> vec;
    int g = 6;
    int vert = 1;
    int n = g * vert;
    for (int i = 0; i<n; i++) {
        vec.push_back(i);
    }
    // [0: g0v0, 1: g1v0, 2:g0v1, 3: g1v1, ...]
    for (int v = 0; v < size(vec); v++) {
        int vGroup = v % g;
        int vVert = (v - vGroup) / g;

        std::clog << "elem " << v << ", group " << vGroup << ", vert " << vVert << '\n';
    }
    return 0;
}