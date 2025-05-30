//
// Created by veetoo on 14/03/25.
//

#include "Filter.h"

#include <queue>

#include "Helpers/VertexEdgeGirthRegular.h"

set<string> writtenGraphs = {};

void filterAndWrite(Graph& graph, int girth, std::ostream& out) {
    if (!graph.isRegular() || !graph.isConnected()) {
        return;
    }

    const int n = size(graph.adjacency);
    const int k = size(graph.adjacency[0]);
    const string g6 = graph.getGraph6(true);

    if (writtenGraphs.contains(g6)) {
        return;
    }
    writtenGraphs.insert(g6);

    const auto [edge_lambda, vertex_lambda] = getBothGirthRegularLambda(graph, girth);

    out << "(k,g)-graph - " << k << " " << girth << " " << n << " - " << g6 << endl;

    if (edge_lambda != -1) {
        out << "egr-graph - " << k << " " << girth << " " << edge_lambda << " " << n << " - " << g6 << endl;
    }

    if (vertex_lambda != -1) {
        out << "vgr-graph - " << k << " " << girth << " " << vertex_lambda << " " << n << " - " << g6 << endl;
    }

    if (girth % 2 == 1 && hasNoCycleOfLength(graph, girth + 1)) {
        out << "(k,g,g+1)-graph - " << k << " " << girth << " " << n << " - " << g6 << endl;
    }
}