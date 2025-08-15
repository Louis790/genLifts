#ifndef VERTEXEDGEGIRTHREGULAR_H
#define VERTEXEDGEGIRTHREGULAR_H

#include "../../graph/Graph.h"
#include "../../group/Group.h"

pair<int, int> getBothGirthRegularLambda(Graph& graph, int girth);
int getEdgeGirthRegularLambda(Graph& graph, int girth, bool useCache = false);
int getVertexGirthRegularLambda(Graph& graph, int girth, bool useCache = false);
bool hasNoCycleOfLength(const Graph& graph, int length);

#endif //VERTEXEDGEGIRTHREGULAR_H
