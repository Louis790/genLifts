#ifndef LIFT_H
#define LIFT_H

#include "../../group/Group.h"
#include "../../graph/Graph.h"


Graph lift(Graph graph, const Group& group);
int girthOfLift(const Graph& graph, const Group& group);

#endif //LIFT_H
