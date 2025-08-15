#ifndef TABUSEARCH_H
#define TABUSEARCH_H

#include "../graph/Graph.h"
#include "../graph/Edge.h"
#include "../group/Group.h"
#include "./Helpers/Lift.h"
#include "./Helpers/CostFunctions.h"

void tabuSearch(
    Graph& graph,
    const Group& group,
    const struct runSetup& setup
    );


#endif //TABUSEARCH_H
