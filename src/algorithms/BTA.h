#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

#include "../graph/Graph.h"
#include "../group/Group.h"
#include "Helpers/Other.h"

bool filteredBTA(Graph& graph, const Group& group, const runSetup& setup);

int expectedNumberOfLifts(Graph& graph, const Group& group, const runSetup& setup);

int main_BruteForce();


#endif //BRUTEFORCE_H
