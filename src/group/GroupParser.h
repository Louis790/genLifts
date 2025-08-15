#ifndef GROUPPARSER_H
#define GROUPPARSER_H
#include <string>

#include "Group.h"
#include "../graph/Graph.h"

using namespace std;


Group fromPath(string path);
Group fromString(string s);
vector<vector<vector<int>>> parseCycleGenerators(const vector<string>& generators);
vector<vector<shared_ptr<Edge>>>  getAutomorphisms(const vector<string>& generators, Graph& graph, int maxAutomorphisms);
Group cyclicGroup(int n);

#endif //GROUPPARSER_H
