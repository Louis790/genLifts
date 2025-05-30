#ifndef OTHER_H
#define OTHER_H

#include <iostream>
#include <variant>

#include "CostFunctions.h"
#include "../../graph/Graph.h"
#include "../../group/Group.h"


struct canonicalStats {
    int notFiltered;
    int automorphismFiltered;
    int graphOrbitFiltered;
};

struct tabuSearchSetup {
    double timeLimit;
    int maxIterations;
    mutable int tabuSize = 1;
    double tabuSizeMult;
    int PerturbAfterNoImprovement;
    int neighbourMinGirth;
    CostFunction costFunction;
};

struct bruteForceSetup {
    double timeLimit;
};

struct mainSetup {
    int seed = 0;
    int k = 3;
    int n = 2;
    mutable int minGirth = 3;
    int logEvery = 10;
    bool useGraphAutomorphisms = true;
    bool useGroupAutomorphisms = true;
};

struct runSetup {
    mainSetup ms;
    enum AlgorithmType { BTA, TabuSearch} algorithm;
    variant<bruteForceSetup, tabuSearchSetup> algorithmSetup;
    std::ostream& out = std::cout;

    runSetup(
        const mainSetup& ms,
        AlgorithmType algo,
        const std::variant<bruteForceSetup, tabuSearchSetup>& setup
    )
        : ms(ms), algorithm(algo), algorithmSetup(setup) {}

    runSetup(
        const mainSetup& ms,
        AlgorithmType algo,
        const std::variant<bruteForceSetup, tabuSearchSetup>& setup,
        std::ostream& out
    )
        : ms(ms), algorithm(algo), algorithmSetup(setup), out(out) {}
};

/**
 * Determines at which index voltage assignment for the given list of edges is no longer canonical for the given group.
 * An assignment is canonical if:
 *  - No group automorphism maps [edge->voltage for edge in edges] to a lexicographically smaller list
 *  - Every edge orbit is monotonic
 *
 * @param edges A vector containing the edges with voltage assignments, assumed to always be given in the same order.
 * @param graphAutomorphisms A vector containing the edge automorphisms of the graph, filtered to only the edges vector above.
 * @param group A group that may contain orbits.
 * @param setup A runSetup struct containing the settings for the run. Indicates whether to use automorphisms and orbits.
 * @return An integer indicating at what index the assignment is no longer canonical. -1 if the entire assignment is canonical
 */
int getNonCanonicalIndex(
    const vector<shared_ptr<Edge>>& edges,
    Graph& graphAutomorphisms,
    const Group& group,
    const runSetup& setup
);


bool cannotAchieveMinGirth(
    const Group& group,
    Graph& graph,
    const shared_ptr<Edge>& edge,
    const runSetup& setup
    );

bool findPath(int current, int target, int depth, int maxLength, const shared_ptr<Edge>& lastTaken, int currentVoltage,
    Graph& graph, const Group& group, bool strict);
bool isKgnoGraph(const Group& group, Graph& graph, const runSetup& setup, const vector<shared_ptr<Edge>>& edges);


vector<shared_ptr<Edge>> filterInverses(const vector<shared_ptr<Edge>>& edges);
vector<vector<shared_ptr<Edge>>> filterOrbitInverses(const Graph& graph, const vector<shared_ptr<Edge>>& legalEdges);
vector<vector<shared_ptr<Edge>>> filterEdgeAutomorphisms(const Graph& graph, const vector<shared_ptr<Edge>>& legalEdges);
vector<int> getVoltagesForEdge(Graph& graph, const Group& group, const shared_ptr<Edge>& edge, int minGirth);

canonicalStats getCanonicalStats();

bool initialAssignment(
    const vector<shared_ptr<Edge>>& edges,
    const Group& group,
    Graph& graph,
    map<shared_ptr<Edge>, vector<int>>& legalEdgeVoltages,
    const runSetup& setup
);

#endif //OTHER_H
