#include "BTA.h"

#include <chrono>
#include <iostream>
#include <set>

#include "Filter.h"
#include "../group/GroupParser.h"
#include "./Helpers/Lift.h"
#include "./Helpers/Other.h"

bool iterBTA(
    const vector<shared_ptr<Edge>>& edges,
    Graph& graph,
    const Group& group,
    map<shared_ptr<Edge>, vector<int>>& legalEdgeVoltages,
    const runSetup& setup
) {
    auto BFSetup = get<bruteForceSetup>(setup.algorithmSetup);
    vector<int> stack;
    stack.reserve(size(edges));
    for (int i = 0; i < size(edges); i++) {
        stack.push_back(-1);
    }

    auto start = chrono::high_resolution_clock::now();

    // Make sure the settable edges are set to -1
    for (const shared_ptr<Edge>& edge : edges) {
        edge->setVoltage(-1, -1);
    }

    int currentIndex = 0;
    while (true) {
        auto seconds = chrono::duration<double>(chrono::high_resolution_clock::now() - start).count();
        if (seconds > BFSetup.timeLimit) {
            // clog << "Time limit reached, Exiting BF" << endl;
            return false;
        }

        // Backtrack if needed
        if (stack[currentIndex] == size(legalEdgeVoltages[edges[currentIndex]]) - 1) {
            stack[currentIndex] = -1;
            edges[currentIndex]->setVoltage(-1, -1);
            currentIndex--;

            if (currentIndex < 0) {
                return true;
            }
            continue;
        }

        // Try next voltage
        stack[currentIndex]++;
        edges[currentIndex]->setVoltage(legalEdgeVoltages[edges[currentIndex]][stack[currentIndex]],
                                        group.inverse[legalEdgeVoltages[edges[currentIndex]][stack[currentIndex]]]);

        // Check if girth is achieved
        if (cannotAchieveMinGirth(group, graph, edges[currentIndex], setup)) {
            continue;
        }

        // Check lift if all voltages are set
        if (currentIndex == size(edges) - 1) {

            // Must be canonical
            int breaksCanon = getNonCanonicalIndex(edges, graph, group, setup);
            if (breaksCanon != -1) {
                // Backtrack to breaksCanon
                for(int i = currentIndex; i > breaksCanon; i--) {
                    edges[i]->setVoltage(-1, -1);
                    stack[i] = -1;
                }
                currentIndex = breaksCanon;
                continue;
            }

            // // Check for kgnog+1
            // if (!isKgnoGraph(group, graph, setup, edges)) {
            //     continue;
            // }

            Graph lifted = lift(graph, group);
            const int girth = lifted.getGirth();

            // for (int i = 0; i < size(edges); i++) {
            //     clog << "arc from " << edges[i]->start << " to " << edges[i]->end << " with voltage "
            //          << edges[i]->voltage << endl;
            // }

            if (girth < setup.ms.minGirth) {
                cerr << "Girth was smaller than expected: " << girth << endl;
                continue;
            }

            filterAndWrite(lifted, girth, setup.out);
        } else {
            currentIndex++;
        }
    }
}

bool filteredBTA(Graph& graph, const Group& group, const runSetup& setup) {
    vector<shared_ptr<Edge>> edges = graph.getEdgesWithoutDefaultVoltage();
    vector<shared_ptr<Edge>> edgesWithoutInverse = filterInverses(edges);
    // const vector<vector<shared_ptr<Edge>>> filteredEdgeOrbits = filterOrbitInverses(graph, edgesWithoutInverse);
    // const vector<vector<shared_ptr<Edge>>> graphAutomorphisms = filterEdgeAutomorphisms(graph, edgesWithoutInverse);

    map<shared_ptr<Edge>, vector<int>> legalEdgeVoltages;
    for (const shared_ptr<Edge>& edge : edgesWithoutInverse) {
        legalEdgeVoltages[edge] = getVoltagesForEdge(graph, group, edge, setup.ms.minGirth);

        // No good edges
        if (size(legalEdgeVoltages[edge]) == 0) {
            return true;
        }
    }

    if(group.involutionCount-1 < graph.maxSemiEdgesPerVertex) {
        // Not all semi edges can get a voltage
        // clog << "Not all semi edges can get a voltage" << endl;
        return true;
    }

    return iterBTA(edgesWithoutInverse, graph, group, legalEdgeVoltages, setup);
}
