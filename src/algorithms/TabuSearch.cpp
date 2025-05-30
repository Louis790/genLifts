#include "TabuSearch.h"

#include <chrono>
#include <deque>
#include <iostream>

#include "Filter.h"
#include "Helpers/Other.h"
#include "Helpers/RNG.h"

struct change {
    shared_ptr<Edge> edge;
    int voltage;

    bool operator==(const change& other) const {
        return edge == other.edge && voltage == other.voltage;
    }

    bool operator<(const change& other) const {
        if (edge == other.edge) {
            return voltage < other.voltage;
        }
        return edge < other.edge;
    }
};

bool perturb(
    const vector<shared_ptr<Edge>>& edges,
    const Group& group,
    Graph& graph,
    map<shared_ptr<Edge>, vector<int>>& legalEdgeVoltages,
    const runSetup& setup
)
{
    auto shuffledVoltages = legalEdgeVoltages;
    for (auto& [edge, voltages]: shuffledVoltages) {
        shuffleVec(voltages);
    }
    return initialAssignment(edges, group, graph, shuffledVoltages, setup);
}

void doTabuSearch(
    const vector<shared_ptr<Edge>>& edges,
    Graph& graph,
    const Group& group,
    map<shared_ptr<Edge>, vector<int>>& legalEdgeVoltages,
    const runSetup& setup
)
{
    tabuSearchSetup tsSetup = get<tabuSearchSetup>(setup.algorithmSetup);
    int originalMinGirth = setup.ms.minGirth;
    int neighbourMinGirth = max(tsSetup.neighbourMinGirth, 3);
    // clog << "Using minGirth: " << neighbourMinGirth << " for neighbour search." << endl;

    auto start = chrono::high_resolution_clock::now();

    setup.ms.minGirth = neighbourMinGirth;
    auto foundInitial = initialAssignment(edges, group, graph, legalEdgeVoltages, setup);
    setup.ms.minGirth = originalMinGirth;

    if (!foundInitial) {
        return;
    }

    int iterations = 0;
    int iterationsSinceLastImprovement = 0;
    double globalBestScore = -1;
    deque<change> tabuList;
    set<change> tabuSet;
    shared_ptr<Edge> chosenEdge = edges[0];
    while (true) {
        iterations++;
        iterationsSinceLastImprovement++;
        if (iterations >= tsSetup.maxIterations) {
            break;
        }

        if (iterationsSinceLastImprovement >= tsSetup.PerturbAfterNoImprovement) {
            // clog << "No improvement for " << tsSetup.PerturbAfterNoImprovement << " iterations, perturbing." << endl;
            auto perturbSuccess = perturb(edges, group, graph, legalEdgeVoltages, setup);
            if (!perturbSuccess) {
                // clog << "Perturbation failed, giving up." << endl;
                return;
            }
            tabuSet.clear();
            tabuList.clear();
        }

        if (tsSetup.timeLimit > 0) {
            auto end = chrono::high_resolution_clock::now();
            if (chrono::duration<double> elapsed = end - start; elapsed.count() > tsSetup.timeLimit) {
                break;
            }
        }

        if (!cannotAchieveMinGirth(group, graph, chosenEdge, setup)) {
            Graph lifted = lift(graph, group);
            int girth = lifted.getGirth();
            if (girth >= setup.ms.minGirth) {
                filterAndWrite(lifted, girth, setup.out);
            }
        }

        vector<change> bestChanges;
        int options = 0;
        double bestScore = -INT_MAX;
        for (auto& edge: edges) {
            for (int voltage: legalEdgeVoltages[edge]) {
                int oldVoltage = edge->voltage;
                if (oldVoltage == voltage) {
                    continue;
                }
                edge->setVoltage(voltage, group.inverse[voltage]);

                setup.ms.minGirth = neighbourMinGirth;
                if (cannotAchieveMinGirth(group, graph, edge, setup)) {
                    setup.ms.minGirth = originalMinGirth;
                    edge->setVoltage(oldVoltage, group.inverse[oldVoltage]);
                    continue;
                }
                setup.ms.minGirth = originalMinGirth;

                int breaksCanon = getNonCanonicalIndex(edges, graph, group, setup);
                if (breaksCanon != -1) {
                    edge->setVoltage(oldVoltage, group.inverse[oldVoltage]);
                    continue;
                }

                change c = {edge, voltage};
                double score;
                options++;
                score = tsSetup.costFunction(graph, group, edge);

                bool isTabu = tabuSet.contains(c);
                bool isGlobalBest = score > globalBestScore && perturb == 0;
                bool isBest = score > bestScore;
                bool isClose = abs(score - bestScore) < 1e-6;
                // clog << "Score: " << score << ", Best: " << bestScore << ", IsBest: " << isBest
                //      << ", IsTabu: " << isTabu << ", IsClose: " << isClose << endl;

                if ((isBest && !isTabu) || isGlobalBest) {
                    bestScore = score;
                    bestChanges.clear();
                    bestChanges.push_back(c);
                    if (isGlobalBest) {
                        iterationsSinceLastImprovement = 0;
                        globalBestScore = score;
                        // clog << "New global best score: " << score << endl;
                    }
                }
                else if (isClose && !isTabu) {
                    bestChanges.push_back(c);
                }
                edge->setVoltage(oldVoltage, group.inverse[oldVoltage]);
            }
        }

        // clog << "Iteration: " << iterations << ", Options: " << options << ", Best score: " << bestScore << endl;

        if (size(bestChanges) == 0) {
            // clog << "No options left, perturbing." << endl;
            auto perturbSuccess = perturb(edges, group, graph, legalEdgeVoltages, setup);
            if (!perturbSuccess) {
                // clog << "Perturbation failed, giving up." << endl;
                return;
            }
            iterationsSinceLastImprovement = 0;
            tabuSet.clear();
            tabuList.clear();
            continue;
        }
        auto bestChange = bestChanges[getRandomInt(0, size(bestChanges)-1)];

        change undoChange = {bestChange.edge, group.inverse[bestChange.edge->voltage]};
        tabuList.push_back(undoChange);
        tabuSet.insert(undoChange);
        if (size(tabuList) > tsSetup.tabuSize) {
            tabuSet.erase(tabuList.front());
            tabuList.pop_front();
        }
        bestChange.edge->setVoltage(bestChange.voltage, group.inverse[bestChange.voltage]);
        chosenEdge = bestChange.edge;
    }
}

void tabuSearch(Graph& graph, const Group& group, const runSetup& setup) {
    vector<shared_ptr<Edge>> edges = graph.getEdgesWithoutDefaultVoltage();
    vector<shared_ptr<Edge>> edgesWithoutInverse = filterInverses(edges);
    // const vector<vector<shared_ptr<Edge>>> filteredEdgeOrbits = filterOrbitInverses(graph, edgesWithoutInverse);
    // const auto graphAutomorphisms = filterEdgeAutomorphisms(graph, edgesWithoutInverse);

    map<shared_ptr<Edge>, vector<int>> legalEdgeVoltages;
    for (const shared_ptr<Edge>& edge : edgesWithoutInverse) {
        legalEdgeVoltages[edge] = getVoltagesForEdge(graph, group, edge, setup.ms.minGirth);

        // No good edges
        if (size(legalEdgeVoltages[edge]) == 0) {
            return;
        }
    }

    tabuSearchSetup tsSetup = get<tabuSearchSetup>(setup.algorithmSetup);
    tsSetup.tabuSize = group.multiplicationTable.size() * tsSetup.tabuSizeMult;

    return doTabuSearch(edgesWithoutInverse, graph, group, legalEdgeVoltages, setup);
}
