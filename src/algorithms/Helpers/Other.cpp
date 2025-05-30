#include "Other.h"

#include <chrono>

int notFiltered = 0;
int automorphismFiltered = 0;
int graphOrbitFiltered = 0;

int getNonCanonicalIndex(
    const vector<shared_ptr<Edge>>& edges,
    Graph& graph,
    const Group& group,
    const runSetup& setup
) {

    if (setup.ms.useGroupAutomorphisms) {
        for (auto& automorphism : group.automorphisms) {
            for (int i = 0; i < size(edges) ; i++) {
                auto& edge = edges[i];
                if (automorphism[edge->voltage] > edge->voltage) {
                    break;
                } if(automorphism[edge->voltage] < edge->voltage) {
                    automorphismFiltered++;
                    return i;
                }
            }
        }
    }

    if (setup.ms.useGraphAutomorphisms) {
        for (auto& automorphism : graph.edgeAutomorphisms) {
            for (int i = 0; i < size(edges) ; i++) {
                auto& edge = edges[i];
                if (automorphism[i]->voltage > edge->voltage) {
                    break;
                } if(automorphism[i]->voltage < edge->voltage) {
                    automorphismFiltered++;
                    int otherEdgeIndex = -1;
                    for (int j = 0; j < size(edges); j++) {
                        if (edges[j] == automorphism[i]) {
                            otherEdgeIndex = j;
                            break;
                        }
                    }
                    return max(i, otherEdgeIndex);
                }
            }
        }
    }

    notFiltered++;
    return -1;
}

bool findPath(int current, int target, int depth, int maxLength, const shared_ptr<Edge>& lastTaken, int currentVoltage, Graph& graph, const Group& group, bool
              strict) {
    if (depth >= maxLength) return false;

    // This skips paths consisting of a single loop or semi edge, but those should be caught by the static checks
    for (const auto& [neighbor, edgesToNeighbor] : graph.neighbourToEdge[current]) {
        for (const auto& edge : edgesToNeighbor) {
            // Skip if voltage not set
            if (edge->voltage == -1) continue;

            // Dont take the reverse edge of the previously taken edge
            if (edge->reverseEdge == lastTaken) continue;

            // Semi edges are their own inverse
            if (edge->start == edge->end && edge->reverseEdge == nullptr && edge == lastTaken) continue;

            currentVoltage = group.multiplicationTable[currentVoltage][edge->voltage];
            if (neighbor == target && currentVoltage == 0 && (!strict || depth == maxLength - 1)) {
                return true;
            }

            bool foundPath = findPath(neighbor, target, depth + 1, maxLength, edge, currentVoltage, graph, group, strict);
            if (foundPath) {
                return true;
            }

            currentVoltage = group.multiplicationTable[currentVoltage][group.inverse[edge->voltage]];
        }
    }
    return false;
}

bool cannotAchieveMinGirth(const Group& group, Graph& graph, const shared_ptr<Edge>& edge,
    const runSetup& setup) {
    shared_ptr<Edge> takenEdge = edge->reverseEdge;
    if (takenEdge == nullptr) {
        takenEdge = edge;
    }
    return findPath(edge->start, edge->end, 1, setup.ms.minGirth-1, takenEdge, takenEdge->voltage, graph, group, false);
}

bool isKgnoGraph(const Group& group, Graph& graph, const runSetup& setup, const vector<shared_ptr<Edge>>& edges) {
    // assume < minGirth was already checked
    int girth = setup.ms.minGirth;

    // Check if girth is exactly minGirth
    bool found = false;
    for (const shared_ptr<Edge>& edge : edges) {
        shared_ptr<Edge> takenEdge = edge->reverseEdge;
        if (takenEdge == nullptr) {
            takenEdge = edge;
        }
        if (findPath(edge->start, edge->end, 1, girth, takenEdge, takenEdge->voltage, graph, group, true)) {
            found = true;
            break;
        }
    }
    if (!found) {
        return false;
    }

    // check if no girth+1 cycles exist
    for (const shared_ptr<Edge>& edge : edges) {
        shared_ptr<Edge> takenEdge = edge->reverseEdge;
        if (takenEdge == nullptr) {
            takenEdge = edge;
        }
        if (findPath(edge->start, edge->end, 1, girth+1, takenEdge, takenEdge->voltage, graph, group, true)) {
            return false;
        }
    }
    return true;
}

vector<shared_ptr<Edge>> filterInverses(const vector<shared_ptr<Edge>>& edges) {
    vector<shared_ptr<Edge>> result;
    set<shared_ptr<Edge>> illegal;

    for (const shared_ptr<Edge>& edge : edges) {
        if (illegal.contains(edge)) {
            continue;
        }

        result.push_back(edge);
        illegal.insert(edge->reverseEdge);
    }

    return result;
}

vector<vector<shared_ptr<Edge>>> filterOrbitInverses(const Graph& graph, const vector<shared_ptr<Edge>>& legalEdges) {
    vector<vector<shared_ptr<Edge>>> result;
    for (const vector<shared_ptr<Edge>>& orbit : graph.edgeOrbits) {
        vector<shared_ptr<Edge>> filteredOrbit;
        for (const shared_ptr<Edge>& edge : orbit) {
            if (find(legalEdges.begin(), legalEdges.end(), edge) != legalEdges.end()) {
                filteredOrbit.push_back(edge);
            }
        }
        // skip if orbit is empty
        if (size(filteredOrbit) == 0) {
            continue;
        }
        result.push_back(filteredOrbit);
    }

    return result;
}

vector<int> getVoltagesForEdge(Graph& graph, const Group& group, const shared_ptr<Edge>& edge, int minGirth) {
    // semi edge
    if (edge->start == edge->end && edge->reverseEdge == nullptr) {
        vector<int> result;
        for (int i = 0; i < size(group.multiplicationTable); i++) {
            if (group.powerToIdentity[i] == 2) {
                result.push_back(i);
            }
        }
        return result;
    }

    int loopSize;
    if (edge-> start == edge->end) {
        loopSize = 1;
    } else {
        loopSize = graph.distance(edge->start, edge->end, true, true).first+1;
    }

    vector<int> result;
    for (int i = 0; i < size(group.multiplicationTable); i++) {
        int loopsNeeded = group.powerToIdentity[i];
        if (loopsNeeded == INT_MAX) {
            result.push_back(i);
            continue;
        }

        if (loopSize * loopsNeeded >= minGirth) {
            result.push_back(i);
        }
    }

    return result;
}

canonicalStats getCanonicalStats() {
    return {notFiltered, automorphismFiltered, graphOrbitFiltered};
}

bool initialAssignment(
    const vector<shared_ptr<Edge>>& edges,
    const Group& group,
    Graph& graph,
    map<shared_ptr<Edge>, vector<int>>& legalEdgeVoltages,
    const runSetup& setup
)
{
    vector<int> stack;
    stack.reserve(size(edges));
    for (int i = 0; i < size(edges); i++) {
        stack.push_back(-1);
    }

    for(const shared_ptr<Edge>& edge : edges) {
        edge->setVoltage(-1, -1);
    }

    auto start = chrono::high_resolution_clock::now();

    int currentIndex = 0;
    while (true) {
        // Backtrack if needed
        if (stack[currentIndex] == size(legalEdgeVoltages[edges[currentIndex]]) - 1) {
            stack[currentIndex] = -1;
            edges[currentIndex]->setVoltage(-1, -1);
            currentIndex--;

            if (currentIndex < 0) {
                break;
            }
            continue;
        }

        // Give up after 2 seconds
        if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - start).count() > 2) {
            clog << "Initial assignment took too long, giving up" << endl;
            return false;
        }

        // Try next voltage
        stack[currentIndex]++;
        edges[currentIndex]->setVoltage(legalEdgeVoltages[edges[currentIndex]][stack[currentIndex]],
                                        group.inverse[legalEdgeVoltages[edges[currentIndex]][stack[currentIndex]]]);

        if (cannotAchieveMinGirth(group, graph, edges[currentIndex], setup)) {
            continue;
        }

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
            return true;
        } else {
            currentIndex++;
        }
    }
    clog << "Could not find an initial assignment" << endl;
    return false;
}
