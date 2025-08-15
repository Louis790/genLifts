#include "Other.h"

#include <chrono>

int notFiltered = 0;
int automorphismFiltered = 0;
int graphOrbitFiltered = 0;

// Returns the largest index to backtrack to if min_tau t(b) < a lexicographically; else -1.
// positionMap[i] = index in 'edges' that b[i] came from.
// For pure group comparison (b == a), positionMap is the identity: positionMap[i] = i and returnOnFirstViolation is true.
int compareLexUnderGroup(
    const vector<int>& a,
    const vector<int>& b,
    const Group& group,
    const vector<int>& positionMap,
    bool returnOnFirstViolation
) {
    const int n = static_cast<int>(a.size());
    if (static_cast<int>(b.size()) != n || static_cast<int>(positionMap.size()) != n) return -1;

    int idx = -1;
    for (const auto& tau : group.automorphisms) {
        for (int i = 0; i < n; ++i) {
            int tb = tau[b[i]];
            if (tb < a[i]) {
                // backtrack to the later-assigned of the pair
                idx = max(idx, max(i, positionMap[i]));
            }
            if (tb > a[i]) break;  // this t can't be lexicographically smaller than a
            if (returnOnFirstViolation && idx != -1) return idx;

        }
    }
    return -1;
}

int getNonCanonicalIndex(
    const vector<shared_ptr<Edge>>& edges,
    vector<vector<shared_ptr<Edge>>>& graphAutomorphisms,
    const Group& group,
    const runSetup& setup
) {
    const int n = static_cast<int>(edges.size());

    // Current voltage vector v in 'edges' order
    vector<int> v(n);
    for (int i = 0; i < n; ++i) v[i] = edges[i]->voltage;

    unordered_map<const Edge*, int> edgeIndex;
    edgeIndex.reserve(n * 2);
    for (int i = 0; i < n; ++i) edgeIndex[edges[i].get()] = i;

    // Pure group automorphisms
    if (setup.ms.useGroupAutomorphisms) {
        vector<int> idMap(n);
        for (int i = 0; i < n; ++i) idMap[i] = i;

        int idx = compareLexUnderGroup(v, v, group, idMap, true);
        if (idx != -1) return idx;
    }

    // Compare min_t t(s(v)) vs v
    if (setup.ms.useGraphAutomorphisms) {
        int idx = -1;
        for (const auto& sigmaEdges : graphAutomorphisms) {
            // s(v) in edges order
            vector<int> sv(n);
            // Maps position i in sv to original index of sigmaEdges[i]
            vector<int> sigmaIndexMap(n);
            for (int i = 0; i < n; ++i) {
                const Edge* e = sigmaEdges[i].get();
                int j = edgeIndex.find(e)->second;
                sv[i] = edges[j]->voltage;
                sigmaIndexMap[i] = j;
            }

            idx = max(compareLexUnderGroup(v, sv, group, sigmaIndexMap, false), idx);
        }
        if (idx != -1) return idx;
    }

    return -1; // canonical
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
            int breaksCanon = getNonCanonicalIndex(edges, graph.edgeAutomorphisms, group, setup);
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
