#include "GroupParser.h"

#include <chrono>

#include "../graph/Graph.h"
#include "../algorithms/Helpers/Other.h"

#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <regex>

Group fromPath(string path) {
    ifstream file(path);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << path << endl;
        abort();
    }

    string input;
    string line;
    // Merge all lines into one string
    while (getline(file, line)) {
        input += line;
    }

    return fromString(input);
}

Group fromString(string s) {
    stringstream ss(s);
    vector<int> currentRow;
    int num;
    char ch;
    vector<vector<int>> multiplicationTable;
    vector<vector<int>> orbits;
    vector<vector<int>> automorphisms;

    // Parse order and number
    int order = 0;
    int number = 0;
    size_t orderPos = s.find("Order:");
    size_t numberPos = s.find("Number:");

    if (orderPos != string::npos && numberPos != string::npos) {
        order = std::stoi(s.substr(orderPos + 6, numberPos - orderPos - 6));
        number = std::stoi(s.substr(numberPos + 7));
    }


    // Parse the multiplication table
    while (ss >> ch) {
        if (ch == '[') {
            currentRow.clear();
        } else if (ch == ']') {
            if (!currentRow.empty()) {
                multiplicationTable.push_back(currentRow);
                currentRow.clear();
            }
        } else if (isdigit(ch) || ch == '-') {
            ss.putback(ch);
            ss >> num;
            currentRow.push_back(num-1);
        } else if (ch == 'A') {
            // Automorphisms start
            break;
        }
    }

    // Parse the automorphisms
    while (ss >> ch) {
        if (ch == '[') {
            currentRow.clear();
        } else if (ch == ']') {
            if (!currentRow.empty()) {
                automorphisms.push_back(currentRow);
                currentRow.clear();
            }
        } else if (isdigit(ch) || ch == '-') {
            ss.putback(ch);
            ss >> num;
            currentRow.push_back(num-1);
        } else if (ch == 'O') {
            // Orbits start
            break;
        }
    }

    while (ss >> ch) {
        if (ch == '[') {
            currentRow.clear();
        } else if (ch == ']') {
            if (!currentRow.empty()) {
                orbits.push_back(currentRow);
                currentRow.clear();
            }
        } else if (isdigit(ch) || ch == '-') {
            ss.putback(ch);
            ss >> num;
            currentRow.push_back(num-1);
        }
    }

    auto g = Group(multiplicationTable, orbits, automorphisms);
    g.name = "SmallGroup(" + to_string(order) + ", " + to_string(number) + ")";
    return g;
}

vector<vector<vector<int>>> parseCycleGenerators(const vector<string>& generators) {
    vector<vector<vector<int>>> parsed;

    regex cycleRegex(R"(\(([^\)]+)\))"); // matches contents inside parentheses

    for (const string& gen : generators) {
        vector<vector<int>> cycles;
        auto begin = sregex_iterator(gen.begin(), gen.end(), cycleRegex);
        auto end = sregex_iterator();

        for (auto it = begin; it != end; ++it) {
            string cycleStr = it->str(1); // get inside of ( )
            istringstream iss(cycleStr);
            vector<int> cycle;
            int num;
            while (iss >> num) {
                cycle.push_back(num);
            }
            if (!cycle.empty())
                cycles.push_back(cycle);
        }

        parsed.push_back(cycles);
    }

    return parsed;
}

vector<int> cyclesToPermutation(const vector<vector<int>>& cycles, int n) {
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0); // identity by default

    for (const auto& cycle : cycles) {
        for (size_t i = 0; i < cycle.size(); ++i) {
            perm[cycle[i]] = cycle[(i + 1) % cycle.size()];
        }
    }

    return perm;
}

vector<int> applyPermutation(const vector<int>& f, const vector<int>& g) {
    const int n = f.size();
    vector<int> result(n);
    for (int i = 0; i < n; ++i)
        result[i] = f[g[i]];
    return result;
}

vector<vector<shared_ptr<Edge>>> getAutomorphisms(const vector<string>& generators, Graph& graph, int maxAutomorphisms) {
    // Parse each generator (a b)(c d) into a vector of vectors
    vector<vector<vector<int>>> parsedGenerators = parseCycleGenerators(generators);

    int graphSize = graph.edges.size();
    vector<shared_ptr<Edge>> edges = graph.getEdgesWithoutDefaultVoltage();
    vector<shared_ptr<Edge>> edgesWithoutInverse = filterInverses(edges);

    // Find index of filtered edge in the original graph edges
    vector<int> edgeIndices;
    for (const auto& edge : edgesWithoutInverse) {
        for (int i = 0; i < graph.edges.size(); ++i) {
            if (graph.edges[i] == edge) {
                edgeIndices.push_back(i);
                break;
            }
        }
    }
    set validIndices(edgeIndices.begin(), edgeIndices.end());

    // Convert generators to permutation vectors in reverse order
    vector<vector<int>> perms;
    for (int i = parsedGenerators.size() - 1; i >= 0; --i) {
        perms.push_back(cyclesToPermutation(parsedGenerators[i], graphSize));
    }

    // Identity permutation
    vector<int> id(graphSize);
    iota(id.begin(), id.end(), 0);

    // BFS to build automorphism group
    set<vector<int>> seen;
    queue<vector<int>> q;
    seen.insert(id);
    q.push(id);

    set<vector<shared_ptr<Edge>>> validEdgeMappings;
    auto start = chrono::high_resolution_clock::now();

    while (!q.empty() && validEdgeMappings.size() < maxAutomorphisms) {
        vector<int> current = q.front(); q.pop();

        if (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() > 2000) {
            break;
        }

        // Check if the current permutation is valid, only mapping edgesWithoutInverse to edgesWithoutInverse
        vector<shared_ptr<Edge>> currentMapping;
        bool isGood = true;
        for (int i = 0; i < edgesWithoutInverse.size(); ++i) {
            int mappedIndex = current[edgeIndices[i]];
            if (validIndices.contains(mappedIndex)) {
                currentMapping.push_back(graph.edges[mappedIndex]);
            } else {
                isGood = false;
                break;
            }
        }
        if (isGood) {
            validEdgeMappings.insert(currentMapping);
        }

        for (const auto& gen : perms) {
            vector<int> next = applyPermutation(gen, current);
            if (seen.insert(next).second) {
                q.push(next);
            }
        }
    }

    vector<vector<shared_ptr<Edge>>> result;
    // Filter out automorphisms that map semi edges to loops and vice versa
    // Also remove the first automorphism (identity)
    for (const auto& mapping : validEdgeMappings) {
        if (mapping == edgesWithoutInverse) {
            continue; // Skip identity
        }
        bool isValid = true;
        for (int i = 0; i < mapping.size(); ++i) {
            auto e1 = mapping[i];
            auto e2 = edgesWithoutInverse[i];
            if (e1->start == e1->end && e2->start == e2->end &&
                    ((e1->reverseEdge != nullptr && e2->reverseEdge == nullptr) ||
                    (e1->reverseEdge == nullptr && e2->reverseEdge != nullptr))){
                isValid = false;
                break;
            }
        }
        if (isValid) {
            result.push_back(mapping);
        }
    }

    // Print edges
    // for (auto e : edgesWithoutInverse) {
    //     clog << e->start << " " << e->end << " " << e->reverseEdge << endl;
    // }

    // Print automorphisms as integers indexed according to graph.edges
    // for (const auto& mapping : result) {
    //     vector<int> automorphism(mapping.size());
    //     for (int i = 0; i < mapping.size(); ++i) {
    //         auto e = mapping[i];
    //         int index = i;
    //         for (int j = 0; j < graph.edges.size(); ++j) {
    //             if (graph.edges[j] == e) {
    //                 index = j;
    //                 break;
    //             }
    //         }
    //         automorphism[i] = index;
    //     }
    //     clog << "Automorphism: ";
    //     for (int i : automorphism) {
    //         clog << i << " ";
    //     }
    //     clog << endl;
    // }

    return result;
}


Group cyclicGroup(int n) {
    vector<vector<int>> multiplicationTable;

    for (int i = 0; i < n; ++i) {
        vector<int> row;
        for (int j = 0; j < n; ++j) {
            row.push_back((i + j) % n);
        }
        multiplicationTable.push_back(row);
    }

    return Group(multiplicationTable);
}