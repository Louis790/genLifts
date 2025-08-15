#include <iostream>
#include <vector>
#include <string>
#include <cstring>   // for std::memcmp
#include <cstdint>   // for fixed-size integer types
#include <fstream>
#include <sstream>
#include <stdexcept> // for std::runtime_error

#include "Graph.h"
#include "GraphParser.h"
#include "../nautyAndMultigraph/multigraph.h"

using std::cout;
using std::cerr;
using std::cin;

int getGraph6NumberOfVertices(string graphString)
{
    if(graphString.empty()){
        printf("Error: String is empty.\n");
        abort();
    }
    else if((graphString[0] < 63 || graphString[0] > 126) && graphString[0] != '>') {
        printf("Error: Invalid start of graphstring.\n");
        abort();
    }

    int index = 0;
    if (graphString[index] == '>') { // Skip >>graph6<< header.
        index += 10;
    }

    if(graphString[index] < 126) { // 0 <= n <= 62
        return (int) graphString[index] - 63;
    }

    else if(graphString[++index] < 126) {
        int number = 0;
        for(int i = 2; i >= 0; i--) {
            number |= (graphString[index++] - 63) << i*6;
        }
        return number;
    }

    else if (graphString[++index] < 126) {
        int number = 0;
        for (int i = 5; i >= 0; i--) {
            number |= (graphString[index++] - 63) << i*6;
        }
        return number;
    }

    else {
        printf("Error: Format only works for graphs up to 68719476735 vertices.\n");
        abort();
    }
}

Graph parseGraph6String(string graphString) {
    int numberOfVertices = getGraph6NumberOfVertices(graphString);
    Graph graph;
    vector<int> emp;
    graph.adjacency.assign(numberOfVertices,emp);
    int startIndex = 0;
    if (graphString[startIndex] == '>') { // Skip >>graph6<< header.
        startIndex += 10;
    }
    if (numberOfVertices <= 62) {
        startIndex += 1;
    }
    else if (numberOfVertices <= MAXVERTICES) {
        startIndex += 4;
    }
    else {
        printf("Error: Program can only handle graphs with %d vertices or fewer.\n",MAXVERTICES);
        abort();
    }

    int currentVertex = 1;
    int sum = 0;
    for (int index = startIndex; index<graphString.size(); index++) {
        int i;
        for (i = custom_prev(graphString[index] - 63, 6); i != -1; i = custom_prev(graphString[index] - 63, i)) {
            while(5-i+(index-startIndex)*6 - sum >= 0) {
                sum += currentVertex;
                currentVertex++;
            }
            sum -= --currentVertex;
            int neighbour = 5-i+(index - startIndex)*6 - sum;
            graph.adjacency[currentVertex].push_back(neighbour);
            graph.adjacency[neighbour].push_back(currentVertex);
        }
    }

    return graph;
}

bool isLittleEndian(string s) {
    // littleEndian if ">>pregraph_code<<" or ">>pregraph_code le<<" is found
    // BigEndian if ">>pregraph_code be<<" is found
    // -> check a single char
    if (size(s) < 16 || s[0] != '>')
    {
        // cout << "Cannot find header, assuming little endian format" << endl;
        return true;
    }

    if (s[16] == 'l' || s[16] == '<') {
        return true;
    }
    if (s[16] == 'b') {
        return false;
    }
    cerr << "Error: Invalid pregraph string." << endl;
    abort();
}

uint8_t readByte(const string& data, size_t& pos) {
    return static_cast<uint8_t>(data[pos++]);
}

// Function to read 2-byte (16-bit) values from a string in little or big endian
uint16_t readTwoBytes(const string& data, size_t& pos, bool littleEndian) {
    uint16_t byte1 = readByte(data, pos);
    uint16_t byte2 = readByte(data, pos);
    if (littleEndian) {
        return (byte2 << 8) | byte1;
    } else {
        return (byte1 << 8) | byte2;
    }
}

vector<Graph> parsePregraphs(string& s) {
    bool littleEndian = isLittleEndian(s);
    vector<Graph> graphs;

    // remove header if present
    if (s[0] == '>') {
        s = s.substr(18);
        if (s[0] == '<') {
            s = s.substr(2);
        }
    }

    size_t pos = 0;
    size_t vertexCount = 0;
    bool largeGraph = false;

    while (pos < s.size()) {
        uint8_t firstByte = readByte(s, pos);
        if (firstByte == 0) {
            // Large graph, use 2-byte vertex counts
            vertexCount = readTwoBytes(s, pos, littleEndian);
            largeGraph = true;
        } else
        {
            // Small graph, use 1-byte vertex counts
            vertexCount = firstByte;
        }

        Graph graph;
        for (size_t i = 0; i < vertexCount; ++i) {
            vector<int> neighbors;
            while (true) {
                int neighbor;
                if (largeGraph) {
                    neighbor = readTwoBytes(s, pos, littleEndian);
                } else {
                    neighbor = readByte(s, pos);
                }

                if ((!largeGraph && neighbor == 32) || (largeGraph && neighbor == 64) || neighbor == 0) {
                    break;
                }

                neighbors.push_back(neighbor - 1); // Adjust for 0-indexing
            }

            if (size(neighbors) == 0) {
                graph.adjacency.emplace_back();
            } else {
                graph.adjacency.push_back(neighbors);
            }
        }

        // Bidirectional graph
        for (int i = 0; i < vertexCount; ++i) {
            for (int neighbor : graph.adjacency[i]) {
                if (neighbor > i) {
                    graph.adjacency[neighbor].push_back(i);
                }
            }
        }
        graphs.push_back(graph);
    }

    return graphs;
}

void insertInto(char** l, const int index, const string& s) {
    const auto cstr = new char[s.length() + 1];
    strcpy(cstr, s.c_str());
    l[index] = cstr;
}

vector<Graph> parseMultiGraphs(vector<int>& ns){

    // if only a single node graph
    int sum = 0;
    for (int n : ns) {
        sum += n;
    }
    if (sum == 1 && ns[ns.size()-1] == 1) {
        Graph g;
        g.adjacency.emplace_back();
        return {g};
    }

    // call 'multigraph ns[0] ns[1] ... ns[n-1] t1 o'
    char* argv[size(ns) + 3];

    insertInto(argv, 0, "multigraph");
    for (int i = 0; i < size(ns); i++) {
        insertInto(argv, i + 1, to_string(ns[i]));
    }
    insertInto(argv, size(ns)+1, "t1");
    insertInto(argv, size(ns)+2, "o");

    reset();
    const auto structRes = run_multigraph_gen(size(ns) + 3, argv);
    if (structRes.graphs == -1) {
        return {};
    }

    string result;
    for (int i = 0; i < structRes.graphs; i++) {
        result += string(reinterpret_cast<char*>(structRes.data + i * structRes.size), structRes.size) + " ";
    }

    // // Write result to file
    // ofstream file("multigraphs.txt");
    // file << result;
    // file.close();

    vector<Graph> graphs = parsePregraphs(result);

    return graphs;
}

vector<vector<int>> parseOrbit(string& s) {
    // Orbits separated by ';'
    // a single orbit can contain multiple nodes separated by spaces or a range of nodes indicated by ':'. Example:
    // 0 1 3; 2; 4:6;
    vector<vector<int>> orbits;
    vector<int> orbit;
    stringstream ss(s);
    string part;

    while (getline(ss, part, ';')) {
        stringstream orbitStream(part);
        string nodeRange;

        orbit.clear();

        while (orbitStream >> nodeRange) {
            size_t colonPos = nodeRange.find(':');

            if (colonPos != string::npos) {
                // Range found, e.g., "4:6"
                int start = stoi(nodeRange.substr(0, colonPos));
                int end = stoi(nodeRange.substr(colonPos + 1));

                for (int i = start; i <= end; ++i) {
                    orbit.push_back(i);
                }
            } else {
                // Single node
                orbit.push_back(stoi(nodeRange));
            }
        }

        if (!orbit.empty()) {
            orbits.push_back(orbit);
        }
    }
    return orbits;
}
