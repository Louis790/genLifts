#include "VertexEdgeGirthRegular.h"

#include <functional>
#include <iostream>
#include <queue>
#include <stack>

map<pair<int, int>, pair<int, int>> calcNumberShortestPathsCache = {};

bool hasCycle(const Graph& graph, vector<bool>& vis, const int current, const int start, const int remaining, vector<int>& cycle) {
    if (remaining == 0) {
        for (const int neigh : graph.adjacency[current]) {
            if (neigh == start) return true;
        }
        return false;
    }
    for (int neigh : graph.adjacency[current]) {
        if (vis[neigh] || neigh < start) continue;
        vis[neigh] = true;
        cycle[remaining] = neigh;
        if (hasCycle(graph, vis, neigh, start, remaining - 1, cycle)) return true;
        vis[neigh] = false;
    }
    return false;
};

bool hasNoCycleOfLength(const Graph& graph, const int length) {
    // Modified from https://github.com/JorikJooken/KGNoGPlus1Graphs/blob/main/Code/filterHasNoCycleOfLengthVeryLargeGraph.cpp
    const int n = graph.adjacency.size();
    vector vis(n, false);
    vector cycle(n, 0);

    for (int i = 0; i + length - 1 < n; i++) {
        vis.assign(n, false);
        vis[i] = true;
        cycle[0] = i;
        if (hasCycle(graph, vis, i, i, length - 1, cycle)) {
            return false;
        }
    }
    return true;
}


pair<int, int> calcNumberShortestPaths(Graph& graph, int u, int v, int girth, bool useCache) {
    // Modified from https://github.com/JorikJooken/edgeGirthRegularGraphs/blob/master/Code/veryLargeGraphs/calcLambdaRegular.cpp
    if (useCache && calcNumberShortestPathsCache.contains({u, v})) {
        return calcNumberShortestPathsCache[{u, v}];
    }

    int qStart=0;
    int qFinish=0;
    int n = size(graph.adjacency);
    vector<int> vertexQ(size(graph.adjacency));
    vector<int> distBFS1(size(graph.adjacency));
    vector<int> numberShortestPaths(size(graph.adjacency));
    for (int i = 0; i < size(graph.adjacency); i++) {
        vertexQ[i] = 0;
        distBFS1[i] = INT_MAX;
        numberShortestPaths[i] = 0;
    }

    vertexQ[qStart] = u;
    numberShortestPaths[u] = 1;
    distBFS1[u] = 0;

    while(qStart<=qFinish) {
        int now = vertexQ[qStart];
        if (now == v) {
            break;
        }
        qStart++;
        for(int neigh : graph[now]) {
            if(now==u && neigh==v) continue;
            if(now==v && neigh==u) continue;
            if(distBFS1[neigh]>n) {
                distBFS1[neigh] = distBFS1[now]+1;
                numberShortestPaths[neigh] = numberShortestPaths[now];
                qFinish++;
                vertexQ[qFinish]=neigh;
            }
            else if(distBFS1[neigh] == distBFS1[now]+1) {
                numberShortestPaths[neigh] += numberShortestPaths[now];
            }
        }
    }

    pair p = {distBFS1[v] , numberShortestPaths[v]};
    if (useCache) {
        calcNumberShortestPathsCache[{u, v}] = p;
    }
    return p;
}

pair<int, int> getBothGirthRegularLambda(Graph& graph, int girth) {
    int edgeLambda = getEdgeGirthRegularLambda(graph, girth, true);
    int vertexLambda = getVertexGirthRegularLambda(graph, girth, true);
    calcNumberShortestPathsCache.clear();
    return {edgeLambda, vertexLambda};
}

int getEdgeGirthRegularLambda(Graph& graph, int girth, bool useCache) {
    int currentLambda = -1;
    for (int start = 0; start < size(graph.adjacency); start++) {
        for (int end : graph.adjacency[start]) {
            auto p = calcNumberShortestPaths(graph, start, end, girth, useCache);
            if (p.first + 1 != girth) {
                return -1;
            }
            if (currentLambda == -1) {
                currentLambda = p.second;
            }
            if (currentLambda != p.second) {
                return -1;
            }
        }
    }

    return currentLambda;
}

int getVertexGirthRegularLambda(Graph& graph, int girth, bool useCache) {
    // modified from https://github.com/JorikJooken/vertexGirthRegularGraphs/blob/master/Code/veryLargeGraphs/calcLambdaRegular.cpp
    int lambda = -2;
    for (int start=0; start<size(graph.adjacency) && lambda != -1; start++) {
        long long sum = 0;
        vector<int>& neighsOfStart = graph.adjacency[start];
        for (const int end : neighsOfStart) {
            auto [pathLen, pathCount] = calcNumberShortestPaths(graph, start, end, girth, useCache);
            if(pathLen + 1 == girth) {
                sum += pathCount;
            }
        }
        if (sum % 2 != 0){
            cerr << "sum is odd, this should never happen" << endl;
        }
        sum /= 2;
        if (lambda == -2) {
            lambda = sum;
        } else if (lambda != sum) {
            lambda = -1;
        }
    }
    return lambda;
}
