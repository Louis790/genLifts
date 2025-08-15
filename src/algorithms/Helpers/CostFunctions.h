#ifndef COSTFUNCTIONS_H
#define COSTFUNCTIONS_H

#include <functional>

#include "../../graph/Graph.h"
#include "../../group/Group.h"
#include "./Lift.h"
#include "RNG.h"

struct CostFunction {
    std::function<double(Graph&, const Group&, const shared_ptr<Edge>& edge)> func;
    std::string name;
    std::vector<std::string> paramNames;
    std::vector<double> paramValues;

    double operator()(Graph& graph, const Group& group, const shared_ptr<Edge>& edge) const {
        return func(graph, group, edge);
    }

    std::string getName() const {
        return name;
    }

    std::string getParams() const {
        std::string result = "{";
        for (size_t i = 0; i < paramNames.size(); ++i) {
            result += paramNames[i] + "=" + std::to_string(paramValues[i]);
            if (i + 1 < paramNames.size()) result += ", ";
        }
        result += "}";
        return result;
    }
};

double averageSimplifiedShortestCycles(Graph& graph, const Group& group, int amountOfCycles, const vector<int>& weights = {});
double averageCycleLength(Graph& graph, const Group& group, int ignoreCount);
double averageFundamentalCycleLength(Graph& graph, const Group& group, int ignoreCount);
double monteCarloAverageCycleLength(Graph& graph, const Group& group, int N);
double voltageDiversity(Graph& graph, const Group& group);
double walkSampler(const Group& group, Graph& graph, const shared_ptr<Edge>& edge, int minGirth, int N);
double walkRegularity(Graph &graph, const Group &group);

inline CostFunction createSimplifiedShortestCyclesCostFunction(int amountOfCycles, const std::vector<int>& weights = {}) {
    return CostFunction{
        [amountOfCycles, weights](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return averageSimplifiedShortestCycles(graph, group, amountOfCycles, weights);
        },
        "SimplifiedShortestCycles",
        {"amountOfCycles", "weights"},
        {static_cast<double>(amountOfCycles), static_cast<double>(weights.size())} // Simplify weights to size as an example
    };
}

inline CostFunction createAverageCycleLengthCostFunction(int ignoreCount) {
    return CostFunction{
        [ignoreCount](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return averageCycleLength(graph, group, ignoreCount);
        },
        "AverageCycleLength",
        {"ignoreCount"},
        {static_cast<double>(ignoreCount)}
    };
}

inline CostFunction createAverageFundamentalCycleLengthCostFunction(int ignoreCount) {
    return CostFunction{
        [ignoreCount](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return averageFundamentalCycleLength(graph, group, ignoreCount);
        },
        "AverageFundamentalCycleLength",
        {"ignoreCount"},
        {static_cast<double>(ignoreCount)}
    };
}

inline CostFunction createMonteCarloCycleLengthCostFunction(int N) {
    return CostFunction{
        [N](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return monteCarloAverageCycleLength(graph, group, N);
        },
        "MonteCarloCycleLength",
        {"N"},
        {static_cast<double>(N)}
    };
}

inline CostFunction createRandomCostFunction() {
    return CostFunction{
        [](Graph&, const Group&, const shared_ptr<Edge>&) {
            return -getRandomInt(200000, 300000);
        },
        "Random",
        {},
        {}
    };
}

inline CostFunction createVoltageDiversityCostFunction() {
    return CostFunction{
        [](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return voltageDiversity(graph, group);
        },
        "VoltageDiversity",
        {},
        {}
    };
}

inline CostFunction createWalkRegularityCostFunction() {
    return CostFunction{
        [](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return walkRegularity(graph, group);
        },
        "WalkRegularity",
        {},
        {}
    };
}

inline CostFunction createWalkSamplerCostFunction(int N, int minGirth) {
    return CostFunction{
        [N, minGirth](Graph& graph, const Group& group, const shared_ptr<Edge>& edge) {
            return walkSampler(group, graph, edge, minGirth, N);
        },
        "WalkSampler",
        {"N", "minGirth"},
        {static_cast<double>(N), static_cast<double>(minGirth)}
    };
}

set<set<pair<int, int>>> getCycles(Graph& graph);
vector<vector<int>> getFundamentalCycles(Graph& graph);



#endif //COSTFUNCTIONS_H
