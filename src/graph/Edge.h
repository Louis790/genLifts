#ifndef EDGE_H
#define EDGE_H
#include <memory>


class Edge {
public:
    Edge(const int start, const int end) {
        this -> start = start;
        this -> end = end;
    }

    void setVoltage(const int voltage, const int reverseVoltage) {
        this -> voltage = voltage;
        if (reverseEdge != nullptr) {
            reverseEdge -> voltage = reverseVoltage;
        }
    }

    int start;
    int end;
    int voltage = -1;
    bool isPartOfTree = false;
    std::shared_ptr<Edge> reverseEdge = nullptr;
};



#endif //EDGE_H
