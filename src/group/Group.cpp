#include "Group.h"

#include <climits>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "GroupParser.h"

Group::Group(const vector<vector<int>>& multiplicationTable, const vector<vector<int>>& orbits, const vector<vector<int>>& automorphisms) {
    this->multiplicationTable = multiplicationTable;
    this->orbits = orbits;
    this->automorphisms = automorphisms;

    // orbitsSets
    for (auto &orbit : orbits) {
        set<int> orbitSet;
        for (int element : orbit) {
            orbitSet.insert(element);
        }
        this->orbitsSets.push_back(orbitSet);
    }

    // inverse
    for (int element = 0; element < size(multiplicationTable); element++) {
        for (int i = 0; i < multiplicationTable.size(); i++) {
            if (multiplicationTable[element][i] == 0) {
                this->inverse.push_back(i);
            }
        }
    }

    // powerToIdentity
    for (int element = 0; element < size(multiplicationTable); element++) {
        set<int> seen;
        int current = element;
        int power = 1;
        while (current != 0) {
            power++;
            seen.insert(current);
            int next = multiplicationTable[current][element];
            if (seen.contains(next)) {
                this->powerToIdentity.push_back(INT_MAX);
                cerr << "Power to identity is infinite for element " << element << " (this should never happen)" << endl;
                break;
            }
            current = next;
        }
        this->powerToIdentity.push_back(power);
        if (power == 2 || power == 1) {
            this->involutionCount++;
        }
    }
}

int Group::main_Group() {
    Group group = fromPath("../Groups/Z2xZ48.txt");

    // Print multiplication table
    clog << "Multiplication table:" << endl;
    for (auto &row : group.multiplicationTable) {
        for (int i : row) {
            clog << i << " ";
        }
        clog << endl;
    }

    return 0;
}

void Group::parseFromStdIn() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    string line;
    string groupline;
    int i = 0;
    while(getline(cin,line))
    {
        if (line == "####")
        {
            i++;
            Group group = fromString(groupline);
            clog << "Group " << i << ":" << group.multiplicationTable.size() << endl;
            groupline = "";
        } else {
            groupline += line;
        }
    }
}
