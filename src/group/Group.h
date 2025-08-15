#ifndef GROUP_H
#define GROUP_H

#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;


class Group {
public:
    explicit Group(const vector<vector<int>>& multiplicationTable, const vector<vector<int>>& orbits = {}, const vector<vector<int>>& automorphisms = {});

    static int main_Group();

    static void parseFromStdIn();

    /**
     * The multiplication table of this group. multiplicationTable[a][b] is the result of a*b.
     */
    vector<vector<int>> multiplicationTable;

    /**
     * The automorphisms of this group. Automorphisms[i][j] is what automorphism i maps the j-th element to.
     */
    vector<vector<int>> automorphisms;

    /**
    * The orbits of this group.
    */
    vector<vector<int>> orbits;

     /**
     * The orbits of this group stored as sets.
     */
     vector<set<int>> orbitsSets;

    /**
     * The inverse of each element. inverse[a] is the element b such that a*b = 0.
     */
    vector<int> inverse;

    /**
     * The power to identity of each element. powerToIdentity[a] is the smallest n such that a^n = 0.
     */
    vector<int> powerToIdentity;

    /**
     * The number of involutions in the group. (g = g^-1)
     */
    int involutionCount = 0;

    string name = "Unnamed Group";
};



#endif //GROUP_H
