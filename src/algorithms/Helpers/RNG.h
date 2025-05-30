#ifndef RNG_H
#define RNG_H

#include <random>
#include <vector>

using std::vector;

inline int seed = 7;
inline std::mt19937 rng(seed);

int getRandomInt(int min, int max);
template<typename T>
void shuffleVec(std::vector<T>& vec) {
    shuffle(vec.begin(), vec.end(), rng);
}
void resetWithSeed(int newSeed);

#endif //RNGESUS_H
