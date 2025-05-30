#include "RNG.h"

#include <algorithm>
#include <random>

int getRandomInt(const int min, const int max) {
    std::uniform_int_distribution dist(min, max);
    return dist(rng);
}

void resetWithSeed(const int newSeed) {
    seed = newSeed;
    rng = std::mt19937(seed);
}
