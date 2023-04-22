
#include <queue>
#include <unordered_set>

using namespace std;

extern vector<string> parseCsvLine(const string &line);

template<typename S>
extern auto select_random(const S &s, size_t n);

template<typename T>
extern void fisherYatesShuffle(
        vector<T> v,
        size_t numRandom
);

template<typename T, typename H, typename E>
extern void chooseFromSet(
        unordered_set<T, H, E>& s,
        size_t numRandom
);

extern int roundToInt(double num);