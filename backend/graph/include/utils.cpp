#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <unordered_set>

using namespace std;

vector<string> parseCsvLine(const string &line) {
    vector<string> res;
    string word;
    for (char ch: line) {
        if (ch == ',') {
            res.push_back(word);
            word.clear();
        } else {
            word.push_back(ch);
        }
    }
    if (!word.empty()) {
        res.push_back(word);
    }
    return res;
}

template<typename S>
extern auto select_random(const S &s, size_t n) {
    auto it = std::begin(s);
    // 'advance' the iterator n times
    std::advance(it, n);
    return it;
}

template<typename T>
void fisherYatesShuffle(
        vector<T> v,
        size_t numRandom
    ){
    if (v.size() < numRandom){
        throw runtime_error("Number of selection is bigger than the vector size");
    }
    size_t left = v.size();
    if (numRandom < v.size() / 2){
        auto begin = v.begin();
        while (numRandom--) {
            swap(begin++, advance(begin, rand() % left--));
        }
    } else {
        auto end = v.rend();
        numRandom = left - numRandom;
        while (numRandom--){
            swap(end++, advance(end, rand() % left--));
        }
    }
}

template<typename T, typename H, typename E>
void chooseFromSet(
        unordered_set<T, H, E>& s,
        size_t numRandom
        ){
    if (numRandom > s.size())
        throw runtime_error("Number of selection is bigger than the set size");
    while (s.size() > numRandom){
        s.erase(advance(s.begin(), rand() % s.size()));
    }
}

