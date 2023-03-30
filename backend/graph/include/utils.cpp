#include <iostream>
#include <fstream>
#include <vector>
#include <queue>

using namespace std;

vector<string> parseCsvLine(const string& line){
    vector<string> res;
    string word;
    for (char ch: line){
        if (ch == ','){
            res.push_back(word);
            word.clear();
        } else {
            word.push_back(ch);
        }
    }
    if (!word.empty()){
        res.push_back(word);
    }
    return res;
}
