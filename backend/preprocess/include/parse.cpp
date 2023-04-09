#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <pqxx/pqxx>

using namespace std;
using namespace pqxx;

pair<string, string> parsePoint(string& Point){

    size_t start = Point.find('(') + 1;
    size_t end = Point.find(')');

    size_t mid;

    for (int i=start; i<end; i++){
        if (Point[i] == ' '){
            mid = i;
            break;
        }
    }
    return {Point.substr(start, mid-start), Point.substr(mid+1, end-mid-1)};
}

vector<pair<string, string>> parseLineString(string& lineString){

    size_t start = lineString.find('(') + 1;
    size_t end = lineString.find(')');

    vector<pair<string, string>> res;

    pair<string, string> p;
    int prev(start);
    for (int i=start; i<end; i++){
        if (lineString[i] == ' '){
            p.first = lineString.substr(prev, i-prev);
            prev = i + 1;
        } else if (lineString[i] == ','){
            p.second = lineString.substr(prev, i-prev);
            res.push_back(p);
            prev = i + 1;
        }
    }
    p.second = lineString.substr(prev, end-prev);
    res.push_back(p);


    return res;
}
