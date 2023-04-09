
#include <string>
#include <iostream>
#include <pqxx/pqxx>
#include "HTTPRequest.hpp"

using namespace std;
using namespace pqxx;


vector<int> fetch(string& url, vector<pair<string, string>>& v){

    // make body of the post request
    if (v.empty()){
        throw invalid_argument("Empty vector is not allowed.");
    }
    string body = "{ \"locations\": [";

    for(int i=0; i<v.size()-1; i++){
        body += "{ \"longitude\": " + v[i].first + ", \"latitude\": " + v[i].second + "}";
        body += ',';
    }
    body += "{ \"longitude\": " + v.back().first + ", \"latitude\": " + v.back().second + "}";
    body += "]}";

    // make request and get response
    http::Response response;
    try
    {
        http::Request request{"http://192.168.1.20:1000/api/v1/lookup"};
        response = request.send("POST", body, {
                {"Content-Type", "application/json"}
        });
    }
    catch (const std::exception& e)
    {
        std::cerr << "Request failed, error: " << e.what() << '\n';
    }

    string json = std::string{response.body.begin(), response.body.end()};

    // parse the response
    vector<int> elevations;
    size_t i = json.find('[');
    size_t i2 = json.rfind(']');
    int count = 0;
    int num = 0;
    while (++i < i2){
        char ch = json[i];
        if (ch == '{'){
            count = 0;
            num = 0;
        } else if (ch == ':'){
            count++;
        } else if (count == 3 && '0' <= ch && ch <= '9'){
            num = num * 10 + (ch - '0');
        } else if (ch == '}'){
            elevations.push_back(num);
        }
    }

    return elevations;
}


int _processNodes(
        vector<long>& node_ids,
        vector<pair<string, string>>& points,
        string& fetchUrl,
        stream_to& stream2
        ){
    int count = 0;
    tuple<long, int> row2;
    vector<int> elevations = fetch(fetchUrl, points);
    if (node_ids.size() != elevations.size()){
        throw runtime_error("node_ids and elevation does not match in size.");
    }
    if (node_ids.size() == 0){
        return 0;
    }
    for (int i=0; i<node_ids.size(); i++){
        count += 1;
        row2 = {node_ids[i], elevations[i]};
        stream2 << row2;
    }
    return count;
}

pair<int, int> _processEdges(
        vector<long>& rows,
        vector<int>& edgeGeomSizes,
        vector<pair<string, string>>& points,
        string& fetchUrl,
        stream_to& stream2
        ){
    if (rows.size() != edgeGeomSizes.size()){
        throw runtime_error("Number of edges and number of geom sizes doesn't match.");
    }

    if (rows.size() == 0){
        return {0, 0};
    }

    int total(0);
    for (int size: edgeGeomSizes){
        total += size;
    }

    if (total != points.size()){
        throw runtime_error("Number of points does not match sum of geom sizes.");
    }

    tuple<long, string, int, int, string> row2;
    vector<int> elevations = fetch(fetchUrl, points);

    if (total != elevations.size()){
        throw runtime_error("Number of elevations does not match");
    }

    int count = 0;
    for (int i=0; i<edgeGeomSizes.size(); i++){
        int prevElevation = elevations[count++];
        int size = edgeGeomSizes[i];
        int elevationDiff(0);
        for (int j=1; j<size; j++){
            int curElevation = elevations[count++];
            elevationDiff += abs(curElevation - prevElevation);
        }
        tuple<long, int> row2 = {rows[i], elevationDiff};
        stream2 << row2;
    }
    return {edgeGeomSizes.size(), total};
}

