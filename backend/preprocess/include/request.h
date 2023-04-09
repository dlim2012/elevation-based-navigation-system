#include <string>
#include <iostream>
#include <pqxx/pqxx>

#include "HTTPRequest.hpp"


using namespace std;
using namespace pqxx;


vector<int> fetch(string& url, vector<pair<string, string>>& v);
int _processNodes(
        vector<long>& node_ids,
        vector<pair<string, string>>& points,
        string& fetchUrl,
        stream_to& stream2
        );

pair<int, int> _processEdges(
        vector<long>& rows,
        vector<int>& edgeGeomSizes,
        vector<pair<string, string>>& points,
        string& fetchUrl,
        stream_to& stream2
        );