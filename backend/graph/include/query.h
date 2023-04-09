#include <pqxx/pqxx>


#include "crow.h"


using namespace std;
using namespace pqxx;
//

enum MapConfig {
    all,
};


void prepare_nodes_query(connection_base &C, MapConfig mapConfig);

void prepare_edges_query(connection_base &C, MapConfig mapConfig);

void prepare_restrictions_query(connection_base &C, MapConfig mapConfig);

vector<pair<double, double>> parseLineString(string &lineString);

unordered_map<long, vector<pair<double, double>>> allEdges(string connectionString);

