#include <unordered_map>
#include <iomanip>

#include "pqxx/pqxx"
#include "query.h"


using namespace std;
using namespace pqxx;


void prepare_nodes_query(connection_base &C, MapConfig mapConfig) {
    if (mapConfig == all) {
        C.prepare(
                "nodes",
                "SELECT id, ST_X(geom) AS lon, ST_Y(geom) AS lat, elevation FROM nodes"
        );
    }
}


void prepare_edges_query(connection_base &C, MapConfig mapConfig) {
    if (mapConfig == all) {
        C.prepare(
                "edges",
                "SELECT id, u, v, length, elevation, oneway FROM edges"
        );
    }
}

void prepare_restrictions_query(connection_base &C, MapConfig mapConfig) {
    if (mapConfig == all) {
        C.prepare(
                "restrictions",
                "SELECT type_enum, from_id, via_id, to_id FROM restrictions\n"
                "     WHERE type_enum IS NOT NULL"
        );
    }
}

vector<pair<double, double>> parseLineString(string &lineString) {

    size_t start = lineString.find('(') + 1;
    size_t end = lineString.find(')');

    vector<pair<double, double>> res;

    pair<double, double> p;
    int prev(start);
    for (int i = start; i < end; i++) {
        if (lineString[i] == ' ') {
            p.first = stod(lineString.substr(prev, i - prev));
            prev = i + 1;
        } else if (lineString[i] == ',') {
            p.second = stod(lineString.substr(prev, i - prev));
            res.push_back(p);
            prev = i + 1;
        }
    }
    p.second = stod(lineString.substr(prev, end - prev));
    res.push_back(p);

    return res;
}

unordered_map<long, vector<pair<double, double>>> allEdges(string connectionString) {
    pqxx::connection C(connectionString);
    pqxx::work w(C);
    pqxx::result r = w.exec("SELECT id, ST_AsText(geom, 15) from edges");
    unordered_map<long, vector<pair<double, double>>> umap;
    for (pqxx::row row: r) {
        string geom_text = row[1].c_str();
        umap[stol(row[0].c_str())] = parseLineString(geom_text);
    }
    return umap;
};
