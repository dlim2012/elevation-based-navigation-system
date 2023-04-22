#include <unordered_map>
#include <iomanip>

#include "pqxx/pqxx"
#include "query.h"


using namespace std;
using namespace pqxx;

const string motorwayTypes = "(\n"
                             "\t\t'motorway', 'motorway_link'\n"
                             "\t)";

const string trunkTypes = "(\n"
                             "\t\t'motorway', 'motorway_link,'\n"
                          "\t\t'trunk', 'trunk_link'\n"
                             "\t)";

const string primaryTypes = "(\n"
                       "\t\t'motorway', 'motorway_link',\n"
                       "\t\t'trunk', 'trunk_link',\n"
                       "\t\t'primary', 'primary_link'\n"
                       "\t)";

const string secondaryTypes = "\t(\n"
                             "\t\t'motorway', 'motorway_link',\n"
                             "\t\t'trunk', 'trunk_link',\n"
                             "\t\t'primary', 'primary_link',\n"
                             "\t\t'secondary', 'secondary_link'\n"
                             "\t)";

const string tertiaryTypes = "\t(\n"
                             "\t\t'motorway', 'motorway_link',\n"
                             "\t\t'trunk', 'trunk_link',\n"
                             "\t\t'primary', 'primary_link',\n"
                             "\t\t'secondary', 'secondary_link',\n"
                             "\t\t'tertiary', 'tertiary_link'\n"
                             "\t)";

const string localTypes = "\t(\n"
                          "\t\t'motorway', 'motorway_link',\n"
                          "\t\t'trunk', 'trunk_link',\n"
                          "\t\t'primary', 'primary_link',\n"
                          "\t\t'secondary', 'secondary_link',\n"
                          "\t\t'tertiary', 'tertiary_link',\n"
                          "\t\t'unclassified', 'residential\n"
                          "\t)";

const string pedestrianTypes = "\t('footway', 'pedestrian', 'living_street', 'path', 'track')";

const string& highwayConfigToType(HighwayConfig highwayConfig){
    switch (highwayConfig){
        case motorway:
        case motorway_without_restrictions:
            return motorwayTypes;
        case trunk:
            return trunkTypes;
        case primary:
            return primaryTypes;
        case secondary:
            return secondaryTypes;
        case tertiary:
            return tertiaryTypes;
        case local:
            return localTypes;
        case pedestrian:
            return pedestrianTypes;
        default:
            throw runtime_error("highwayConfig not found");
    }
}

void prepare_nodes_query(connection_base &C, HighwayConfig highwayConfig) {
    if (highwayConfig == all_highways) {
        C.prepare(
                "nodes",
                "SELECT id, ST_X(geom) AS lon, ST_Y(geom) AS lat, elevation FROM nodes;"
        );
    } else if (highwayConfig == cycling || highwayConfig == cycling_with_restrictions) {
        C.prepare(
                "nodes",
                "SELECT id, ST_X(geom) AS lon, ST_Y(geom) AS lat, elevation FROM nodes"
                "\tWHERE cycling = 'accept';"
        );
    } else if (highwayConfig == hiking) {
        C.prepare(
                "nodes",
                "SELECT id, ST_X(geom) AS lon, ST_Y(geom) AS lat, elevation FROM nodes"
                "\tWHERE hiking = 'accept';"
        );
    } else if (highwayConfig == pedestrian){
        C.prepare(
                "nodes",
                "SELECT id, ST_X(geom) AS lon, ST_Y(geom) AS lat, elevation FROM nodes"
                "\tWHERE type != 'motorway' AND type != 'motorway_link';"
        );
    } else {
        C.prepare(
                "nodes",
                "SELECT id, ST_X(geom) AS lon, ST_Y(geom) AS lat, elevation FROM nodes"
                "\tWHERE type in " + highwayConfigToType(highwayConfig) + ";"
        );
    }
}


void prepare_edges_query(connection_base &C, HighwayConfig highwayConfig) {
    if (highwayConfig == all_highways){
        C.prepare(
                "edges",
                "SELECT id, u, v, length, elevation, oneway FROM edges;"
        );
    } else if (highwayConfig == cycling || highwayConfig == cycling_with_restrictions) {
        C.prepare(
                "edges",
                "SELECT id, u, v, length, elevation, oneway FROM edges"
                "\tWHERE cycling = 'accept';"
        );
    } else if (highwayConfig == hiking) {
        C.prepare(
                "edges",
                "SELECT id, u, v, length, elevation, oneway FROM edges"
                "\tWHERE hiking = 'accept';"
        );
    } else {
        C.prepare(
                "edges",
                "SELECT id, u, v, length, elevation, oneway FROM edges"
                "\tWHERE type in " + highwayConfigToType(highwayConfig) + ";"
        );
    }
}

void prepare_restrictions_query(connection_base &C, HighwayConfig highwayConfig) {
    if (highwayConfig == hiking || highwayConfig == cycling || highwayConfig == motorway_without_restrictions) {
        return;
    }
    C.prepare(
            "restrictions",
            "SELECT type_enum, from_id, via_id, to_id FROM restrictions\n"
            "     WHERE type_enum IS NOT NULL"
    );
}

void prepare_edges_geom_query(connection_base &C, HighwayConfig highwayConfig) {
    if (highwayConfig == all_highways){
        C.prepare(
                "edges_geom",
                "SELECT id, ST_AsText(geom, 10), u, v from edges;"
        );
    } else if (highwayConfig == cycling || highwayConfig == cycling_with_restrictions) {
        C.prepare(
                "edges_geom",
                "SELECT id, ST_AsText(geom, 10), u, v from edges"
                "\tWHERE cycling = 'accept';"
        );
    } else if (highwayConfig == hiking) {
        C.prepare(
                "edges_geom",
                "SELECT id, ST_AsText(geom, 10), u, v from edges"
                "\tWHERE hiking = 'accept';"
        );
    } else {
        C.prepare(
                "edges_geom",
                "SELECT id, ST_AsText(geom, 10), u, v from edges\n"
                "\tWHERE type in " + highwayConfigToType(highwayConfig) + ";"
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

