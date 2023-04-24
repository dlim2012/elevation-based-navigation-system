#include <pqxx/pqxx>


#include "crow.h"


using namespace std;
using namespace pqxx;
//

enum HighwayConfig {
    all_highways,
    motorway,
    motorway_without_restrictions,
    trunk,
    primary,
    secondary,
    tertiary,
    local,
    cycling,
    hiking,
    pedestrian,
    cycling_with_restrictions,
};

enum LocationConfig {
    US_mainland,
    Alaska,
    Hawaii
};



void prepare_nodes_query(connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig);

void prepare_edges_query(connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig);

void prepare_restrictions_query(connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig);

void prepare_edges_geom_query(connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig);

vector<pair<double, double>> parseLineString(string &lineString);


