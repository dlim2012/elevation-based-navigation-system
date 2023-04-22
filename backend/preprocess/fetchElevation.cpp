#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <pqxx/pqxx>

#include "include/parse.h"
#include "include/request.h"
#include "include/sqlQuery.h"
#include <iomanip>

using namespace std;
using namespace pqxx;

const int maxEdgesRowSelect = 100000;
const vector<string> majorTypes = {
        "motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link",
        "secondary", "secondary_link", "tertiary", "tertiary_link", "unclassified", "unclassified_link",
        "residential", "residential_link"
};
const vector<string> tagFiltering = {
        "reject", "conditional", "accept"
};

void addElevationToNodes(
        connection& C,
        connection& C2,
        string& fetchUrl,
        int maxFetch
){

    long numRows = getNumberOfRows(C, -1, "nodes");
    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS elevation INT;", false);
    timedExecution(C, "DROP TABLE IF EXISTS _nodes", false);
    timedExecution(C, "CREATE TABLE IF NOT EXISTS _nodes ("
                      "    id BIGINT PRIMARY KEY,"
                      "    elevation INT"
                      ")", false);

    cout << "Fetching elevation data for nodes and saving to the new table... (total " << numRows << " rows)" << endl;
    work tx1(C);
    work tx2(C2);

    result r =tx1.exec("SELECT id, ST_AsText(geom) FROM nodes;");
    tx1.commit();
//    stream_from stream1(tx1, "_nodes", vector<string>{"id", "geom_text", "geom"});
    tuple<long, string> row1;

    stream_to stream2(tx2, "_nodes", vector<string>{"id", "elevation"});


    long countRead(0), countWrote(0);
    long percentage(0);
    vector<long> node_ids;
    vector<pair<string, string>> points;
    vector<string> geoms;
    auto t0 = chrono::high_resolution_clock::now();
    for (row row: r){
//    while (stream1 >> row1){
        row1 = {stol(row[0].c_str()), row[1].c_str()};

        countRead += 1;

        node_ids.push_back(get<0>(row1));
        string gemo_text = get<1>(row1);
        pair<string, string> point = parsePoint(gemo_text);
        points.push_back(point);

        if (points.size() >= maxFetch){
            countWrote += _processNodes(node_ids, points, fetchUrl, maxFetch, stream2);
            node_ids.clear();
            points.clear();
            if (countWrote * 100 >= (percentage+1) * numRows){
                while (countWrote * 100 >= (percentage+1) * numRows)
                    percentage++;
                auto t1 = chrono::high_resolution_clock::now();
                auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << setw(8) << right << countWrote << " rows / " << setw(8) << right << numRows << " rows (time: " << right << setw(8) << t << " seconds, expected time left: "
                     << right << setw(8) << t / countWrote * (numRows - countWrote) << +" seconds)" << endl;

            }
        }

    }
    if (countWrote < numRows){
        countWrote += _processNodes(node_ids, points, fetchUrl, maxFetch, stream2);
        if (countWrote * 100 >= (percentage+1) * numRows){
            while (countWrote * 100 > percentage * numRows)
                percentage++;
            auto t1 = chrono::high_resolution_clock::now();
            auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << setw(8) << right << countWrote << " rows / " << setw(8) << right << numRows << " rows (time: " << right << setw(8) << t << " seconds, expected time left: "
                 << right << setw(8) << t / countWrote * (numRows - countWrote) << +" seconds)" << endl;

        }
    }

//    stream1.complete();
    stream2.complete();
    tx2.commit();

    cout << countRead << " " << countWrote << endl;
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "(time: " << t << " seconds)" << endl << endl;

    timedExecution(C, "UPDATE nodes SET elevation = _nodes.elevation FROM _nodes WHERE nodes.id = _nodes.id", false);
    timedExecution(C, "DROP TABLE _nodes", false);
}





void addElevationToEdges(
        connection& C,
        connection& C2,
        string& fetchUrl,
        int maxFetch
){

    long numRows = getNumberOfRows(C, -1, "edges");
    long numPoints = getNumberOfPoints(C, -1, "edges");
    timedExecution(C, "ALTER TABLE edges ADD COLUMN IF NOT EXISTS elevation INT;", false);
//    timedExecution(C, "DROP TABLE IF EXISTS _edges", false);
    timedExecution(C, "CREATE TABLE IF NOT EXISTS _edges (\n"
                      "    id BIGINT PRIMARY KEY,\n"
                      "    elevation INT\n"
                      ");", false);

    cout << "Edges table: " << numRows << " rows, " << numPoints << " points" << endl;

    long countEdgeRead(0), countEdgeWrote(0);
    long countPointRead(0), countPointProcessed(0), countPointCur(0);
    long percentage(0);
    vector<long> ids;
    vector<int> edgeGeomSizes;
    vector<pair<string, string>> points;
    long id = getLastId(C, "_edges");

    long numRowsLeft = getNumberOfRows(C, id, "edges");
    long numPointsLeft = getNumberOfPoints(C, id, "edges");

    cout << "Fetching elevation data for edges and saving to the new table... (todo: " << numRowsLeft << " rows, " << numPointsLeft << " points)" << endl;


    auto t0 = chrono::high_resolution_clock::now();
    while (true) {
        work w(C);
        result r = w.exec("SELECT count(*) FROM edges WHERE id > " + to_string(id) + ";");
        w.commit();

        work tx1(C);
        if (stol(r[0][0].c_str()) == 0)
            break;
        r = tx1.exec("SELECT id, ST_AsText(geom) FROM edges WHERE id > " + to_string(id) +
                " ORDER BY id limit " + to_string(maxEdgesRowSelect) + ";");
        tx1.commit();
        //    stream_from stream1(tx1, "edges", vector<string>{"id", "geom_text"});
        //    tuple<long, string> row1;

        work tx2(C2);
        stream_to stream2(tx2, "_edges", vector < string > {"id", "elevation"});

        for (row row: r) {
            //    while (stream1 >> row1){
            countEdgeRead++;
            id = stol(row[0].c_str());
            //        long id = get<0>(row1);
            string geom_text = row[1].c_str();
            //        string geom_text = get<1>(row1);

            vector <pair<string, string>> linePoints = parseLineString(geom_text);
            if (linePoints.size() < 2) {
                cout << "Encountered an edge with less than two points. Ignoring this edge." << endl;
                continue;
            }

            if (countPointCur + linePoints.size() > maxFetch) {
                pair<int, int> counts = _processEdges(ids, edgeGeomSizes, points, fetchUrl, maxFetch, stream2);
                countEdgeWrote += counts.first;
                countPointProcessed += counts.second;
                ids.clear();
                edgeGeomSizes.clear();
                points.clear();
                countPointCur = 0;

                if (countEdgeWrote * 100 >= (percentage+1) * numRowsLeft) {
                    while (countEdgeWrote * 100 >= (percentage+1) * numRowsLeft)
                        percentage++;
                    auto t1 = chrono::high_resolution_clock::now();
                    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                    cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << right << setw(8)
                         << countEdgeWrote << " rows / " << setw(8) << right << numRowsLeft << " rows, "
                         << setw(8) << right << countPointProcessed << " points / " << setw(8) << right << numPointsLeft
                         << " points (time: " << right
                         << setw(8) << t << " seconds, expected time left: " << right
                         << setw(8) << t / countPointProcessed * (numPointsLeft - countPointProcessed) << +" seconds)"
                         << endl;

                }
            }

            ids.push_back(id);
            edgeGeomSizes.push_back(linePoints.size());
            points.insert(points.end(), make_move_iterator(linePoints.begin()), make_move_iterator(linePoints.end()));
            countPointRead += linePoints.size();
            countPointCur += linePoints.size();

        }
        //    stream1.complete();
        stream2.complete();
        tx2.commit();
    }
    if (!ids.empty()) {
        work tx3(C2);
        stream_to stream3(tx3, "_edges", vector < string > {"id", "elevation"});
        pair<int, int> counts = _processEdges(ids, edgeGeomSizes, points, fetchUrl, maxFetch, stream3);
        countEdgeWrote += counts.first;
        countPointProcessed += counts.second;
        if (countEdgeWrote * 100 > (percentage+1) * numRowsLeft) {
            while (countEdgeWrote * 100 > (percentage+1) * numRowsLeft)
                percentage++;
            auto t1 = chrono::high_resolution_clock::now();
            auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << right << setw(8)
                 << countEdgeWrote << " rows / " << setw(8) << right << numRowsLeft << " rows, "
                 << setw(8) << right << countPointProcessed << " points / " << setw(8) << right << numPointsLeft
                 << " points (time: " << right
                 << setw(8) << t << " seconds, expected time left: " << right
                 << setw(8) << t / countPointProcessed * (numPointsLeft - countPointProcessed) << +" seconds)" << endl;

        }
        stream3.complete();
        tx3.commit();
    }

    cout << countEdgeRead << " " << countEdgeWrote << " " << countPointRead << " " << countPointProcessed << endl;

    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "(time: " << t << " seconds)" << endl << endl;

    timedExecution(C, "UPDATE edges SET elevation = _edges.elevation FROM _edges WHERE edges.id = _edges.id", false);
    timedExecution(C, "DROP TABLE _edges;", false);
}

void processNodes(connection& C, connection& C2, string& fetchUrl, int maxFetch){
    // Add geom_text temporarily to use the stream functionality of libpqxx

    timedExecution(C, "CREATE TABLE IF NOT EXISTS nodes (\n"
                      "     id bigint primary KEY,\n"
                      "     geom geometry(POINT, 4326),\n"
                      "     type text\n"
                      ");", false);
    timedExecution(C, "INSERT INTO nodes(id, geom)\n"
                      "    SELECT nodes[ARRAY_LOWER(nodes, 1)], ST_PointN(geom, 1) FROM edges\n"
                      "    ON CONFLICT DO NOTHING;", false);
    timedExecution(C, "INSERT INTO nodes(id, geom)\n"
                      "    SELECT nodes[ARRAY_UPPER(nodes, 1)], ST_PointN(geom, -1) FROM edges\n"
                      "    ON CONFLICT DO NOTHING;", false);

    for (auto it=majorTypes.rbegin(); it != majorTypes.rend(); it++) {
        string type = *it;
        timedExecution(C, "UPDATE nodes SET type = '" + *it + "' WHERE id IN (\n"
                          "\tSELECT nodes[ARRAY_LOWER(nodes, 1)] FROM edges WHERE type = '" + *it + "'\n"
                          "\tUNION\n"
                          "\tSELECT nodes[ARRAY_UPPER(nodes, 1)] FROM edges WHERE type = '" + *it + "'\n"
                          "\t);", false);
    }

    // Fetch elevation from Open-Elevation API and save to the new table;
    addElevationToNodes(C, C2, fetchUrl, maxFetch);

    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS relation_ids bigint[];", false);
    timedExecution(C, "UPDATE nodes\n"
                      "SET relation_ids = t.relation_ids\n"
                      "FROM (\n"
                      "         SELECT via_id as node_id, array_agg(relation_id) as relation_ids\n"
                      "         FROM restrictions\n"
                      "         WHERE via_type = 'n'\n"
                      "         group by via_id\n"
                      "     ) t\n"
                      "WHERE nodes.id = t.node_id;", false);
    timedExecution(C, "DROP TABLE IF EXISTS _nodes", false);
}

void processEdges(connection& C, connection& C2, string& fetchUrl, int maxFetch){

    timedExecution(C, "ALTER TABLE edges RENAME COLUMN way_id TO id", true);

    addElevationToEdges(C, C2, fetchUrl, maxFetch);

    timedExecution(C, "ALTER TABLE edges\n"
                      "    ADD COLUMN IF NOT EXISTS u bigint,\n"
                      "    ADD COLUMN IF NOT EXISTS v bigint;", false);

    timedExecution(C, "UPDATE edges\n"
                      "SET u = nodes[ARRAY_LOWER(nodes, 1)];", false);
    timedExecution(C, "UPDATE edges\n"
                      "SET v = nodes[ARRAY_UPPER(nodes, 1)];", false);
    timedExecution(C, "UPDATE edges\n"
                      "SET oneway = 0 WHERE oneway IS NULL;", false);

    timedExecution(C, "ALTER TABLE edges ADD COLUMN IF NOT EXISTS length real;", false);
    timedExecution(C, "UPDATE edges\n"
                      "SET length = ST_Length(geom::geography);", false);

    // relation_ids for edges are already populated through osm2pgsql
//    timedExecution(C, "ALTER TABLE edges ADD COLUMN IF NOT EXISTS relation_ids bigint[];");
//    timedExecution(C, "UPDATE edges\n"
//                      "SET relation_ids = t.relation_ids\n"
//                      "FROM (\n"
//                      "         SELECT via_id as edge_id, array_agg(relation_id) as relation_ids\n"
//                      "         FROM restrictions\n"
//                      "         WHERE via_type = 'w'\n"
//                      "         group by via_id\n"
//                      "     ) t\n"
//                      "WHERE edges.way_id = t.edge_id;");

}


void processRestrictions(connection& C){
    timedExecution(C, "CREATE TYPE restriction_type AS ENUM ('_', 'nn', 'no', 'en', 'eo');", true);
    timedExecution(C, "ALTER TABLE restrictions ADD COLUMN IF NOT EXISTS type_enum restriction_type;", false);
    timedExecution(C, "UPDATE restrictions SET type_enum = 'nn' WHERE via_type = 'n' AND restriction LIKE 'no%';", false);
    timedExecution(C, "UPDATE restrictions SET type_enum = 'no' WHERE via_type = 'n' AND restriction LIKE 'only%'", false);
    timedExecution(C, "UPDATE restrictions SET type_enum = 'en' WHERE via_type = 'e' AND restriction LIKE 'no%'", false);
    timedExecution(C, "UPDATE restrictions SET type_enum = 'eo' WHERE via_type = 'e' AND restriction LIKE 'only%';", false);
}

void addCycling(connection& C){
    timedExecution(C, "CREATE TYPE tag_filtering AS ENUM ('reject', 'accept', 'conditional');", true);
    try {
        timedExecution(C, "ALTER TABLE edges ADD COLUMN cycling tag_filtering;", false);
    } catch (const std::exception &e) {
        timedExecution(C, "UPDATE edges SET cycling = NULL;", false);
    }

    // update whether ways can be used for cycling
    // reference: https://giscience.github.io/openrouteservice/documentation/Tag-Filtering.html
    timedExecution(C, "UPDATE edges SET cycling = 'accept' WHERE \n"
                      "\t( \n"
                      "\t\ttags->>'man_made' = 'pier' \n"
                      "\t \tOR tags->>'railway' = 'platform' \n"
                      "\t\tOR (tags->>'route' in ('shuttle_train', 'ferry') \n"
                      "\t\t\tAND ( tags->>'bicycle' = 'yes' \n"
                      "\t\t\t\t OR (tags->>'bicycle' is not null \n"
                      "\t\t\t\t\t AND tags->>'foot' is not null\n"
                      "\t\t\t\t )\n"
                      "\t\t\t )\n"
                      "\t\t)\n"
                      "\t)\n"
                      "\tAND COALESCE(tags->>'bicycle', '') NOT IN ('private', 'no', 'restricted', 'military', 'emergency')\n"
                      "\tAND COALESCE(tags->>'vehicle', '') NOT IN ('private', 'no', 'restricted', 'military', 'emergency')\n"
                      "\tAND COALESCE(tags->>'access', '') NOT IN ('private', 'no', 'restricted', 'military', 'emergency');", false);
    timedExecution(C, "UPDATE edges SET cycling = 'reject' WHERE\n"
                      "\tCOALESCE(type, '') NOT IN (\n"
                      "\t\t'cycleway', 'path', 'footway', 'pedestrian', 'track', 'service', 'residential', 'living_street', \n"
                      "\t\t'steps', 'unclassified', 'road', 'trunk', 'trunk_link', 'primary', 'primary_link', 'secondary',\n"
                      "\t\t'secondary_link', 'tertiary', 'tertiary_link'\n"
                      "\t);", false);
    timedExecution(C, "UPDATE edges SET cycling = 'accept' WHERE\n"
                      "\ttags->>'sac_scale' = 'hiking' AND type = 'cycleway' AND cycling IS NULL;", false);
    timedExecution(C, "UPDATE edges SET cycling = 'reject' WHERE\n"
                      "\ttags->>'sac_scale' NOT IN (\n"
                      "\t\t'hiking', 'mountain_hiking', 'demanding_mountain_hiking', 'alpine_hiking'\n"
                      "\t) AND cycling IS NULL;", false);
    timedExecution(C, "UPDATE edges SET cycling = 'accept' WHERE\n"
                      "\t(\n"
                      "\t\ttags->>'bicycle' IN (\n"
                      "\t\t\t'yes', 'designated', 'official', 'permissive', 'dismount'\n"
                      "\t\t) OR type = 'cycleway'\n"
                      "\t\tOR tags->>'bicycle_road' = 'yes'\n"
                      "\t) AND cycling IS NULL;", false);
    timedExecution(C, "UPDATE edges SET cycling = 'reject' WHERE\n"
                      "\t(\n"
                      "\t\ttype IN ( 'motorway', 'motorway_link' )\n"
                      "\t\tOR tags->>'motorroad' = 'yes'\n"
                      "\t) AND cycling IS NULL;", false);
    timedExecution(C, "UPDATE edges SET cycling = 'conditional' WHERE\n"
                      "\t(\n"
                      "\t\ttype = 'ford' OR tags->>'ford' IS NOT NULL\n"
                      "\t) AND cycling IS NOT NULL;", false);
    timedExecution(C, "UPDATE edges SET cycling = 'reject' WHERE\n"
                      "\t(\t\n"
                      "\t\tCOALESCE(tags->>'bicycle', '') IN ('private', 'no', 'restricted', 'military', 'emergency')\n"
                      "\t\tOR COALESCE(tags->>'vehicle', '') IN ('private', 'no', 'restricted', 'military', 'emergency')\n"
                      "\t\tOR COALESCE(tags->>'access', '') IN ('private', 'no', 'restricted', 'military', 'emergency')\n"
                      "\t) AND cycling IS NULL;", false);
    timedExecution(C, "UPDATE edges SET cycling = 'accept' WHERE\n"
                      "\tcycling is NULL;", false);

    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS cycling tag_filtering;", true);
    timedExecution(C, "", false);

    for (auto it=tagFiltering.begin(); it != tagFiltering.end(); it++) {
        timedExecution(C, "UPDATE nodes SET cycling = '" + *it + "' WHERE id IN (\n"
                          "\tSELECT nodes[ARRAY_LOWER(nodes, 1)] FROM edges WHERE cycling = '" + *it + "'\n"
                          "\tUNION\n"
                          "\tSELECT nodes[ARRAY_UPPER(nodes, 1)] FROM edges WHERE cycling = '" + *it + "'\n"
                          "\t);", false);
    }
}

void addHiking(connection& C) {
    timedExecution(C, "CREATE TYPE tag_filtering AS ENUM ('reject', 'accept', 'conditional');", true);
    try {
        timedExecution(C, "ALTER TABLE edges ADD COLUMN hiking tag_filtering;", false);
    } catch (const std::exception &e) {
        timedExecution(C, "UPDATE edges SET hiking = NULL;", false);
    }

    // update whether ways can be used for hiking
    // reference: https://giscience.github.io/openrouteservice/documentation/Tag-Filtering.html
    timedExecution(C, "UPDATE edges SET hiking = 'reject' WHERE\n"
                      "\ttags->>'sac_scale' NOT IN (\n"
                      "\t\t'hiking', 'mountain_hiking', 'demanding_mountain_hiking', 'alpine_hiking'\n"
                      "\t);", false);
    timedExecution(C, "UPDATE edges SET hiking = 'accept' WHERE\n"
                      "\ttags->>'foot' IN (\n"
                      "\t\t'yes', 'designated', 'official', 'permissive'\n"
                      "\t) AND hiking IS NULL;", false);
    timedExecution(C, "UPDATE edges SET hiking = 'conditional' WHERE\n"
                      "\ttags->>'foot' IN (\n"
                      "\t\t'private', 'no', 'restricted', 'military', 'emergency'\n"
                      "\t) AND hiking iS NULL;", false);
    timedExecution(C, "UPDATE edges SET hiking = 'accept' WHERE\n"
                      "\ttags->>'sidewalk' IN (\n"
                      "\t\t'yes', 'both', 'left', 'right'\n"
                      "\t) AND hiking IS NULL;", false);
    timedExecution(C, "UPDATE edges SET hiking = 'reject' WHERE\n"
                      "\tCOALESCE(type, '') NOT IN (\n"
                      "\t\t'footway', 'path', 'steps', 'pedestrian', 'living_street', 'track', 'residential', 'service', \n"
                      "\t\t'trunk', 'trunk_link', 'primary', 'primary_link', 'secondary', 'secondary_link', 'tertiary',\n"
                      "\t\t'tertiary_link', 'cycleway', 'unclassified', 'road'\n"
                      "\t) AND hiking IS NULL;", false);
    timedExecution(C, "UPDATE edges SET hiking = 'reject' WHERE\n"
                      "\ttags->>'motorroad' = 'yes' AND hiking IS NULL;", false);
    timedExecution(C, "UPDATE edges SET hiking = 'conditional' WHERE\n"
                      "\t(\n"
                      "\t\ttype = 'ford' OR tags->>'ford' IS NOT NULL\n"
                      "\t) AND hiking IS NULL;", false);

    timedExecution(C, "UPDATE edges SET hiking = 'accept' WHERE hiking IS NULL;", false);
    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS hiking tag_filtering;", true);
    for (auto it=tagFiltering.begin(); it != tagFiltering.end(); it++) {
        timedExecution(C, "UPDATE nodes SET hiking = '" + *it + "' WHERE id IN (\n"
                          "\tSELECT nodes[ARRAY_LOWER(nodes, 1)] FROM edges WHERE hiking = '" + *it + "'\n"
                          "\tUNION\n"
                          "\tSELECT nodes[ARRAY_UPPER(nodes, 1)] FROM edges WHERE hiking = '" + *it + "'\n"
                          "\t);", false);
    }
}


//(man_made = pier OR railway = platform OR ( route = [shuttle_train, ferry] AND (bicycle = yes OR ( bicycle != * AND foot != *)))) AND !(restrictions = restrictedValues)
//restrictions = [bicycle, vehicle, access]
//restrictedValues = [private, no, restricted, military, emergency] + [ dismount, discouraged ]
//intendedValues = [yes, designated, official, permissive]


// highway != [cycleway, path, footway, pedestrian, track, service, residential, living_street, steps, unclassified,
//  road, trunk, trunk_link, primary, primary_link, secondary, secondary_link, tertiary, tertiary_link]

int main(int argc, char *argv[]){
    if (argc < 2){
        cout << "No argument is given for db name." << endl;
        return 1;
    }
    string dbname = argv[1];
    string sqlConnection = "dbname = " + dbname + " user = postgres password = postgres hostaddr = 192.168.1.20 port = 5432";
//    string sqlConnection = "dbname = " + dbname + " user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432";
    string fetchUrl = "http://192.168.1.20:1002/api/v1/lookup";
    int maxFetch = 1000;

    cout << "Database name: " << dbname << endl;
    cout << "URL to fetch elevation data from: " << fetchUrl << endl;

    try
    {
        connection C(sqlConnection);
        connection C2(sqlConnection);

        processNodes(C, C2, fetchUrl, maxFetch);
        processEdges(C, C2, fetchUrl, maxFetch);
        processRestrictions(C);
        addCycling(C);
        addHiking(C);

    }
    catch (const std::exception &e)
    {
        cerr << e.what() << std::endl;
    }

    return 0;
}
