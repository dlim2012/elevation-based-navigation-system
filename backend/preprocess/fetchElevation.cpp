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

void addElevationToNodes(
        connection& C,
        connection& C2,
        string& fetchUrl,
        int maxFetch
){

    long numRows = getNumberOfRows(C, "nodes");
    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS elevation INT;");
    timedExecution(C, "DROP TABLE IF EXISTS _nodes");
    timedExecution(C, "CREATE TABLE IF NOT EXISTS _nodes ("
                      "    id BIGINT PRIMARY KEY,"
                      "    elevation INT"
                      ")");

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
            countWrote += _processNodes(node_ids, points, fetchUrl, stream2);
            node_ids.clear();
            points.clear();
            if (countWrote * 100 > percentage * numRows){
                while (countWrote * 100 > percentage * numRows)
                    percentage++;
                auto t1 = chrono::high_resolution_clock::now();
                auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << setw(8) << right << countWrote << " rows / " << setw(8) << right << numRows << " rows (time: " << right << setw(8) << t << " seconds, expected time left: "
                     << right << setw(8) << t / countWrote * (numRows - countWrote) << +" seconds)" << endl;

            }
        }

    }
    if (countWrote < numRows){
        countWrote += _processNodes(node_ids, points, fetchUrl, stream2);
        if (countWrote * 100 > percentage * numRows){
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

    timedExecution(C, "UPDATE nodes SET elevation = _nodes.elevation FROM _nodes WHERE nodes.id = _nodes.id");
    timedExecution(C, "DROP TABLE _nodes");
}





void addElevationToEdges(
        connection& C,
        connection& C2,
        string& fetchUrl,
        int maxFetch
){

    long numRows = getNumberOfRows(C, "edges");
    long numPoints = getNumberOfPoints(C, "edges");
    timedExecution(C, "ALTER TABLE edges ADD COLUMN IF NOT EXISTS elevation INT;");
    timedExecution(C, "DROP TABLE IF EXISTS _edges");
    timedExecution(C, "CREATE TABLE IF NOT EXISTS _edges (\n"
                      "    id BIGINT PRIMARY KEY,\n"
                      "    elevation INT\n"
                      ");");

    cout << "Fetching elevation data for edges and saving to the new table... (total " << numRows << " rows, " << numPoints << " points)" << endl;
    work tx1(C);
    work tx2(C2);


    long countEdgeRead(0), countEdgeWrote(0);
    long countPointRead(0), countPointProcessed(0), countPointCur(0);
    long percentage(0);
    vector<long> rows;
    vector<int> edgeGeomSizes;
    vector<pair<string, string>> points;
    auto t0 = chrono::high_resolution_clock::now();


    result r = tx1.exec("SELECT way_id, ST_AsText(geom) FROM edges;");
    tx1.commit();
//    stream_from stream1(tx1, "_edges", vector<string>{"id", "geom_text"});
    long row1;

    stream_to stream2(tx2, "_edges", vector<string>{"id", "elevation"});

    for (row row: r){
//    while (stream1 >> row1){
//        countEdgeRead++;
        row1 = stol(row[0].c_str());

        string geom_text = row[1].c_str();
        vector<pair<string, string>> linePoints = parseLineString(geom_text);
        if (linePoints.size() < 2){
            cout << "Encountered an edge with less than two points. Ignoring this edge." << endl;
            continue;
        }

        if (countPointCur + linePoints.size() > maxFetch){
            pair<int, int> counts = _processEdges(rows, edgeGeomSizes, points, fetchUrl, stream2);
            countEdgeWrote += counts.first;
            countPointProcessed += counts.second;
            rows.clear();
            edgeGeomSizes.clear();
            points.clear();
            countPointCur = 0;

            if (countEdgeWrote * 100 > percentage * numRows){
                while (countEdgeWrote * 100 > percentage * numRows)
                    percentage++;
                auto t1 = chrono::high_resolution_clock::now();
                auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << right << setw(8) << countEdgeWrote << " rows / " << setw(8) <<  right << numRows << " rows, "
                        << setw(8)<<  right << countPointProcessed << " points / " << setw(8) <<  right << numPoints << " points (time: " << right
                        << setw(8) << t << " seconds, expected time left: " << right
                     << setw(8) << t / countPointProcessed * (numPoints - countPointProcessed) << +" seconds)" << endl;

            }
        }

        rows.push_back(row1);
        edgeGeomSizes.push_back(linePoints.size());
        points.insert(points.end(), make_move_iterator(linePoints.begin()),make_move_iterator(linePoints.end()));
        countPointRead += linePoints.size();
        countPointCur += linePoints.size();

    }
    pair<int, int> counts = _processEdges(rows, edgeGeomSizes, points, fetchUrl, stream2);
    countEdgeWrote += counts.first;
    countPointProcessed += counts.second;
    if (countEdgeWrote * 100 > percentage * numRows){
        while (countEdgeWrote * 100 > percentage * numRows)
            percentage++;
        auto t1 = chrono::high_resolution_clock::now();
        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "[" << setw(3) << percentage << "%] " << fixed << setprecision(2) << right << setw(8) << countEdgeWrote << " rows / " << setw(8) <<  right << numRows << " rows, "
             << setw(8)<<  right << countPointProcessed << " points / " << setw(8) <<  right << numPoints << " points (time: " << right
             << setw(8) << t << " seconds, expected time left: " << right
             << setw(8) << t / countPointProcessed * (numPoints - countPointProcessed) << +" seconds)" << endl;

    }
//    stream1.complete();
    stream2.complete();
    tx2.commit();

    cout << countEdgeRead << " " << countEdgeWrote << " " << countPointRead << " " << countPointProcessed << endl;

    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
    cout << "(time: " << t << " seconds)" << endl << endl;

    timedExecution(C, "UPDATE edges SET elevation = _edges.elevation FROM _edges WHERE edges.way_id = _edges.id");
    timedExecution(C, "DROP TABLE _edges;");
}

void processNodes(connection& C, connection& C2, string& fetchUrl, int maxFetch){
    // Add geom_text temporarily to use the stream functionality of libpqxx

    timedExecution(C, "CREATE TABLE IF NOT EXISTS nodes (\n"
                      "     id bigint primary KEY,\n"
                      "     geom geometry(POINT, 4326)\n"
                      ");");
    timedExecution(C, "INSERT INTO nodes(id, geom)\n"
                      "    SELECT nodes[array_lower(nodes, 1)], ST_PointN(geom, 1) FROM edges\n"
                      "    ON CONFLICT DO NOTHING;");
    timedExecution(C, "INSERT INTO nodes(id, geom)\n"
                      "    SELECT nodes[array_upper(nodes, 1)], ST_PointN(geom, -1) FROM edges\n"
                      "    ON CONFLICT DO NOTHING;");

    // Fetch elevation from Open-Elevation API and save to the new table;
    addElevationToNodes(C, C2, fetchUrl, maxFetch);

    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS longitude real;");
    timedExecution(C, "UPDATE nodes SET longitude = ST_X(geom);");

    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS latitude real;");
    timedExecution(C, "UPDATE nodes SET latitude = ST_Y(geom);");

    timedExecution(C, "ALTER TABLE nodes ADD COLUMN IF NOT EXISTS relation_ids bigint[];");
    timedExecution(C, "UPDATE nodes\n"
                      "SET relation_ids = t.relation_ids\n"
                      "FROM (\n"
                      "         SELECT via_id as node_id, array_agg(relation_id) as relation_ids\n"
                      "         FROM restrictions\n"
                      "         WHERE via_type = 'n'\n"
                      "         group by via_id\n"
                      "     ) t\n"
                      "WHERE nodes.id = t.node_id;");
    timedExecution(C, "DROP TABLE IF EXISTS _nodes");
}

void processEdges(connection& C, connection& C2, string& fetchUrl, int maxFetch){

    addElevationToEdges(C, C2, fetchUrl, maxFetch);

    timedExecution(C, "ALTER TABLE edges\n"
                      "    ADD COLUMN IF NOT EXISTS u bigint,\n"
                      "    ADD COLUMN IF NOT EXISTS v bigint;");

    timedExecution(C, "UPDATE edges\n"
                      "SET u = nodes[array_lower(nodes, 1)];");
    timedExecution(C, "UPDATE edges\n"
                      "SET v = nodes[array_upper(nodes, 1)];");
    timedExecution(C, "UPDATE edges\n"
                      "SET oneway = 0 WHERE oneway IS NULL;");

    timedExecution(C, "ALTER TABLE edges ADD COLUMN IF NOT EXISTS length real;");
    timedExecution(C, "UPDATE edges\n"
                      "SET length = ST_Length(geom::geography);");

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

    timedExecution(C, "ALTER TABLE edges RENAME COLUMN way_id TO id");

}


void processRestrictions(connection& C){
    string createTypeEnumQuery = "CREATE TYPE restriction_type AS ENUM ('_', 'nn', 'no', 'en', 'eo');";
    try {
        timedExecution(C, createTypeEnumQuery);
    } catch (const std::exception &e) {
        cout << "\"" << createTypeEnumQuery << "\" failed. (skipping)" << endl;
    }
    timedExecution(C, "ALTER TABLE restrictions ADD COLUMN IF NOT EXISTS type_enum restriction_type;");
    timedExecution(C, "UPDATE restrictions SET type_enum = 'nn' WHERE via_type = 'n' AND restriction LIKE 'no%';");
    timedExecution(C, "UPDATE restrictions SET type_enum = 'no' WHERE via_type = 'n' AND restriction LIKE 'only%'");
    timedExecution(C, "UPDATE restrictions SET type_enum = 'en' WHERE via_type = 'e' AND restriction LIKE 'no%'");
    timedExecution(C, "UPDATE restrictions SET type_enum = 'eo' WHERE via_type = 'e' AND restriction LIKE 'only%';");
}

int main(int argc, char *argv[]){
    if (argc < 2){
        cout << "No argument is given for db name." << endl;
        return 1;
    }
    string dbname = argv[1];
    string sqlConnection = "dbname = " + dbname + " user = postgres password = postgres hostaddr = 192.168.1.20 port = 5432";
//    string sqlConnection = "dbname = " + dbname + " user = postgres password = postgres hostaddr = 127.0.0.1 port = 5432";
    string fetchUrl = "http://192.168.1.20:1000/api/v1/lookup";
    int maxFetch = 1000;

    try
    {
        connection C(sqlConnection);
        connection C2(sqlConnection);

        processNodes(C, C2, fetchUrl, maxFetch);
        processEdges(C, C2, fetchUrl, maxFetch);
        processRestrictions(C);

    }
    catch (const std::exception &e)
    {
        cerr << e.what() << std::endl;
    }

    return 0;
}
