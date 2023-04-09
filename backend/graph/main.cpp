

#include <iostream>
#include <vector>
#include <chrono>
#include <pqxx/pqxx>
#include <float.h>

#include "include/restApi.h"

#ifndef PATHFINDING_H
#define PATHFINDING_H
#include "include/pathfinding.h"
#endif

#ifndef QUERY_H
#define QUERY_H
#include "include/query.h"
#endif

//#include "include/restApi.h"

using namespace std;
using namespace pqxx;

void runRandom(Graph* graph){
    // todo: getnegative elevation
    cout << graph->maxElevation << endl;
    int n = 100000;

    long start_id(10059067058), end_id(6302552417);
    for (int i=0; i<n; i++) {
//        start_id = graph->getRandomNodeId();
//        end_id = graph->getRandomNodeId();
        cout << endl;
        cout << "[" << i << "] " << start_id << " " << end_id << endl;

        Node *start(graph->nodes[start_id]), *end(graph->nodes[end_id]);

        auto t0 = chrono::high_resolution_clock::now();
//        Edges *aStarPath = aStarAlgorithm(graph, start, end, Node::Edge::getLength);
        auto t1 = chrono::high_resolution_clock::now();
//        auto t_sl = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//        cout << "(a-star path: length)      length: " << aStarPath->length
//             << " elevation: " << aStarPath->elevation
//             << " (time: " << t_sl << ")" << endl;

//        t0 = chrono::high_resolution_clock::now();
//        Path *dijkstraPath = edgeBasedDijkstraAlgorithm(graph, start, end, Node::Edge::getLength, false, nullptr);
//        t1 = chrono::high_resolution_clock::now();
//        auto t_dl = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//        cout << "(dijkstra path: length)    length: " << dijkstraPath->length
//             << " elevation: " << dijkstraPath->elevation
//             << " (time: " << t_dl << ")" << endl;

//        t0 = chrono::high_resolution_clock::now();
//        Edges *aStarPathE = aStarAlgorithm(graph, start, end, Node::Edge::getElevation);
//        t1 = chrono::high_resolution_clock::now();
//        auto t_se = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//        cout << "(a-star path: elevation)   length: " << aStarPathE->length
//             << " elevation: " << aStarPathE->elevation
//             << " (time: " << t_se << ")" << endl;
//
        t0 = chrono::high_resolution_clock::now();
        Path *dijkstraPathE = edgeBasedDijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, false, nullptr);
        t1 = chrono::high_resolution_clock::now();
        auto t_de = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "(dijkstra path: elevation) length: " << dijkstraPathE->length
             << " elevation: " << dijkstraPathE->elevation
             << " (time: " << t_de << ")" << endl;

        t0 = chrono::high_resolution_clock::now();
        Path *elenaPath = elenaPathFindMinUsingDijkstra(graph, start, end, 1.5);
        t1 = chrono::high_resolution_clock::now();
        auto t_e1 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "(elena path)               length: " << elenaPath->length
             << " elevation: " << elenaPath->elevation
             << " (time: " << t_e1 << ")" << endl;

//        t0 = chrono::high_resolution_clock::now();
//        t1 = chrono::high_resolution_clock::now();
//        auto t_e2 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//        cout << "(elena path)               length: " << elenaPath2->length
//             << " elevation: " << elenaPath2->elevation
//             << " (time: " << t_e2 << ")" << endl;

    }
}

//void runUsingLocalMap(string map_name, bool preprocessed){
//    //
//
////    auto* graph = new Graph("maps/test", false);
////    auto* graph = new Graph("maps/Helsinki", false);
////    auto* graph = new Graph("maps/test.maxGroup", true);
////    auto* graph = new Graph("maps/Helsinki.maxGroup", true);
//
//    auto* graph = new Graph(map_name, preprocessed, true);
//    runRandom(graph);
//}

void runUsingDBMap(string connectionString){

    Graph* graph = new Graph(true);
//    try {
    pqxx::connection C(connectionString);
//        pqxx::connection C("dbname = osm-mass user = postgres password = postgres \
//             hostaddr = 192.168.1.20 port = 5432");
    prepare_nodes_query(C, all);
    prepare_edges_query(C, all);
    prepare_restrictions_query(C, all);
    graph->addAllNodesFromDB(C, false);
    graph->addAllEdgesFromDB(C, false);
    graph->addAllRestrictionsFromDB(C);
//    graph->maxGroup();
//    graph->edgeBasedMaxGroup();
    graph->maxGroup();
//    } catch  (const std::exception &e){
//        cerr << e.what() << std::endl;
//        throw e;
//    }

    runRandom(graph);



    // experiments for genetic algorithm
//    auto t0 = chrono::high_resolution_clock::now();
//    int n = 20;
//    int elevation(0);
//    int maxElevation(0);
//    for (int i=0; i<n; i++) {
//        Node *start = graph->getRandomNode();
//        Node *end = graph->getRandomNode();
//        start = graph->nodes[4405347363];
//        end = graph->nodes[3616168673];
////        cout << start->id << " " << end->id << endl;
//        PathEdges* pathEdges = geneticAlgorithm(graph, start, end, 1.5, 50, 20, 1000);
//        cout << pathEdges->getElevation() << endl;
//        elevation += pathEdges->getElevation();
//        maxElevation = max(maxElevation, pathEdges->getElevation());
//    }
//    auto t1 = chrono::high_resolution_clock::now();
//    auto t_e1 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//    cout << t_e1 << " seconds" << endl;
//    cout << elevation / n << " " << maxElevation << endl;
    // front / back (165.504 seconds)
    // 26087 32176
    // 25425 30677 165.504 seconds
    // random
    // 24359 32252
    // 24350 32189 173.75 seconds
    // square weighted
    // 24875 32281 185.347 seconds
}

int main(int argc, char *argv[]){
//    runUsingLocalMap("maps/test", false);
    if (argc < 3){
        cout << "No argument is given for db name." << endl;
        return 1;
    }
    int port = stoi(argv[1]);
    string dbname = argv[2];

    string connectionString = "dbname = " + dbname + " user = postgres password = postgres \
             hostaddr = 192.168.1.20 port = 5432";


//    runUsingDBMap(connectionString);
    runServer(connectionString, port, all);

//    Graph* graph = new Graph(connectionString);
//    cout << graph->nearestNode({-71.0505, 42.1017})->id << endl;








//    Graph* graph = new Graph(true);
//    int numNodes = 10;
//    vector<Node*> v;
//    for (long i=0; i<numNodes; i++){
//        v.push_back(graph->addNode(i, 0.0, 0.0, 0.0));
//    }
//    cout << endl;
//
//    graph->addEdge(0, v[0], v[1], 100, 100, oneway);
//    graph->addEdge(1, v[1], v[2], 50, 100, oneway);
//    graph->addEdge(2, v[2], v[3], 50, 100, oneway);
//    graph->addEdge(3, v[0], v[2], 50, 100, oneway);
//    graph->addEdge(4, v[0], v[4], 50, 100, oneway);
//    graph->addEdge(5, v[4], v[2], 50, 100, oneway);
//    Node *start(v[0]), *end(v[3]);
//
//    double curMinWeight = DBL_MAX;
//    Node::Edge* edge = nullptr;
//    double maxWeightRatio = 1.5;
//    vector<Node::Edge*> edges;
//    unordered_map<Node::Edge*, Node::Edge*> prevEdgeMap;
//
//    bool reversed = false;
//    unordered_map<Node::Edge*, double> minWeightStart = edgeBasedDijkstraAlgorithm(
//            graph, start, end, curMinWeight, maxWeightRatio, Node::Edge::getLength, reversed, &prevEdgeMap, edge, nullptr
//    );
////    cout << edge->id << endl;
//
//    if (edge != nullptr){
//        cout << "found " << "minWeight: " << curMinWeight << endl;
//        double length = 0.0;
//        int elevation = 0;
//        while (edge->id != -1){
////            edges.push_back(edge);
//            cout << edge->id << " ";
//            length += edge->getLength();
//            elevation += edge->getElevation();
//            edge = prevEdgeMap[edge];
//        }
//        cout << endl;
//    }
//
//    prevEdgeMap.clear();
//    edge = nullptr;
//
//    reversed = true;
//    unordered_map<Node::Edge*, double> minWeightEnd = edgeBasedDijkstraAlgorithm(
//            graph, start, end, curMinWeight, maxWeightRatio, Node::Edge::getLength, reversed, &prevEdgeMap, edge, nullptr
//    );
//
//    if (edge != nullptr){
//        cout << "found " << "minWeight: " << curMinWeight << endl;
//        double length = 0.0;
//        int elevation = 0;
//        while (edge->id != -1){
////            edges.push_back(edge);
//            cout << edge->id << " ";
//            length += edge->getLength();
//            elevation += edge->getElevation();
//            edge = prevEdgeMap[edge];
//        }
//        cout << endl;
//    }


    return 0;
}
