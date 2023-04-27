

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
    // todo: get negative elevation
    int n = 100000;
    int epoch1 = 10;
    int epoch2 = 100;

    pair<Node*, Node*> p = graph->getTwoNearNodes(100000);
    cout << p.first->id << " " << p.second->id << " " << endl;

//    long start_id(10059067058), end_id(6302552417);
    Node *start(graph->getRandomNode()), *end(graph->getRandomNode());
    double maxLengthRatio = 1.5;
    for (int i=0; i<n; i++) {
        cout << endl;

        start = graph->getRandomNode();
        end = graph->getRandomNode();

        cout << "[" << i << "] " << start->id << " " << end->id << endl;
        auto t0 = chrono::high_resolution_clock::now();
        auto t1 = chrono::high_resolution_clock::now();


        cout << endl;
//
        t0 = chrono::high_resolution_clock::now();
        unordered_map<Node *, Node::Edge *> prevEdgeMap;
        int curMinWeight = INT_MAX;
        unordered_map<Node *, int> minWeightEnd = dijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                true, &prevEdgeMap, nullptr
        );
        Path *shortestPath = pathFromPrevEdgeMap(start, end, prevEdgeMap);
        t1 = chrono::high_resolution_clock::now();
        auto t_dl2 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "(shortest Path: length)    length: " << shortestPath->length
             << " elevation: " << shortestPath->elevation
             << " (time: " << t_dl2 << ")" << endl;


        unordered_map<Node*, int> minWeightStart = dijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                false, nullptr, nullptr
        );

        int maxWeight = getMaxWeight(curMinWeight, maxLengthRatio);

        // for reference
        unordered_set<Node*> validNodes;
        for (pair<Node *const, int> &pair1: minWeightEnd) {
            auto it = minWeightStart.find(pair1.first);
            if (it == minWeightStart.end())
                continue;
            if (pair1.second + it->second > maxWeight)
                continue;
            validNodes.insert(pair1.first);
        }
        Path* path = dijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, &validNodes);
        cout << "(possibly optimal)    length: " << path->length
             << " elevation: " << path->elevation << endl;


        for (int i=0; i<1; i++) {
            cout << "--------------------------------------------------" << endl;
            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath4 = geneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStart, minWeightEnd, new PathEdges(shortestPath),
                    100, 25, epoch1, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, true
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_n4 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_nM)               length: " << elenaPath4->length
                 << " elevation: " << elenaPath4->elevation
                 << " (time: " << t_n4 << ")" << endl;

            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath41 = geneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStart, minWeightEnd, new PathEdges(shortestPath),
                    100, 25, epoch2, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, true
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_n41 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_nM)               length: " << elenaPath41->length
                 << " elevation: " << elenaPath41->elevation
                 << " (time: " << t_n41 << ")" << endl;


            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath42 = geneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStart, minWeightEnd, new PathEdges(shortestPath),
                    100, 25, epoch1, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, false
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_n42 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_nm)               length: " << elenaPath42->length
                 << " elevation: " << elenaPath42->elevation
                 << " (time: " << t_n42 << ")" << endl;



            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath43 = geneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStart, minWeightEnd, new PathEdges(shortestPath),
                    100, 25, epoch2, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, false
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_n43 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_nm)               length: " << elenaPath43->length
                 << " elevation: " << elenaPath43->elevation
                 << " (time: " << t_n43 << ")" << endl;


//            //
//            t0 = chrono::high_resolution_clock::now();
//            Path *elenaPath5 = elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
//                    graph, start, end, 2.0, 100, 30, 10, unallowDuplicateUndirectedEdges,
//                    INT_MAX, 10
//            )->toPath();
//            t1 = chrono::high_resolution_clock::now();
//            auto t_e5 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//            cout << "(elena path_eM)               length: " << elenaPath5->length
//                 << " elevation: " << elenaPath5->elevation
//                 << " (time: " << t_e5 << ")" << endl;
        }

        t0 = chrono::high_resolution_clock::now();
        Node::Edge* lastEdge = nullptr;
        curMinWeight = INT_MAX;
        unordered_map<Node::Edge *, Node::Edge *> prevEdgeMapE;
        unordered_map<Node::Edge *, int> minWeightEndE = edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                true, &prevEdgeMapE, lastEdge, nullptr
        );
        Path* shortestPathE = edgeBasedPathFromPrevEdgeMap(lastEdge, prevEdgeMapE);
        t1 = chrono::high_resolution_clock::now();
        auto t_dle2 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "(shortest path (edge): length)    length: " << shortestPathE->length
             << " elevation: " << shortestPathE->elevation
             << " (time: " << t_dle2 << ")" << endl;


        unordered_map<Node::Edge*, int> minWeightStartE = edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                false, nullptr, lastEdge, nullptr
        );


        for (int i=0; i<1; i++) {
            cout << "--------------------------------------------------" << endl;
            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath4 = edgeBasedGeneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStartE, minWeightEndE, new PathEdges(shortestPath),
                    100, 25, epoch1, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, true
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_e4 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_eM)               length: " << elenaPath4->length
                 << " elevation: " << elenaPath4->elevation
                 << " (time: " << t_e4 << ")" << endl;

            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath41 = edgeBasedGeneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStartE, minWeightEndE, new PathEdges(shortestPath),
                    100, 25, epoch2, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, true
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_e41 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_eM)               length: " << elenaPath41->length
                 << " elevation: " << elenaPath41->elevation
                 << " (time: " << t_e41 << ")" << endl;

            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath42 = edgeBasedGeneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStartE, minWeightEndE, new PathEdges(shortestPath),
                    100, 25, epoch1, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, false
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_e42 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_em)               length: " << elenaPath42->length
                 << " elevation: " << elenaPath42->elevation
                 << " (time: " << t_e42 << ")" << endl;

            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath43 = edgeBasedGeneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStartE, minWeightEndE, new PathEdges(shortestPath),
                    100, 25, epoch2, unallowDuplicateUndirectedEdges,
                    INT_MAX, 1000, false
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_e43 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_em)               length: " << elenaPath43->length
                 << " elevation: " << elenaPath43->elevation
                 << " (time: " << t_e43 << ")" << endl;


//            //
//            t0 = chrono::high_resolution_clock::now();
//            Path *elenaPath5 = elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
//                    graph, start, end, 2.0, 100, 30, 10, unallowDuplicateUndirectedEdges,
//                    INT_MAX, 10
//            )->toPath();
//            t1 = chrono::high_resolution_clock::now();
//            auto t_e5 = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//            cout << "(elena path_eM)               length: " << elenaPath5->length
//                 << " elevation: " << elenaPath5->elevation
//                 << " (time: " << t_e5 << ")" << endl;
        }
//        Path* path = elenaPath4;
//        unordered_set<Node::Edge*> uset(path->edges.begin(), path->edges.end());
//        unordered_set<long> uset2;
//        for (Node::Edge* edge: path->edges){
//            uset2.insert(edge->id);
//        }
//        cout << path->edges.size() << " " << uset2.size() << " " << uset.size() << endl;
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

void runUsingDBMap(string connectionString, HighwayConfig highwayConfig, LocationConfig locationConfig){

    Graph* graph = new Graph(true);
    pqxx::connection C(connectionString);
    graph->addAllNodesFromDB(C, highwayConfig, locationConfig, false);
    graph->addAllEdgesFromDB(C, highwayConfig, locationConfig, false);
    graph->maxGroup();
    graph->addAllRestrictionsFromDB(C, highwayConfig, locationConfig);
    graph->createBallTree();
    runRandom(graph);
}

int main(int argc, char *argv[]){
//    runUsingLocalMap("maps/test", false);
    if (argc < 3){
        cout << "Argument missing." << endl;
        return 1;
    }
    int port = stoi(argv[1]);
    string dbname = argv[2];

    string connectionString = "dbname = " + dbname + " user = postgres password = postgres "
                                                     "hostaddr = 192.168.1.20 port = 5432";
//
    cout << connectionString << endl;
    runServer(connectionString, port);

    if (argc < 4){
        cout << "Argument missing." << endl;
        return 1;
    }

    string configName = argv[3];
    HighwayConfig mapConfig;
    if (configName == "all") {
        mapConfig = all_highways;
    } else if (configName == "motorway"){
        mapConfig = motorway;
    } else if (configName == "trunk"){
        mapConfig = trunk;
    } else if (configName == "primary"){
        mapConfig = primary;
    } else if (configName == "secondary"){
        mapConfig = secondary;
    } else if (configName == "tertiary"){
        mapConfig = tertiary;
    } else if (configName == "local") {
        mapConfig = local;
    } else if (configName == "cycling") {
        mapConfig = cycling;
    } else if (configName == "hiking") {
        mapConfig = hiking;
    } else {
        cout << "invalid configuration name: " << configName << endl;
        return 1;
    }

    runUsingDBMap(connectionString, mapConfig, US_mainland);
    return 0;
}
