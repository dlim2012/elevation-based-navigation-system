

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

#include "include/restApi.h"

using namespace std;
using namespace pqxx;

void testForbidDuplicateEdges(Path* path, DuplicateEdge duplicateEdge, bool edgeBased){
    unordered_set<Node::Edge*> uset(path->getEdges().begin(), path->getEdges().end());
    unordered_set<long> uset2;
    for (Node::Edge* edge: path->getEdges()){
        uset2.insert(edge->getId());
    }
    cout << path->getEdges().size() << " " << uset.size() << " " << uset2.size() << endl;
    if (edgeBased){
        return;
    }
    if (duplicateEdge == minimizeDuplicateDirectedEdges){
        assert (path->getEdges().size() == uset.size());
    } else if (duplicateEdge == minimizeDuplicateUndirectedEdges){
        assert (path->getEdges().size() == uset.size());
        assert (path->getEdges().size() == uset2.size());
    }
}

void testRandomRun(Graph* graph, DuplicateEdge duplicateEdge, int n, vector<int> num_epochs){


//    long start_id(10059067058), end_id(6302552417);
    Node *start(graph->getRandomNode()), *end(graph->getRandomNode());
    double maxLengthRatio = 1.5;
    for (int i=0; i<n; i++) {
        cout << endl;

        start = graph->getRandomNode();
        end = graph->getRandomNode();



        cout << "[" << i << "] " << start->getId() << " " << end->getId() << endl;
        auto t0 = chrono::high_resolution_clock::now();
        auto t1 = chrono::high_resolution_clock::now();

        int curMinWeight = INT_MAX;
        int maxWeight;

        t0 = chrono::high_resolution_clock::now();
        unordered_map<Node *, Node::Edge *> prevEdgeMap;
        unordered_map<Node *, int> minWeightEnd = dijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                true, &prevEdgeMap, nullptr
        );
        Path *shortestPath = pathFromPrevEdgeMap(start, end, prevEdgeMap);
        t1 = chrono::high_resolution_clock::now();
        auto t_s = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "(shortest Path: length)    length: " << shortestPath->getLength()
             << " elevation: " << shortestPath->getElevation()
             << " (time: " << t_s << ")" << endl;
        assert (start == end ? shortestPath->getLength() == 0: shortestPath->getLength() > 0);

        unordered_map<Node*, int> minWeightStart = dijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                false, nullptr, nullptr
        );

        maxWeight = getMaxWeight(curMinWeight, maxLengthRatio);

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
        cout << "(possibly optimal)    length: " << path->getLength()
             << " elevation: " << path->getElevation() << endl;
        delete path;


        for (int num_epoch: num_epochs){
            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath_nM = geneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStart, minWeightEnd, shortestPath,
                    100, 25, num_epoch, duplicateEdge,
                    INT_MAX, 1000, true
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_nM = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_nM)               length: " << elenaPath_nM->getLength()
                 << " elevation: " << elenaPath_nM->getElevation()
                 << " (time: " << t_nM << ")" << endl;


            assert (elenaPath_nM->getLength() >= shortestPath->getLength());
            assert (elenaPath_nM->getLength() < maxWeight);
            assert (elenaPath_nM->getElevation() >= shortestPath->getElevation());
            testForbidDuplicateEdges(elenaPath_nM, duplicateEdge, false);
            delete elenaPath_nM;

            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath_nm = geneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStart, minWeightEnd, shortestPath,
                    100, 25, num_epoch, duplicateEdge,
                    INT_MAX, 1000, false
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_nm = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_nm)               length: " << elenaPath_nm->getLength()
                 << " elevation: " << elenaPath_nm->getElevation()
                 << " (time: " << t_nm << ")" << endl;

            assert (elenaPath_nm->getLength() >= shortestPath->getLength());
            assert (elenaPath_nm->getLength() < maxWeight);
            assert (elenaPath_nm->getElevation() <= shortestPath->getElevation());
            assert (elenaPath_nm->getElevation() >= 0);
            testForbidDuplicateEdges(elenaPath_nm, duplicateEdge, false);
            delete elenaPath_nm;
        }
        delete shortestPath;

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
        auto t_sE = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "(shortest path (edge): length)    length: " << shortestPathE->getLength()
             << " elevation: " << shortestPathE->getElevation()
             << " (time: " << t_sE << ")" << endl;
        assert (start == end ? shortestPathE->getLength() == 0: shortestPathE->getLength() > 0);

        maxWeight = getMaxWeight(curMinWeight, maxLengthRatio);

        unordered_map<Node::Edge*, int> minWeightStartE = edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
                false, nullptr, lastEdge, nullptr
        );

        for (int num_epoch: num_epochs){
            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath_eM = edgeBasedGeneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStartE, minWeightEndE, shortestPathE,
                    100, 25, num_epoch, duplicateEdge,
                    INT_MAX, 1000, true
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_eM = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_eM)               length: " << elenaPath_eM->getLength()
                 << " elevation: " << elenaPath_eM->getElevation()
                 << " (time: " << t_eM << ")" << endl;

            assert (elenaPath_eM->getLength() >= shortestPathE->getLength());
            assert (elenaPath_eM->getLength() < maxWeight);
            assert (elenaPath_eM->getElevation() >= shortestPathE->getElevation());
            testForbidDuplicateEdges(elenaPath_eM, duplicateEdge, true);
            delete elenaPath_eM;

            t0 = chrono::high_resolution_clock::now();
            Path *elenaPath_em = edgeBasedGeneticAlgorithm(
                    graph, start, end, maxWeight, minWeightStartE, minWeightEndE, shortestPathE,
                    100, 25, num_epoch, duplicateEdge,
                    INT_MAX, 1000, false
            );
            t1 = chrono::high_resolution_clock::now();
            auto t_em = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
            cout << "(elena path_em)               length: " << elenaPath_em->getLength()
                 << " elevation: " << elenaPath_em->getElevation()
                 << " (time: " << t_em << ")" << endl;

            assert (elenaPath_em->getLength() >= shortestPathE->getLength());
            assert (elenaPath_em->getLength() < maxWeight);
            assert (elenaPath_em->getElevation() <= shortestPathE->getElevation());
            assert (elenaPath_em->getElevation() >= 0);
            testForbidDuplicateEdges(elenaPath_em, duplicateEdge, true);
            delete elenaPath_em;

        }
        delete shortestPathE;
    }
}


void test(string connectionString, HighwayConfig highwayConfig, LocationConfig locationConfig){

    Graph* graph = new Graph(true);
    pqxx::connection C(connectionString);
    graph->addAllNodesFromDB(C, highwayConfig, locationConfig, false);
    graph->addAllEdgesFromDB(C, highwayConfig, locationConfig, false);
    graph->maxGroup();
    graph->addAllRestrictionsFromDB(C, highwayConfig, locationConfig);
    graph->createBallTree();

    pair<Node*, Node*> p = graph->getTwoNearNodes(100000);

    int n = 1;
    vector<int> num_epochs = {10};
    testRandomRun(graph, minimizeDuplicateUndirectedEdges, n, num_epochs);
    testRandomRun(graph, minimizeDuplicateDirectedEdges, n, num_epochs);
    testRandomRun(graph, ignoreDuplicateEdges, n, num_epochs);

    delete graph;
}

int main(int argc, char *argv[]){
    if (argc < 3){
        cout << "Argument missing." << endl;
        return 1;
    }
    int port = stoi(argv[1]);
    string dbname = argv[2];

    string connectionString = "dbname = " + dbname + " user = postgres password = postgres "
                                                     "hostaddr = 192.168.1.20 port = 5432";

    cout << connectionString << endl;
    Graph* graph;
    if (port >= 0){
        runServer(connectionString, port);
    } else {
        if (argc < 4){
            cout << "Argument missing." << endl;
            return 1;
        }

        string configName = argv[3];
        HighwayConfig highwayConfig;
        LocationConfig locationConfig = US_mainland;
        if (configName == "all") {
            highwayConfig = all_highways;
        } else if (configName == "motorway"){
            highwayConfig = motorway;
        } else if (configName == "trunk"){
            highwayConfig = trunk;
        } else if (configName == "primary"){
            highwayConfig = primary;
        } else if (configName == "secondary"){
            highwayConfig = secondary;
        } else if (configName == "tertiary"){
            highwayConfig = tertiary;
        } else if (configName == "local") {
            highwayConfig = local;
        } else if (configName == "cycling") {
            highwayConfig = cycling;
        } else if (configName == "hiking") {
            highwayConfig = hiking;
        } else {
            cout << "invalid configuration name: " << configName << endl;
            return 1;
        }
//        graph = new Graph(true);
//        pqxx::connection C(connectionString);
//        graph->addAllNodesFromDB(C, highwayConfig, locationConfig, false);
//        graph->addAllEdgesFromDB(C, highwayConfig, locationConfig, false);
//        graph->maxGroup();
//        graph->addAllRestrictionsFromDB(C, highwayConfig, locationConfig);
//        graph->createBallTree();
        test(connectionString, highwayConfig, US_mainland);

    }
    return 0;
}
