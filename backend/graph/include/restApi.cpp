
#include <iostream>
#include <chrono>
#include <deque>

#include "crow.h"
#include "restApi.h"
#include "crow/middlewares/cors.h"
#include "RequestRateLimiter.h"

using namespace crow;
using namespace std;

typedef std::chrono::time_point<std::chrono::system_clock> tp;

const int LOCK_COUNT = 20;

unordered_map<long, Geometry*>* sharedEdgesGeometries;

Graph *main_hiking_graph;
Graph *main_cycling_graph;
Graph *main_motorway_graph;

RequestRateLimiter* requestRateLimiter;

Graph* graphTypeToGraph(string graphType){
    // todo: add Alaska and islands and find the graph based on the nearest node from the start location
    if (graphType == "main_motorway"){
        return main_motorway_graph;
    } else if (graphType == "hiking") {
        return main_hiking_graph;
    } else if (graphType == "cycling"){
        return main_cycling_graph;
    }
    throw runtime_error("Invalid graph type");
}

const int MAXIMUM_MAX_LENGTH_RATIO = 10;

crow::json::wvalue pathToGeoJson(
        Path *path, unordered_map<long, Geometry*> *umap) {
    vector<Node::Edge *> &edges = path->edges;
    crow::json::wvalue w;
    w["GeoJSON"]["type"] = "LineString";

    if (edges.empty()) {
        return w;
    }

    int index = 0;
    Node* start = path->getStart();
    long nodeId = start->id;
    w["GeoJSON"]["coordinates"][index][0] = start->lat;
    w["GeoJSON"]["coordinates"][index][1] = start->lon;
    index++;
    int count = 0;
    for (Node::Edge *edge: edges) {
        Geometry* geom = (*umap)[edge->id];
        if (geom->startId == nodeId){
            for (int i = 1; i < geom->points.size(); i++) {
                w["GeoJSON"]["coordinates"][index][0] = geom->points[i].second;
                w["GeoJSON"]["coordinates"][index][1] = geom->points[i].first;
                index++;
            }
            nodeId = geom->endId;
        } else if (geom->endId == nodeId){
            for (int i = geom->points.size()-2; i >=0; i--) {
                w["GeoJSON"]["coordinates"][index][0] = geom->points[i].second;
                w["GeoJSON"]["coordinates"][index][1] = geom->points[i].first;
                index++;
            }
            nodeId = geom->startId;
        } else {
            throw runtime_error("");
        }
    }
    // sum of elevation diff to elevation gain
    w["elevation"] = (path->elevation - path->getStart()->elevation + path->getEnd()->elevation) / 2;
    w["length"] = (double) path->length / 100000;  // cm to km
    w["result"] = "o";
    return w;
}


crow::json::wvalue jsonFromNode(Node* node){

    json::wvalue w;
    w["id"] = node->id;
    w["lon"] = node->lon;
    w["lat"] = node->lat;
    w["elevation"] = node->elevation;

    if (node->restrictions != nullptr) {
        for (pair<const Node::Edge *, unordered_set<Node::Edge *>> pair1: *(node->restrictions)) {
            int index = 0;
            for (auto it: pair1.second){
                w["restrictions"][to_string(pair1.first->id)][index++] = it->id;
            }
        }
    }
    return w;
}


void runServer(string& connectionString, int port) {
    requestRateLimiter = new RequestRateLimiter(LOCK_COUNT);

    sharedEdgesGeometries = new unordered_map<long, Geometry*>();
//    cout << "Simple driveway" << endl;
//    main_motorway_graph = new Graph(connectionString, motorway);
    cout << "Hiking" << endl;
    main_hiking_graph = new Graph(connectionString, hiking, sharedEdgesGeometries);
//    cout << endl << "Cycling" << endl;
//    main_cycling_graph = new Graph(connectionString, cycling_with_restrictions);
    cout << endl << "Cycling with restrictions" << endl;
    main_cycling_graph = new Graph(connectionString, cycling_with_restrictions, sharedEdgesGeometries);

    unordered_map<string, double> recentUsagePerIPAddress;
    double recentUsageTotal = 0.0;
    deque<pair<string, double>> recentUsage;

    chrono::system_clock::time_point t = chrono::system_clock::now();
    t += chrono::seconds(60);

//    SimpleApp  app;

    crow::App<crow::CORSHandler> app;
    auto& cors = app.get_middleware<crow::CORSHandler>();
    cors
            .global()
            .headers("X-Custom-Header", "Upgrade-Insecure-Requests")
                .methods("POST"_method, "GET"_method)
            .prefix("/cors")
                .origin("http://76.23.247.67")
            .prefix("/nocors")
                .ignore();


    CROW_ROUTE(app, "/api/v1/elena/random-node")
            (
                    []() {
                        json::wvalue w = jsonFromNode(main_hiking_graph->getRandomNode());
                        return response(200, w);
                    }
            );

    CROW_ROUTE(app, "/api/v1/elena/two-near-nodes").methods("POST"_method)
            (
                    [](const request &req) {

                        if (!requestRateLimiter->acquireLock()){
                            return response(503);
                        }

                        auto t0 = chrono::high_resolution_clock::now();
                        if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
                            requestRateLimiter->releaseLock();
                            return response(429);
                        }

                        auto r = json::load(req.body);
                        if (!r) {
                            requestRateLimiter->releaseLock();
                            return response(400);
                        }
                        Graph *graph;
                        int edgeBased, maxDistance;
                        try {
                            graph = graphTypeToGraph(r["graph"].s());
                            edgeBased = r["edgeBased"].i();
                            if (edgeBased < 0 || edgeBased > 2){
                                requestRateLimiter->releaseLock();
                                return response(400);
                            }
                            maxDistance = r["maxDistance"].i(); // in cm
                        } catch (const std::exception &e) {
                            requestRateLimiter->releaseLock();
                            return response(400);
                        }
                        pair<Node*, Node*> p = graph->getTwoNearNodes(maxDistance);
                        json::wvalue w;
                        w["lon1"] = p.first->lon;
                        w["lat1"] = p.first->lat;
                        w["lon2"] = p.second->lon;
                        w["lat2"] = p.second->lat;

                        auto t1 = chrono::high_resolution_clock::now();


                        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                        cout << "two near nodes " << t << endl;

                        requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
                        requestRateLimiter->releaseLock();
                        requestRateLimiter->releaseLock();
                        return response(200, w);
                    }
            );

    CROW_ROUTE(app, "/api/v1/elena/shortest-path").methods("POST"_method)
            (
                    [](const request &req) {

                        if (!requestRateLimiter->acquireLock()){
                            return response(503);
                        }

                        auto t0 = chrono::high_resolution_clock::now();
                        if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
                            requestRateLimiter->releaseLock();
                            return response(429);
                        }


                        auto r = json::load(req.body);
                        if (!r) {
                            requestRateLimiter->releaseLock();
                            return response(400);
                        }

                        Graph *graph;
                        double lon1, lat1, lon2, lat2;
                        int edgeBased;
                        try {
                            graph = graphTypeToGraph(r["graph"].s());
                            lon1 = r["lon1"].d();
                            lat1 = r["lat1"].d();
                            lon2 = r["lon2"].d();
                            lat2 = r["lat2"].d();
                            edgeBased = r["edgeBased"].i();
                            if (edgeBased < 0 || edgeBased > 2){
                                requestRateLimiter->releaseLock();
                                return response(400);
                            }
                        } catch (const std::exception &e) {
                            requestRateLimiter->releaseLock();
                            return response(400);
                        }

                        Node *start = graph->nearestNode({lon1, lat1});
                        Node *end = graph->nearestNode({lon2, lat2});



                        if (start == end){
                            json::wvalue w;
                            w["result"] = "_";
                            requestRateLimiter->releaseLock();
                            return response(200, w);
                        }
                        Path *path;
                        if (edgeBased == 1) {
                            path = edgeBasedDijkstraAlgorithm(
                                    graph, start, end, Node::Edge::getLength, nullptr);
                        } else {
                            path = dijkstraAlgorithm(
                                    graph, start, end, Node::Edge::getLength, nullptr);
                        }
                        json::wvalue w = pathToGeoJson(path, graph->edgeGeometries);
                        auto t1 = chrono::high_resolution_clock::now();


                        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                        cout << "shortest path " << t << endl;

                        requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
                        requestRateLimiter->releaseLock();
                        return response(200, w);
                    }
            );



    CROW_ROUTE(app, "/api/v1/elena/elena-minimize").methods("POST"_method)
        (
            [](const request &req) {

                if (!requestRateLimiter->acquireLock()){
                    return response(503);
                }

                auto t0 = chrono::high_resolution_clock::now();
                if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
                    requestRateLimiter->releaseLock();
                    return response(429);
                }


                auto r = json::load(req.body);
                if (!r) {
                    requestRateLimiter->releaseLock();
                    return response(400);
                }
                Graph *graph;
                double lon1, lat1, lon2, lat2, maxLengthRatio;
                int maxLength, edgeBased;
                try {
                    graph = graphTypeToGraph(r["graph"].s());
                    lon1 = r["lon1"].d();
                    lat1 = r["lat1"].d();
                    lon2 = r["lon2"].d();
                    lat2 = r["lat2"].d();
                    maxLengthRatio = r["maxLengthRatio"].d();
                    if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
                        requestRateLimiter->releaseLock();
                        return response(400, "Too big value for \"maxLengthRatio\".");
                    } else if (0 < maxLengthRatio && maxLengthRatio < 1){
                        requestRateLimiter->releaseLock();
                        return response(400, "Too small value for \"maxLengthRatio\".");
                    }
                    double _maxLength = r["maxLength"].d();
                    maxLength = _maxLength < 0 ? INT_MAX : (int) (_maxLength * 1000);
                    edgeBased = r["edgeBased"].i();
                    if (edgeBased < 0 || edgeBased > 2){
                        requestRateLimiter->releaseLock();
                        return response(400);
                    }
                } catch (const std::exception &e) {
                    requestRateLimiter->releaseLock();
                    return response(400);
                }

                Node *start = graph->nearestNode({lon1, lat1});
                Node *end = graph->nearestNode({lon2, lat2});
                if (start == end){
                    json::wvalue w;
                    w["result"] = "_";
                    requestRateLimiter->releaseLock();
                    return response(200, w);
                }
                Path *path;


                if (edgeBased == 1) {
                    path = elenaPathFindMinUsingEdgeBasedDijkstra(graph, start, end, maxLengthRatio, maxLength);
                } else {
                    path = elenaPathFindMinUsingDijkstra(graph, start, end, maxLengthRatio, maxLength);
                }
                if (path->empty()){
                    json::wvalue w;
                    w["result"] = "x";
                    requestRateLimiter->releaseLock();
                    return response(200, w);
                }
                json::wvalue w = pathToGeoJson(path, graph->edgeGeometries);
                auto t1 = chrono::high_resolution_clock::now();


                auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                cout << "elena-minimize path " << t << " seconds" << endl;

                requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
                requestRateLimiter->releaseLock();
                return response(200, w);
            }
        );


    CROW_ROUTE(app, "/api/v1/elena/elena-maximize").methods("POST"_method)(
                    [](const request &req) {

                        if (!requestRateLimiter->acquireLock()){
                            return response(503);
                        }

                        auto t0 = chrono::high_resolution_clock::now();
                        if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
                            requestRateLimiter->releaseLock();
                            return response(429);
                        }


                        auto r = json::load(req.body);
                        if (!r) {
                            requestRateLimiter->releaseLock();
                            return response(400);
                        }
                        Graph *graph;
                        double lon1, lat1, lon2, lat2, maxLengthRatio;
                        long numProduce, numMaxSelect, numEpoch;
                        int maxLength, edgeBased;
                        DuplicateEdge duplicateEdge;
                        int seconds;
                        try {
                            graph = graphTypeToGraph(r["graph"].s());
                            lon1 = r["lon1"].d();
                            lat1 = r["lat1"].d();
                            lon2 = r["lon2"].d();
                            lat2 = r["lat2"].d();
                            maxLengthRatio = r["maxLengthRatio"].d();
                            if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
                                requestRateLimiter->releaseLock();
                                return response(400, "Too big value for \"maxLengthRatio\".");
                            } else if (0 < maxLengthRatio && maxLengthRatio < 1){
                                requestRateLimiter->releaseLock();
                                return response(400, "Too small value for \"maxLengthRatio\".");
                            }
                            numProduce = r["numProduce"].i();
                            numMaxSelect = r["numMaxSelect"].i();
                            numEpoch = r["numEpoch"].i();
                            int duplicateEdgeInt = r["duplicateEdges"].i();
                            if (duplicateEdgeInt < 0 || duplicateEdgeInt >= static_cast<int>(duplicateEdgeEnumCount)){
                                requestRateLimiter->releaseLock();
                                return response(400);
                            }
                            duplicateEdge = static_cast<DuplicateEdge>(duplicateEdgeInt);
                            double _maxLength = r["maxLength"].d();
                            maxLength = _maxLength < 0 ? INT_MAX : (int) (_maxLength * 1000);
                            edgeBased = r["edgeBased"].i();
                            if (edgeBased < 0 || edgeBased > 2){
                                requestRateLimiter->releaseLock();
                                return response(400);
                            }
                            seconds = r["seconds"].i();
                        } catch (const std::exception &e) {
                            requestRateLimiter->releaseLock();
                            return response(400);
                        }

                        Node *start = graph->nearestNode({lon1, lat1});
                        Node *end = graph->nearestNode({lon2, lat2});
                        if (start == end){
                            json::wvalue w;
                            w["result"] = "_";
                            requestRateLimiter->releaseLock();
                            return response(200, w);
                        }

                        Path* path;
                        if (edgeBased == 1){
                            path = elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
                                    graph, start, end, maxLengthRatio, numProduce, numMaxSelect,
                                    numEpoch, duplicateEdge, maxLength, seconds * 900
                            )->toPath();
                        } else {
                            path = elenaPathSearchMaxUsingGeneticAlgorithm(
                                    graph, start, end, maxLengthRatio, numProduce, numMaxSelect,
                                    numEpoch, duplicateEdge, maxLength, seconds * 900
                            )->toPath();
                        }

                        if (path->empty()){
                            cout << "path not found within maxLength " << maxLength << endl;
                            json::wvalue w;
                            w["result"] = "x";
                            requestRateLimiter->releaseLock();
                            return response(200, w);
                        }
                        unordered_set<Node::Edge*> uset(path->edges.begin(), path->edges.end());
                        unordered_set<long> uset2;
                        for (Node::Edge* edge: path->edges){
                            uset2.insert(edge->id);
                        }
                        cout << duplicateEdge << " " << path->edges.size() << " " << uset2.size() << " " << uset.size() << endl;

                        json::wvalue w = pathToGeoJson(path, graph->edgeGeometries);
                        auto t1 = chrono::high_resolution_clock::now();


                        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
                        cout << "elena-maximize path " << t << " seconds" << endl;

                        requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
                        requestRateLimiter->releaseLock();
                        return response(200, w);
                    }
            );


    app.port(port).multithreaded().run();
}