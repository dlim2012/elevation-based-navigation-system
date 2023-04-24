
#include <iostream>
#include <chrono>
#include <deque>
#include <thread>
#include <future>

#include "crow.h"
#include "restApi.h"
#include "crow/middlewares/cors.h"
#include "RequestRateLimiter.h"


using namespace crow;
using namespace std;

typedef std::chrono::time_point<std::chrono::system_clock> tp;

const int LOCK_COUNT = 16;
const int MAX_REQUESTS = 4;

const int DEFAULT_NUM_PRODUCE = 100;
const int DEFAULT_NUM_MAX_SELECT = 30;
const int DEFAULT_NUM_EPOCH = 1000;

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

crow::json::wvalue emptyGeoJson(Node* start){
    crow::json::wvalue w;
    w["GeoJSON"]["type"] = "Node";
    w["GeoJSON"]["coordinates"][0][0] = start->lat;
    w["GeoJSON"]["coordinates"][0][1] = start->lon;
    w["elevation"] = 0;
    w["length"] = 0;  // cm to km
    return w;
}

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

vector<int> readVectorInt(const json::rvalue& r, string key, int minVal, int maxVal){
    auto _res = r[key];
    vector<json::rvalue> rvalues = _res.lo();
    vector<int> res;
    for (json::rvalue rv: rvalues){
        int val = (int) rv.i();
        if (minVal > val || maxVal < val){
            res.clear();
            return res;
        }
        res.push_back(val);
    }
    return res;
}

vector<double> readVectorDouble(const json::rvalue& r, string key, double minVal, double maxVal){
    auto _res = r[key];
    vector<json::rvalue> rvalues = _res.lo();
    vector<double> res;
    for (json::rvalue rv: rvalues){
        double val = rv.d();
        if (minVal > val || maxVal < val){
            res.clear();
            return res;
        }
        res.push_back(val);
    }
    return res;
}

void searchPath(promise<json::wvalue>&& p, Graph* graph, int nodeId, double lon1, double lat1, double lon2,
                double lat2, int pathType, int edgeBased, double maxLengthRatio,
                DuplicateEdge duplicateEdge, int seconds){
    json::wvalue w;
    w["nodeId"] = nodeId;

    Node *start = graph->nearestNode({lon1, lat1});
    Node *end = graph->nearestNode({lon2, lat2});

    Path *path;
    if (start != end){
        if (pathType == 0){
            if (edgeBased == 1) {
                path = edgeBasedDijkstraAlgorithm(
                        graph, start, end, Node::Edge::getLength, nullptr);
            } else {
                path = dijkstraAlgorithm(
                        graph, start, end, Node::Edge::getLength, nullptr);
            }
        } else if (pathType == 1){
            if (edgeBased == 1) {
                path = elenaPathFindMinUsingEdgeBasedDijkstra(graph, start, end, maxLengthRatio);
            } else {
                path = elenaPathFindMinUsingDijkstra(graph, start, end, maxLengthRatio);
            }
        } else {
            if (edgeBased == 1){
                path = elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
                        graph, start, end, maxLengthRatio, DEFAULT_NUM_PRODUCE, DEFAULT_NUM_MAX_SELECT,
                        DEFAULT_NUM_EPOCH, duplicateEdge, seconds * 900, 10
                )->toPath();
            } else {
                path = elenaPathSearchMaxUsingGeneticAlgorithm(
                        graph, start, end, maxLengthRatio, DEFAULT_NUM_PRODUCE, DEFAULT_NUM_MAX_SELECT,
                        DEFAULT_NUM_EPOCH, duplicateEdge, seconds * 900, 10
                )->toPath();
            }
        }
        w["path"] = pathToGeoJson(path, graph->edgeGeometries);
    } else {
        w["path"] = emptyGeoJson(start);
    }
    p.set_value(w);
}



void runServer(string& connectionString, int port) {
    requestRateLimiter = new RequestRateLimiter(LOCK_COUNT, MAX_REQUESTS);

    sharedEdgesGeometries = new unordered_map<long, Geometry*>();
//    cout << "Simple driveway" << endl;
//    main_motorway_graph = new Graph(connectionString, motorway);
    cout << "Hiking" << endl;
    main_hiking_graph = new Graph(connectionString, hiking, US_mainland, sharedEdgesGeometries);
//    cout << endl << "Cycling" << endl;
//    main_cycling_graph = new Graph(connectionString, cycling_with_restrictions);
    cout << endl << "Cycling with restrictions" << endl;
    main_cycling_graph = new Graph(connectionString, cycling_with_restrictions, US_mainland, sharedEdgesGeometries);

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
//                .origin("http://192.168.1.20")
//                .origin("http://localhost")
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



                        auto r = json::load(req.body);
                        if (!r) {
                            return response(400);
                        }
                        Graph *graph;
                        int edgeBased, maxDistance;
                        try {
                            graph = graphTypeToGraph(r["graph"].s());
                            edgeBased = r["edgeBased"].i();
                            if (edgeBased < 0 || edgeBased > 2){
                                return response(400);
                            }
                            maxDistance = r["maxDistance"].i(); // in cm
                        } catch (const std::exception &e) {
                            return response(400);
                        }

                        int requestId = requestRateLimiter->acquireLock(req.remote_ip_address, 1);

                        if (requestId == -3){
                            return response(400);
                        } else if (requestId == -2){
                            return response(429);
                        } else if (requestId == -1) {
                            return response(503);
                        }

                        pair<Node*, Node*> p = graph->getTwoNearNodes(maxDistance);
                        json::wvalue w;
                        w["lon1"] = p.first->lon;
                        w["lat1"] = p.first->lat;
                        if (p.second != nullptr){
                            w["lon2"] = p.second->lon;
                            w["lat2"] = p.second->lat;
                        } else {
                            w["lon2"] = p.first->lon;
                            w["lat2"] = p.first->lat;
                        }


                        requestRateLimiter->releaseLock(req.remote_ip_address, requestId);
                        return response(200, w);
                    }
            );

//    CROW_ROUTE(app, "/api/v1/elena/shortest-path").methods("POST"_method)
//            (
//                    [](const request &req) {
//
//                        if (!requestRateLimiter->acquireLock(1)){
//                            return response(503);
//                        }
//
//                        auto t0 = chrono::high_resolution_clock::now();
//                        if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
//                            requestRateLimiter->releaseLock(1);
//                            return response(429);
//                        }
//
//
//                        auto r = json::load(req.body);
//                        if (!r) {
//                            requestRateLimiter->releaseLock(1);
//                            return response(400);
//                        }
//
//                        Graph *graph;
//                        double lon1, lat1, lon2, lat2;
//                        int edgeBased;
//                        try {
//                            graph = graphTypeToGraph(r["graph"].s());
//                            lon1 = r["lon1"].d();
//                            lat1 = r["lat1"].d();
//                            lon2 = r["lon2"].d();
//                            lat2 = r["lat2"].d();
//                            edgeBased = r["edgeBased"].i();
//                            if (edgeBased < 0 || edgeBased > 2){
//                                throw runtime_error("");
//                            }
//                        } catch (const std::exception &e) {
//                            requestRateLimiter->releaseLock(1);
//                            return response(400);
//                        }
//
//                        Node *start = graph->nearestNode({lon1, lat1});
//                        Node *end = graph->nearestNode({lon2, lat2});
//
//
//
//                        if (start == end){
//                            json::wvalue w;
//                            w["result"] = "_";
//                            requestRateLimiter->releaseLock(1);
//                            return response(200, w);
//                        }
//                        Path *path;
//                        if (edgeBased == 1) {
//                            path = edgeBasedDijkstraAlgorithm(
//                                    graph, start, end, Node::Edge::getLength, nullptr);
//                        } else {
//                            path = dijkstraAlgorithm(
//                                    graph, start, end, Node::Edge::getLength, nullptr);
//                        }
//                        json::wvalue w = pathToGeoJson(path, graph->edgeGeometries);
//                        auto t1 = chrono::high_resolution_clock::now();
//
//
//                        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//                        cout << "shortest path " << t << endl;
//
//                        requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
//                        requestRateLimiter->releaseLock(1);
//                        return response(200, w);
//                    }
//            );
//
//
//
//    CROW_ROUTE(app, "/api/v1/elena/elena-minimize").methods("POST"_method)
//        (
//            [](const request &req) {
//
//                if (!requestRateLimiter->acquireLock(1)){
//                    return response(503);
//                }
//
//                auto t0 = chrono::high_resolution_clock::now();
//                if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
//                    requestRateLimiter->releaseLock(1);
//                    return response(429);
//                }
//
//
//                auto r = json::load(req.body);
//                if (!r) {
//                    requestRateLimiter->releaseLock(1);
//                    return response(400);
//                }
//                Graph *graph;
//                double lon1, lat1, lon2, lat2, maxLengthRatio;
//                int maxLength, edgeBased;
//                try {
//                    graph = graphTypeToGraph(r["graph"].s());
//                    lon1 = r["lon1"].d();
//                    lat1 = r["lat1"].d();
//                    lon2 = r["lon2"].d();
//                    lat2 = r["lat2"].d();
//                    maxLengthRatio = r["maxLengthRatio"].d();
//                    if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
//                        throw runtime_error("");
//                    } else if (0 < maxLengthRatio && maxLengthRatio < 1){
//                        throw runtime_error("");
//                    }
//                    double _maxLength = r["maxLength"].d();
//                    maxLength = _maxLength < 0 ? INT_MAX : (int) (_maxLength * 1000);
//                    edgeBased = r["edgeBased"].i();
//                    if (edgeBased < 0 || edgeBased > 2){
//                        throw runtime_error("");
//                    }
//                } catch (const std::exception &e) {
//                    requestRateLimiter->releaseLock(1);
//                    return response(400);
//                }
//
//                Node *start = graph->nearestNode({lon1, lat1});
//                Node *end = graph->nearestNode({lon2, lat2});
//                if (start == end){
//                    json::wvalue w;
//                    w["result"] = "_";
//                    requestRateLimiter->releaseLock(1);
//                    return response(200, w);
//                }
//                Path *path;
//
//
//                if (edgeBased == 1) {
//                    path = elenaPathFindMinUsingEdgeBasedDijkstra(graph, start, end, maxLengthRatio);
//                } else {
//                    path = elenaPathFindMinUsingDijkstra(graph, start, end, maxLengthRatio);
//                }
//                if (path->empty()){
//                    json::wvalue w;
//                    w["result"] = "x";
//                    requestRateLimiter->releaseLock(1);
//                    return response(200, w);
//                }
//                json::wvalue w = pathToGeoJson(path, graph->edgeGeometries);
//                auto t1 = chrono::high_resolution_clock::now();
//
//
//                auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//                cout << "elena-minimize path " << t << " seconds" << endl;
//
//                requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
//                requestRateLimiter->releaseLock(1);
//                return response(200, w);
//            }
//        );
//
//
//    CROW_ROUTE(app, "/api/v1/elena/elena-maximize").methods("POST"_method)(
//                    [](const request &req) {
//
//                        if (!requestRateLimiter->acquireLock(1)){
//                            return response(503);
//                        }
//
//                        auto t0 = chrono::high_resolution_clock::now();
//                        if (!requestRateLimiter->isAllowed(req.remote_ip_address, t0)){
//                            requestRateLimiter->releaseLock(1);
//                            return response(429);
//                        }
//
//
//                        auto r = json::load(req.body);
//                        if (!r) {
//                            requestRateLimiter->releaseLock(1);
//                            return response(400);
//                        }
//                        Graph *graph;
//                        double lon1, lat1, lon2, lat2, maxLengthRatio;
//                        long numProduce, numMaxSelect, numEpoch;
//                        int edgeBased;
//                        DuplicateEdge duplicateEdge;
//                        int seconds;
//                        try {
//                            graph = graphTypeToGraph(r["graph"].s());
//                            lon1 = r["lon1"].d();
//                            lat1 = r["lat1"].d();
//                            lon2 = r["lon2"].d();
//                            lat2 = r["lat2"].d();
//                            maxLengthRatio = r["maxLengthRatio"].d();
//                            if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
//                                throw runtime_error("");
//                            } else if (0 < maxLengthRatio && maxLengthRatio < 1){
//                                throw runtime_error("");
//                            }
//                            numProduce = r["numProduce"].i();
//                            numMaxSelect = r["numMaxSelect"].i();
//                            numEpoch = r["numEpoch"].i();
//                            int duplicateEdgeInt = r["duplicateEdges"].i();
//                            if (duplicateEdgeInt < 0 || duplicateEdgeInt >= static_cast<int>(duplicateEdgeEnumCount)){
//                                throw runtime_error("");
//                            }
//                            duplicateEdge = static_cast<DuplicateEdge>(duplicateEdgeInt);
//                            edgeBased = r["edgeBased"].i();
//                            if (edgeBased < 0 || edgeBased > 2){
//                                throw runtime_error("");
//                            }
//                            seconds = r["seconds"].i();
//                        } catch (const std::exception &e) {
//                            requestRateLimiter->releaseLock(1);
//                            return response(400);
//                        }
//
//                        Node *start = graph->nearestNode({lon1, lat1});
//                        Node *end = graph->nearestNode({lon2, lat2});
//                        if (start == end){
//                            json::wvalue w;
//                            w["result"] = "_";
//                            requestRateLimiter->releaseLock(1);
//                            return response(200, w);
//                        }
//
//                        Path* path;
//                        if (edgeBased == 1){
//                            path = elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
//                                    graph, start, end, maxLengthRatio, numProduce, numMaxSelect,
//                                    numEpoch, duplicateEdge,  seconds * 900, 10
//                            )->toPath();
//                        } else {
//                            path = elenaPathSearchMaxUsingGeneticAlgorithm(
//                                    graph, start, end, maxLengthRatio, numProduce, numMaxSelect,
//                                    numEpoch, duplicateEdge, seconds * 900, 10
//                            )->toPath();
//                        }
//
//                        if (path->empty()){
//                            json::wvalue w;
//                            w["result"] = "x";
//                            requestRateLimiter->releaseLock(1);
//                            return response(200, w);
//                        }
//                        unordered_set<Node::Edge*> uset(path->edges.begin(), path->edges.end());
//                        unordered_set<long> uset2;
//                        for (Node::Edge* edge: path->edges){
//                            uset2.insert(edge->id);
//                        }
//
//                        json::wvalue w = pathToGeoJson(path, graph->edgeGeometries);
//                        auto t1 = chrono::high_resolution_clock::now();
//
//
//                        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
//                        cout << "elena-maximize path " << t << " seconds" << endl;
//
//                        requestRateLimiter->addRecord(req.remote_ip_address, t0, t1);
//                        requestRateLimiter->releaseLock(1);
//                        return response(200, w);
//                    }
//            );





    // Handle multiple requests at once

    CROW_ROUTE(app, "/api/v1/elena/paths").methods("POST"_method)
            (
                    [](const request &req) {

                        auto r = json::load(req.body);
                        if (!r) {
                            return response(400);
                        }
                        Graph *graph;
                        int edgeBased, seconds;
                        DuplicateEdge duplicateEdge;
                        vector<int> nodeId, pathTypes;
                        vector<double> lon1, lat1, lon2, lat2, maxLengthRatio;
                        try {
                            graph = graphTypeToGraph(r["graph"].s());
                            edgeBased = r["edgeBased"].i();
                            seconds = r["seconds"].i();
                            int duplicateEdgeInt = r["duplicateEdges"].i();
                            if (duplicateEdgeInt < 0 || duplicateEdgeInt >= static_cast<int>(duplicateEdgeEnumCount)){
                                throw runtime_error("");
                            }
                            duplicateEdge = static_cast<DuplicateEdge>(duplicateEdgeInt);

                            nodeId = readVectorInt(r, "nodeId", 0, INT_MAX);
                            pathTypes = readVectorInt(r, "pathType", 0, 2);
                            lon1 = readVectorDouble(r, "lon1", -180.0, 180.0);
                            lat1 = readVectorDouble(r, "lat1", -90.0, 90.0);
                            lon2 = readVectorDouble(r, "lon2", -180.0, 180.0);
                            lat2 = readVectorDouble(r, "lat2", -90.0, 90.0);
                            maxLengthRatio = readVectorDouble(r, "maxLengthRatio", 1, INT_MAX);

                            if (nodeId.empty()
                                || nodeId.size() != pathTypes.size()
                                || lon1.size() != nodeId.size()
                                || lat1.size() != nodeId.size()
                                || lon2.size() != nodeId.size()
                                || lat2.size() != nodeId.size()
                                || maxLengthRatio.size() != nodeId.size()
                                || edgeBased < 0 || edgeBased > 2
                                || nodeId.size() > MAX_REQUESTS
                                ){
                                throw runtime_error("");
                            }
                        } catch (const std::exception &e) {
                            return response(400);
                        }

                        int numPaths = nodeId.size();


                        int requestId = requestRateLimiter->acquireLock(req.remote_ip_address, numPaths);

                        if (requestId == -3){
                            return response(400);
                        } else if (requestId == -2){
                            return response(429);
                        } else if (requestId == -1) {
                            return response(503);
                        }

                        json::wvalue w;

                        vector<future<json::wvalue>> futures;
                        vector<thread> threads;


                        for (int i=0; i<numPaths; i++){
                            promise<json::wvalue> _p;
                            future<json::wvalue> _f = _p.get_future();
                            thread _t(&searchPath, std::move(_p), graph, nodeId[i], lon1[i], lat1[i], lon2[i], lat2[i],
                                      pathTypes[i], edgeBased, maxLengthRatio[i], duplicateEdge, seconds);
                            futures.push_back(std::move(_f));
                            threads.push_back(std::move(_t));
                        }

                        for (int i=0; i<numPaths; i++){
                            threads[i].join();
                            w[i] = futures[i].get();
                        }

                        requestRateLimiter->releaseLock(req.remote_ip_address, requestId);
                        return response(200, w);
                    }
            );


    thread requestCleanUpThread(&RequestRateLimiter::repeatCleanRequest, requestRateLimiter, 3600);

    auto _a = app.port(port).multithreaded().run_async();

}
