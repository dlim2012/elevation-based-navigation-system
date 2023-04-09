
#include <iostream>

#include "crow.h"
#include "restApi.h"

using namespace crow;
using namespace std;

Graph *graph;
unordered_map<long, vector<pair<double, double>>> umap;

const int MAXIMUM_MAX_LENGTH_RATIO = 10;

crow::json::wvalue jsonFromPath(
        Path *path, unordered_map<long, vector<pair<double, double>>> &umap) {
    vector<Node::Edge *> &edges = path->edges;
    crow::json::wvalue w;
    if (edges.empty()) {
        return w;
    }

    int index = 0;
    Node* start = path->getStart();
    w["path"][index]["lon"] = start->lon;
    w["path"][index]["lat"] = start->lat;
    index++;
    int count = 0;
    for (Node::Edge *edge: edges) {
        vector<pair<double, double>> &v = umap[edge->id];
        for (int i = 1; i < v.size(); i++) {
            w["path"][index]["lon"] = v[i].first;
            w["path"][index]["lat"] = v[i].second;
            index++;
        }
    }
    w["elevation"] = path->elevation;
    w["length"] = path->length;
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

void runServer(string connectionString, int port, MapConfig mapConfig) {

    graph = new Graph(connectionString, mapConfig);
    umap = allEdges(connectionString);


    SimpleApp app;


    CROW_ROUTE(app, "/api/v1/elena/random-node")
            (
                    []() {
                        json::wvalue w = jsonFromNode(graph->getRandomNode());
                        return response(200, w);
                    }
            );

    CROW_ROUTE(app, "/api/v1/elena/nearest-node").methods("POST"_method)
        (
            [](const request &req) {
                auto r = json::load(req.body);
                if (!r)
                    return response(400);
                double lon, lat;
                try {
                    lon = r["lon"].d();
                    lat = r["lat"].d();
                } catch(const std::exception& e){
                    return response(400);
                }

                Node *node = graph->nearestNode({lon, lat});

                json::wvalue w = jsonFromNode(node);
                return response(200, w);
            }
    );

    CROW_ROUTE(app, "/api/v1/elena/shortest-path").methods("POST"_method)
            (
                    [](const request &req) {
                        auto r = json::load(req.body);
                        if (!r)
                            return response(400);
                        double lon1, lat1, lon2, lat2;
                        try {
                            lon1 = r["lon1"].d();
                            lat1 = r["lat1"].d();
                            lon2 = r["lon2"].d();
                            lat2 = r["lat2"].d();
                        } catch (const std::exception &e) {
                            return response(400);
                        }

                        Node *start = graph->nearestNode({lon1, lat1});
                        Node *end = graph->nearestNode({lon2, lat2});
                        Path *path = edgeBasedDijkstraAlgorithm(
                                graph, start, end, Node::Edge::getLength, false, nullptr);
                        json::wvalue w = jsonFromPath(path, umap);
                        return response(200, w);
                    }
            );

    CROW_ROUTE(app, "/api/v1/elena/shortest-path/random").methods("GET"_method)
                    (
                            [](){
                                Node *start = graph->getRandomNode();
                                Node *end = graph->getRandomNode();
                                Path *path = edgeBasedDijkstraAlgorithm(
                                        graph, start, end, Node::Edge::getLength, false, nullptr);
                                cout << path->edges.size() << endl;
                                json::wvalue w = jsonFromPath(path, umap);
                                return response(200, w);
                            }
                    );

    CROW_ROUTE(app, "/api/v1/elena/elena-minimize/").methods("POST"_method)
        (
            [](const request &req) {
                auto r = json::load(req.body);
                if (!r)
                    return response(400);
                double lon1, lat1, lon2, lat2, maxLengthRatio;
                try {
                    lon1 = r["lon1"].d();
                    lat1 = r["lat1"].d();
                    lon2 = r["lon2"].d();
                    lat2 = r["lat2"].d();
                    maxLengthRatio = r["maxLengthRatio"].d();
                } catch (const std::exception &e) {
                    return response(400);
                }

                if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
                    return response(400, "Too big value for \"x\".");
                }
                if (maxLengthRatio < 1){
                    return response(400, "Too small value for \"x\".");
                }

                Node *start = graph->nearestNode({lon1, lat1});
                Node *end = graph->nearestNode({lon2, lat2});
                Path *path = elenaPathFindMinUsingDijkstra(graph, start, end, maxLengthRatio);
                json::wvalue w = jsonFromPath(path, umap);
                return response(200, w);
            }
        );

    CROW_ROUTE(app, "/api/v1/elena/elena-minimize/random/max-length-ratio=<double>").methods("GET"_method)
            (
                    [](double maxLengthRatio) {
//                        auto r = json::load(req.body);
//                        if (!r)
//                            return response(400);
//                        double maxLengthRatio;

//                        try {
//                            maxLengthRatio = r["maxLengthRatio"].d();
//                        } catch (const std::exception &e) {
//                            return response(400);
//                        }

                        if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
                            return response(400, "Too big value for max length ratio.");
                        }
                        if (maxLengthRatio < 1){
                            return response(400, "Too small value for max length ratio.");
                        }

                        Node *start = graph->getRandomNode();
                        Node *end = graph->getRandomNode();

                        Path *path = elenaPathFindMinUsingDijkstra(graph, start, end, maxLengthRatio);

                        json::wvalue w = jsonFromPath(path, umap);
                        return response(200, w);
                    }
            );

    CROW_ROUTE(app, "/api/v1/elena/elena-maximize").methods("POST"_method)(
                    [](const request &req) {
                        auto r = json::load(req.body);
                        if (!r)
                            return response(400);
                        double lon1, lat1, lon2, lat2, maxLengthRatio;
                        long numProduce, numMaxSelect, numEpoch;
                        try {
                            lon1 = r["lon1"].d();
                            lat1 = r["lat1"].d();
                            lon2 = r["lon2"].d();
                            lat2 = r["lat2"].d();
                            maxLengthRatio = r["maxLengthRatio"].d();
                            numProduce = r["numProduce"].i();
                            numMaxSelect = r["numMaxSelect"].i();
                            numEpoch = r["numEpoch"].i();
                        } catch (const std::exception &e) {
                            return response(400);
                        }

                        if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
                            return response(400, "Too big value for \"maxLengthRatio\".");
                        }
                        if (maxLengthRatio < 1){
                            return response(400, "Too small value for \"maxLengthRatio\".");
                        }

                        Node *start = graph->nearestNode({lon1, lat1});
                        Node *end = graph->nearestNode({lon2, lat2});
                        Path* path = elenaPathSearchMaxUsingGeneticAlgorithm(
                                graph, start, end, maxLengthRatio, numProduce, numMaxSelect, numEpoch)->toPath();
                        json::wvalue w = jsonFromPath(path, umap);
                        return response(200, w);
                    }
            );


    CROW_ROUTE(app, "/api/v1/elena/elena-maximize/random/max-length-ratio=<double>").methods("GET"_method)(
                    [](double maxLengthRatio) {

                        if (maxLengthRatio > MAXIMUM_MAX_LENGTH_RATIO){
                            return response(400, "Too big value for \"maxLengthRatio\".");
                        }
                        if (maxLengthRatio < 1){
                            return response(400, "Too small value for \"maxLengthRatio\".");
                        }

                        Node *start = graph->getRandomNode();
                        Node *end = graph->getRandomNode();
                        Path* path = elenaPathSearchMaxUsingGeneticAlgorithm(graph, start, end, maxLengthRatio, 30, 10, 10)->toPath();
                        json::wvalue w = jsonFromPath(path, umap);
                        return response(200, w);
                    }
            );

    app.port(port).multithreaded().run();
}