#include <string>
#include <iostream>
#include <vector>
#include <chrono>
#include "include/pathfinding.h"

extern "C"
{

    __declspec(dllexport) Graph* createGraph(char* _map_name, int* size){

        string map_name = _map_name;
        cout << map_name << endl;
        Graph* graph = new Graph(map_name, false);

        unordered_set<long> visited;
        for (pair<long, Node*> pair1: graph->nodes){
            Node* node = pair1.second;
            for (pair<long, Node::Edge*> pair2: node->edges){
                Node::Edge* edge = pair2.second;
                if (visited.find(edge->id) == visited.end()){
                    visited.insert(edge->id);
                }
            }
        }
        *size = visited.size();
        return graph;
    }

    __declspec(dllexport) int getAllEdges(Graph* graph, double* edges){
        unordered_set<long> visited;

        for (pair<long, Node*> pair1: graph->nodes){
            Node* node = pair1.second;
            for (pair<long, Node::Edge*> pair2: node->edges){
                Node::Edge* edge = pair2.second;
                if (visited.find(edge->id) == visited.end()){
                    edges[visited.size()*4] = edge->u->lon;
                    edges[visited.size()*4+1] = edge->u->lat;
                    edges[visited.size()*4+2] = edge->v->lon;
                    edges[visited.size()*4+3] = edge->v->lat;
                    visited.insert(edge->id);
                }
            }
        }

        return visited.size();
    }

    __declspec(dllexport) long getRandomNodeId(Graph* graph, int seed){
        srand(seed);
        return graph->getRandomNodeId();
    }

    __declspec(dllexport) long getNearestNode(Graph* graph, double* coords){
        pair<double, double> v = {coords[0], coords[1]};
        Node* res = graph->nearestNode(v);
        coords[0] = res->lon;
        coords[1] = res->lat;
        return res->id;
    }

    __declspec(dllexport) size_t aStarAlgorithm(
            Graph* graph,
            long start_id,
            long end_id,
            long* array_out,
            double* coords,
            double* length,
            double* elevation
            ){
        auto t0 = chrono::high_resolution_clock::now();

        Node *start(graph->nodes[start_id]), *end(graph->nodes[end_id]);
        Path* res = aStarAlgorithm(graph, start, end, Node::Edge::getLength);
        vector<Path::PathEdge*> pathEdges = res->pathEdges;
        array_out[0] = start_id;
        coords[0] = start->lon;
        coords[1] = start->lat;
        for (int i=0; i<pathEdges.size(); i++){
            Node* node = pathEdges[i]->edge->v;
            array_out[i+1] = node->id;
            coords[2*i+2] = node->lon;
            coords[2*i+3] = node->lat;
        }
        length[0] = res->length;
        elevation[0] = res->elevation;

        auto t1 = chrono::high_resolution_clock::now();
        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "A-star algorithm (time: " << t << ")" << endl;
        return pathEdges.size();
    }

    __declspec(dllexport) size_t dijkstraAlgorithm(
            Graph* graph,
            long start_id,
            long end_id,
            long* array_out,
            double* coords,
            double* length,
            double* elevation
    ){
        auto t0 = chrono::high_resolution_clock::now();

        Node *start(graph->nodes[start_id]), *end(graph->nodes[end_id]);
        Path* res = dijkstraAlgorithm(graph, start, end, Node::Edge::getLength);
        vector<Path::PathEdge*> pathEdges = res->pathEdges;
        array_out[0] = start_id;
        coords[0] = start->lon;
        coords[1] = start->lat;
        for (int i=0; i<pathEdges.size(); i++){
            Node* node = pathEdges[i]->edge->v;
            array_out[i+1] = node->id;
            coords[2*i+2] = node->lon;
            coords[2*i+3] = node->lat;
        }
        length[0] = res->length;
        elevation[0] = res->elevation;


        auto t1 = chrono::high_resolution_clock::now();
        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "Dijkstra algorithm (time: " << t << ")" << endl;

        return pathEdges.size();
    }

    __declspec(dllexport) size_t findElevationBasedPath(
            Graph* graph,
            long start_id,
            long end_id,
            double maxWeightRatio,
            bool minimize,
            long* array_out,
            double* coords,
            double* length,
            double* elevation
    ){
        auto t0 = chrono::high_resolution_clock::now();

        Node *start(graph->nodes[start_id]), *end(graph->nodes[end_id]);
        Path* res = findElevationBasedPath(graph, start, end, maxWeightRatio, minimize);
        vector<Path::PathEdge*> pathEdges = res->pathEdges;
        array_out[0] = start_id;
        coords[0] = start->lon;
        coords[1] = start->lat;
        for (int i=0; i<pathEdges.size(); i++){
            Node* node = pathEdges[i]->edge->v;
            array_out[i+1] = node->id;
            coords[2*i+2] = node->lon;
            coords[2*i+3] = node->lat;
        }
        length[0] = res->length;
        elevation[0] = res->elevation;

        auto t1 = chrono::high_resolution_clock::now();
        auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;
        cout << "Elena algorithm (time: " << t << ")" << endl;

        return pathEdges.size();
    }
};

