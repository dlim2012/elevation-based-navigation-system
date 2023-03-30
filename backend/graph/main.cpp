#include <iostream>
#include <vector>


#include "include/graph.h"

using namespace std;

extern "C"
{
    __declspec(dllexport) void dummy(char* s, int n){
        string str = s;
        cout << str << " " << n << endl;
    }

    __declspec(dllexport) Graph* createGraph(char* _version){

        string version = _version;
        cout << version << endl;
        return new Graph(version, true);
    }

    __declspec(dllexport) size_t shortestPath(
            Graph* graph,
            long start_id,
            long end_id,
            long* array_out,
            double* coords,
            double* length
    ){
        cout << start_id << " " << end_id << endl;
        cout << graph->nodes.size() << endl;

        pair<double, vector<Node*>> res = graph->shortestPath(start_id, end_id);
        vector<Node*> path = res.second;
        for (int i=0; i<path.size(); i++){
            Node* node = path[i];
            array_out[i] = node->id;
            coords[2*i] = node->x;
            coords[2*i+1] = node->y;
        }
        length[0] = res.first;
        return path.size();
    }
    __declspec(dllexport) long getRandomNodeId(Graph* graph, int seed){
        srand(seed);
        auto it = graph->nodes.cbegin();
        int random = rand() % graph->nodes.size();
        std::advance(it, random);
        return it->first;
    }

    __declspec(dllexport) long getNearestNode(Graph* graph, double* coords){
        pair<double, double> target = {coords[0], coords[1]};
        Node* res = graph->nearestNode(target);
        coords[0] = res->x;
        coords[1] = res->y;
        return res->id;
    }

    __declspec(dllexport) size_t findRandomPath(
            Graph* graph,
            long start_id,
            long end_id,
            long* array_out,
            double* coords,
            double* length,
            int seed
    ){

        Node* start(graph->nodes[start_id]), *end(graph->nodes[end_id]);
        pair<double, vector<Node*>> res = graph->findRandomPath(start, end, *length, seed);
        vector<Node*> path = res.second;
        for (int i=0; i<path.size(); i++){
            Node* node = path[i];
            array_out[i] = node->id;
            coords[2*i] = node->x;
            coords[2*i+1] = node->y;
        }
        length[0] = res.first;
        return path.size();
    }
};

int main() {
    cout << "main" << endl;
    auto* g = new Graph("maps/test", false);
//    auto* g = new Graph("test.maxGroup", true);

    for (auto it: g->nodes){
        cout << it.first << endl;
    }

    srand(time(0));

    long start_id(1364765716), end_id(3350088311);
    start_id = getRandomNodeId(g, rand());
    end_id = getRandomNodeId(g, rand());
    pair<double, vector<Node*>> res = g->shortestPath(start_id, end_id);
//    for (Node* node: res.second){
//        cout << node->id << " " ;
//    }
//    cout << endl;
//
    cout << "(shortest path) length: " << res.first << endl;


    Node* start(g->nodes[start_id]), *end(g->nodes[end_id]);
    res = g->findRandomPath(start, end, 99999999999.0, 0);
//    for (Node* node: res.second){
//        cout << node->id << " " ;
//    }
    cout << endl;

    cout << "(random path) length: " << res.first << endl;

    return 0;
}

//https://stephenscotttucker.medium.com/interfacing-python-with-c-using-ctypes-classes-and-arrays-42534d562ce7
