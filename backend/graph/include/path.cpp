#include <iostream>
#include <vector>
#include <algorithm>

#include "path.h"

using namespace std;


PathEdges::PathEdge::PathEdge(Node::Edge *edge, double length, int elevation) {
    this->edge = edge;
    this->length = length;
    this->elevation = elevation;
}

PathEdges::PathEdges(Node* start) {
    this->start = start;
    this->length = 0.0;
    this->elevation = 0.0;
}

PathEdges::PathEdges(Node* start, vector<Node::Edge *>& edges) {
    this->start = start;
    this->length = 0.0;
    this->elevation = 0.0;
    for (Node::Edge *edge: edges) {
        this->addEdge(edge);
    }
}

PathEdges::PathEdges(Node* start, vector<PathEdge *> &pathEdges, size_t index) {
    this->start = start;
    this->pathEdges = vector<PathEdge*>(pathEdges.begin(), pathEdges.begin() + min(pathEdges.size(), index));

    if (this->pathEdges.empty()) {
        this->length = 0.0;
        this->elevation = 0.0;
    } else {
        this->length = this->pathEdges.back()->length;
        this->elevation = this->pathEdges.back()->elevation;
    }
}

void PathEdges::addEdge(Node::Edge *edge) {
    this->length += edge->length;
    this->elevation += edge->elevation;
    this->pathEdges.push_back(new PathEdge(edge, this->length, this->elevation));
}

int PathEdges::size() const{
    return this->pathEdges.size();
}

double PathEdges::getLength() const{
    return this->length;
}

int PathEdges::getElevation() const{
    return this->elevation;
}

Node* PathEdges::getStart(){
    return this->start;
}

Node* PathEdges::getEnd(){
    if (this->pathEdges.empty()){
        return this->start;
    }
    return this->pathEdges.back()->edge->v;
}

PathEdges::PathEdge* PathEdges::at(size_t index){
    if (index > this->pathEdges.size()){
        throw runtime_error("Segmentation Fault.");
    }
    return this->pathEdges.at(index);
}

Node::Edge* PathEdges::lastEdge(){
    if (this->pathEdges.empty())
        return nullptr;
    return this->pathEdges.back()->edge;
}

PathEdges* PathEdges::cutAfter(size_t index){
    return new PathEdges(this->start, this->pathEdges, index + 1);
}

PathEdges* PathEdges::randomCut(){
    int index = rand() % (this->size() + 1);
    return new PathEdges(this->start, this->pathEdges, index);
}

Path* PathEdges::toPath(){
    Path* path = new Path();
    for (auto it: this->pathEdges){
        path->edges.push_back(it->edge);
    }
    path->length = this->length;
    path->elevation = this->elevation;
    return path;
}


Path::Path(){
    this->length = 0.0;
    this->elevation = 0.0;
}

Path::Path(vector<Node::Edge*>& edges, double length, int elevation){
    this->edges = edges;
    this->length = length;
    this->elevation = elevation;
}

void Path::addEdge(Node::Edge* edge){
    this->edges.push_back(edge);
    this->length += edge->length;
    this->elevation += edge->elevation;
}

void Path::revert(){
    reverse(this->edges.begin(), this->edges.end());
}

Node* Path::getStart(){
    if (this->edges.empty()){
        return nullptr;
    }
    return this->edges[0]->u;
}

Node* Path::getEnd(){
    if (this->edges.empty()){
        return nullptr;
    }
    return this->edges.back()->v;
}