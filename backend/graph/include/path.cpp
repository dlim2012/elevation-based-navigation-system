#include <iostream>
#include <vector>
#include <algorithm>

#include "path.h"

using namespace std;


PathEdges::PathEdge::PathEdge(Node::Edge *edge, int length, int elevation) {
    this->edge = edge;
    this->length = length;
    this->elevation = elevation;
}

PathEdges::PathEdges(Node* start) {
    this->start = start;
    this->end = start;
    this->length = 0;
    this->elevation = 0;
}

PathEdges::PathEdges(Node* start, vector<Node::Edge *>& edges) {
    this->start = start;
    if (edges.empty()){
        this->end = start;
    } else {
        this->end = edges.back()->v;
    }
    this->length = 0;
    this->elevation = 0;
    for (Node::Edge *edge: edges) {
        this->addEdge(edge);
    }
}

PathEdges::PathEdges(Node* start, vector<PathEdge *> &pathEdges, int index) {
    this->start = start;
    this->pathEdges = vector<PathEdge*>(pathEdges.begin(), pathEdges.begin() + min((int) pathEdges.size(), index));
    this->end = this->pathEdges.empty() ? start : this->pathEdges.back()->edge->v;

    if (this->pathEdges.empty()) {
        this->length = 0;
        this->elevation = 0;
    } else {
        this->length = this->pathEdges.back()->length;
        this->elevation = this->pathEdges.back()->elevation;
    }
}

PathEdges::PathEdges(Path* path){
    this->start = this->end = path->getStart();
    for (Node::Edge* edge: path->edges){
        this->addEdge(edge);
    }
}

void PathEdges::addEdge(Node::Edge *edge) {
    if (edge->u != end){
        throw runtime_error("Start of edge does not match the end of path edges.");
    }
    this->length += edge->length;
    this->elevation += edge->elevation;
    this->pathEdges.push_back(new PathEdge(edge, this->length, this->elevation));
    this->end = edge->v;
    if (this->size() == 1){
        this->start = edge->u;
    }
}

bool PathEdges::empty() const{
    return this->pathEdges.empty();
}

int PathEdges::size() const{
    return this->pathEdges.size();
}

int PathEdges::getLength() const{
    return this->length;
}

int PathEdges::getElevation() const{
    return this->elevation;
}

Node* PathEdges::getStart(){
    return this->start;
}

Node* PathEdges::getEnd(){
    return this->end;
}

PathEdges::PathEdge* PathEdges::at(size_t index){
    if (index > this->pathEdges.size()){
        throw runtime_error("Segmentation Fault.");
    }
    return this->pathEdges.at(index);
}

Node::Edge* PathEdges::firstEdge(){
    if (this->pathEdges.empty())
        return nullptr;
    return this->pathEdges.front()->edge;
}

Node::Edge* PathEdges::lastEdge(){
    if (this->pathEdges.empty())
        return nullptr;
    return this->pathEdges.back()->edge;
}

PathEdges* PathEdges::cutBefore(int index){
    return new PathEdges(this->start, this->pathEdges, index);
}

PathEdges* PathEdges::cutAfter(int index){
    return new PathEdges(this->start, this->pathEdges, index+1);
}

PathEdges* PathEdges::randomCutReturnFront(){
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

Path::Path(vector<Node::Edge*>& edges, int length, int elevation){
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

bool Path::empty(){
    return this->edges.empty();
}