#include <iostream>
#include <vector>
#include <algorithm>

#include "path.h"

using namespace std;


PathEdges::PathEdge::PathEdge(Node::Edge *edge, int length, int elevation) {
    this->edge = edge;
    this->length = length;
    this->elevation = elevation;
    this->referCount = 1;
}

PathEdges::PathEdge::~PathEdge(){
}

Node::Edge* PathEdges::PathEdge::getEdge(){
    return this->edge;
}

int PathEdges::PathEdge::getLength(){
    return this->length;
}

int PathEdges::PathEdge::getElevation(){
    return this->elevation;
}

void PathEdges::PathEdge::increaseReferCount(){
    this->referCount++;
}

bool PathEdges::PathEdge::decreaseReferCount(){
    this->referCount--;
    return this->referCount == 0;
}

PathEdges::PathEdge *PathEdges::PathEdge::clone(){
    return new PathEdges::PathEdge(this->edge, this->length, this->elevation);
}

PathEdges::PathEdges(Node* start) {
    this->start = start;
    this->end = start;
    this->length = 0;
    this->elevation = 0;
    this->confirmedlyOptimal = false;
}

PathEdges::PathEdges(Node* start, vector<Node::Edge *>& edges) {
    this->start = start;
    if (edges.empty()){
        this->end = start;
    } else {
        this->end = edges.back()->getV();
    }
    this->length = 0;
    this->elevation = 0;
    for (Node::Edge *edge: edges) {
        this->addEdge(edge);
    }
    this->confirmedlyOptimal = false;
}

PathEdges::PathEdges(Node* start, vector<PathEdge *> &pathEdges, int index) {
    this->start = start;
//    this->pathEdges = vector<PathEdge*>(pathEdges.begin(), pathEdges.begin() + min((int) pathEdges.size(), index));
    this->pathEdges = vector<PathEdge*>();
    index = (int) min((int) pathEdges.size(), index);
    for (int i=0; i<index; i++){
        this->pathEdges.push_back(pathEdges[i]->clone());
    }
    this->end = this->pathEdges.empty() ? start : this->pathEdges.back()->getEdge()->getV();


    if (this->pathEdges.empty()) {
        this->length = 0;
        this->elevation = 0;
    } else {
        this->length = this->pathEdges.back()->getLength();
        this->elevation = this->pathEdges.back()->getElevation();
    }

//    for(PathEdges::PathEdge *pathEdge: this->pathEdges){
//        pathEdge->increaseReferCount();
//    }
    this->confirmedlyOptimal = false;
}

PathEdges::PathEdges(Path* path){
    this->start = this->end = path->getStart();
    this->length = 0;
    this->elevation = 0;
    for (Node::Edge* edge: path->getEdges()){
        this->addEdge(edge);
    }
}

PathEdges::~PathEdges(){
    for (PathEdges::PathEdge * pathEdge: this->pathEdges){
        if (pathEdge->decreaseReferCount()) {
            delete pathEdge;
        }
    }
}

void PathEdges::addEdge(Node::Edge *edge) {
    if (edge->getU() != end){
        throw runtime_error("Start of edge does not match the end of path edges.");
    }
    this->length += edge->getLength();
    this->elevation += edge->getElevation();
    this->pathEdges.push_back(new PathEdge(edge, this->length, this->elevation));
    this->end = edge->getV();
    if (this->size() == 1){
        this->start = edge->getU();
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

vector<PathEdges::PathEdge *>& PathEdges::getPathEdges(){
    return this->pathEdges;
};

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

Node::Edge* PathEdges::getFirstEdge(){
    if (this->pathEdges.empty())
        return nullptr;
    return this->pathEdges.front()->getEdge();
}

Node::Edge* PathEdges::getLastEdge(){
    if (this->pathEdges.empty())
        return nullptr;
    return this->pathEdges.back()->getEdge();
}

PathEdges* PathEdges::cutBefore(int index){
    return new PathEdges(this->start, this->pathEdges, index);
}

bool PathEdges::isConfirmedOptimal() {
    return this->confirmedlyOptimal;
}

void PathEdges::confirmOptimal() {
    this->confirmedlyOptimal = true;
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
        path->getEdges().push_back(it->getEdge());
    }
    path->setLength(this->length);
    path->setElevation(this->elevation);

    return path;
}


Path::Path(){
    this->length = 0.0;
    this->elevation = 0.0;
    this->confirmedlyOptimal = false;
}

Path::~Path(){}

Path::Path(vector<Node::Edge*>& edges, int length, int elevation){
    this->edges = edges;
    this->length = length;
    this->elevation = elevation;
    this->confirmedlyOptimal = false;
}

void Path::addEdge(Node::Edge* edge){
    this->edges.push_back(edge);
    this->length += edge->getLength();
    this->elevation += edge->getElevation();
}

void Path::revert(){
    reverse(this->edges.begin(), this->edges.end());
}

Node* Path::getStart(){
    if (this->edges.empty()){
        return nullptr;
    }
    return this->edges[0]->getU();
}

Node* Path::getEnd(){
    if (this->edges.empty()){
        return nullptr;
    }
    return this->edges.back()->getV();
}


vector<Node::Edge*>& Path::getEdges(){
    return this->edges;
}

int Path::getLength(){
    return this->length;
}

int Path::getElevation(){
    return this->elevation;
}

bool Path::isConfirmedOptimal() {
    return this->confirmedlyOptimal;
}

void Path::setLength(int length){
    this->length = length;
}

void Path::setElevation(int elevation){
    this->elevation = elevation;
}

void Path::confirmOptimal(){
    this->confirmedlyOptimal = true;
}

bool Path::empty() const{
    return this->edges.empty();
}

size_t pathEdgesHash::operator() (PathEdges *pathEdges) const {
    size_t hash = 0;
    for (PathEdges::PathEdge* pathEdge: pathEdges->getPathEdges()){
        hash = (hash * HASHMULT + ((int) (pathEdge->getEdge()->getV()->getId() % HASHMOD))) % HASHMOD;
    }
    return hash;
}

bool pathEdgesEqual::operator() (PathEdges *pathEdges1, PathEdges* pathEdges2) const {
    if (pathEdges1->size() != pathEdges2->size()){
        return false;
    }
    size_t size = pathEdges1->size();
    for (size_t i=0; i<size; i++){
        if (pathEdges1->getPathEdges()[i]->getEdge() != pathEdges2->getPathEdges()[i]->getEdge()){
            return false;
        }
    }
    return true;
}