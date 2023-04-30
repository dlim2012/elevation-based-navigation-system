
#include <math.h>
#include <unordered_map>
#include <algorithm>
#include <iostream>

#include "node.h"

#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

using namespace std;

Node::Edge::Edge(Node *v) {
    this->id = -1;
    this->u = nullptr;
    this->v = v;
    this->length = -1.0;
    this->elevation = -1;
}

Node::Edge::Edge(Node *node, bool reversed) {
    this->id = -1;
    this->u = reversed ? node : nullptr;
    this->v = reversed ? nullptr : node;
    this->length = -1.0;
    this->elevation = -1;
}

Node::Edge::Edge(long id, Node *u, Node *v,
                 int length, int elevation) {
    this->id = id;
    this->u = u;
    this->v = v;
    this->length = length;
    this->elevation = elevation;
}

Node::Edge::~Edge(){
}

long Node::Edge::getId() const { return this->id; }

int Node::Edge::getLength() { return this->length; };

int Node::Edge::getElevation() { return this->elevation; };

unordered_map<long, Node::Edge *>& Node::getEdges(){
    return this->edges;
}

unordered_map<long, Node::Edge *>& Node::getReversedEdges(){
    return this->reversedEdges;
}

bool Node::hasEdge(long id) {
    return this->edges.find(id) != this->edges.end();
}

bool Node::hasReversedEdge(long id) {
    return this->reversedEdges.find(id) != this->reversedEdges.end();
}

void Node::setEdge(long id, Node::Edge* edge){
    this->edges[id] = edge;
}

void Node::setReversedEdge(long id, Node::Edge* edge){
    this->reversedEdges[id] = edge;
}

unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *Node::getRestrictions(){
    return this->restrictions;
}

unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *Node::getReversedRestrictions(){
    return this->reversedRestrictions;
}

void Node::setRestrictions(unordered_map<Node::Edge *, unordered_set<Node::Edge *>> * restrictions){
    this->restrictions = restrictions;
}

void Node::setReversedRestrictions(unordered_map<Node::Edge *, unordered_set<Node::Edge *>> * reversedRestrictions){
    this->reversedRestrictions = reversedRestrictions;
}

Node* Node::Edge::getU() { return this->u; }

Node* Node::Edge::getV() { return this->v; }

int Node::Edge::getLength(Edge *edge) { return edge->getLength(); }

int Node::Edge::getElevation(Edge *edge) { return edge->getElevation(); }

Node::Node(long id, double lon, double lat, int elevation) {
    this->id = id;
    this->lon = lon;
    this->lat = lat;
    this->elevation = elevation;
    this->restrictions = nullptr;
    this->reversedRestrictions = nullptr;
}

Node::Node(long id, double lon, double lat, int elevation, unordered_map<long, Edge *> &edges) {
    this->id = id;
    this->lon = lon;
    this->lat = lat;
    this->elevation = elevation;
    this->edges = edges;
    this->restrictions = nullptr;
    this->reversedRestrictions = nullptr;
}

Node::~Node(){
    delete this->restrictions;
    delete this->reversedRestrictions;
    for (pair<const long, Node::Edge*> pair1: this->edges){
        pair1.second->getV()->getReversedEdges().erase(pair1.second->getId());
        delete pair1.second;
    }
    for (pair<const long, Node::Edge*> pair1: this->reversedEdges){
        pair1.second->getU()->getEdges().erase(pair1.second->getId());
        delete pair1.second;
    }
}

Node *Node::cloneWithoutEdges() {
    return new Node(this->id, this->lon, this->lat, this->elevation);
}

//void Node::addEdge(long id, Node* v, double length, double elevation, enum direction_no direction){
//    this->edges[id] = new Edge(id, this, v, length, elevation, direction);
//}
//
//void Node::addClonedEdge(Edge* edge, Node* newV){
//    this->edges[edge->id] = new Edge(edge->id, this, newV, edge->length, edge->elevation, edge->direction);
//}

bool Node::operator==(const Node *node2) const {
    return this->id == node2->id;
}


int distance(double lon1, double lon2, double lat1, double lat2) {
    double lon1_rad = lon1 * M_PI / 180;
    double lon2_rad = lon2 * M_PI / 180;
    double lat1_rad = lat1 * M_PI / 180;
    double lat2_rad = lat2 * M_PI / 180;
    return roundToInt(acos(sin(lat1_rad) * sin(lat2_rad) + cos(lat1_rad) * cos(lat2_rad) * cos(lon2_rad - lon1_rad)) * 6371 * 1000 * 100);
}

int distance(Node *node1, Node *node2) {
    return distance(node1->getLon(), node2->getLon(), node1->getLat(), node2->getLat());
    //    return sqrt(pow(node1->lon - node2->lon, 2) + pow(node1->lat - node2->lat, 2));
}

int distance(Node *node1, pair<double, double> &p) {
    return distance(node1->getLon(), p.first, node1->getLat(), p.second);
//    return sqrt(pow(node1->lon - p.first, 2) + pow(node1->lat - p.second, 2));
}

long Node::getId(){
    return this->id;
}

double Node::getLon() const{
    return this->lon;
}

double Node::getLat() const{
    return this->lat;
}

int Node::getElevation() const{
    return this->elevation;
}