
#include <math.h>
#include <unordered_map>
#include <algorithm>

#include "node.h"

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
                 double length, int elevation) {
    this->id = id;
    this->u = u;
    this->v = v;
    this->length = length;
    this->elevation = elevation;
}


double Node::Edge::getLength() { return this->length; };

double Node::Edge::getElevation() { return this->elevation; };

double Node::Edge::getLength(Edge *edge) { return edge->getLength(); }

double Node::Edge::getElevation(Edge *edge) { return edge->getElevation(); }


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


double distance(double lon1, double lon2, double lat1, double lat2) {
    double lon1_rad = lon1 * M_PI / 180;
    double lon2_rad = lon2 * M_PI / 180;
    double lat1_rad = lat1 * M_PI / 180;
    double lat2_rad = lat2 * M_PI / 180;
    return acos(sin(lat1_rad) * sin(lat2_rad) + cos(lat1_rad) * cos(lat2_rad) * cos(lon2_rad - lon1_rad)) * 6371 * 1000;
}

double distance(Node *node1, Node *node2) {
    return distance(node1->lon, node2->lon, node1->lat, node2->lat);
    //    return sqrt(pow(node1->lon - node2->lon, 2) + pow(node1->lat - node2->lat, 2));
}

double distance(Node *node1, pair<double, double> &p) {
    return distance(node1->lon, p.first, node1->lat, p.second);
//    return sqrt(pow(node1->lon - p.first, 2) + pow(node1->lat - p.second, 2));
}