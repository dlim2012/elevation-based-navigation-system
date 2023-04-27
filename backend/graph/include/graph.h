#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>
#include <algorithm>
#include <pqxx/pqxx>
#include <stack>


#include "utils.h"
#include "path.h"

#ifndef QUERY_H
#define QUERY_H
#include "query.h"
#endif

using namespace std;




struct BallTreeNode {
    Node *node;
    BallTreeNode *left = nullptr;
    BallTreeNode *right = nullptr;
    int radius = 0; // unit: centimeter
    int count;

    BallTreeNode(Node *node, int count);
};

struct Geometry {
    long startId;
    long endId;
    vector<pair<double, double>> points;

    Geometry(long startId, long endId, const vector<pair<double, double>> &points);
};

class Graph {

public:
    // todo some fields can be private
    string version;
    unordered_map<long, Node *> nodes;
    BallTreeNode *ballTree = nullptr;
    default_random_engine rng;
    string nodesCsvHeader;
    string edgesCsvHeader;
    unordered_map<long, pair<Node::Edge *, Node::Edge *>> edges; // first edge: u->id < v->id
    int edgeCount;
    bool populateEdgesField;
    unordered_map<long, Geometry*> *edgeGeometries;

private:

    // add restrictions
    void addRestriction(
            unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions,
            Node::Edge *fromEdge,
            Node::Edge *toEdge);
    void addRestrictions(Node::Edge *fromEdge, Node::Edge *toEdge, string &type_enum);

    BallTreeNode *createBallTree(vector<Node *> data);
    void nearestNode(BallTreeNode *ballTreeNode, pair<double, double> &target, Node *&node, int &minDistance);

public:

    // initialize
    Graph(bool populateEdgesField);
    Graph(string& connectionString, HighwayConfig highwayConfig, LocationConfig locationConfig);
    Graph(string& connectionString,
          HighwayConfig highwayConfig, LocationConfig locationConfig,
          unordered_map<long, Geometry*>* sharedEdgesGeometries);
    // todo: make graph from csv

    // add to graph
    Node *addNode(long id, double lon, double lat, int elevation);
    void addNode(Node *node);
    void addEdge(long id, Node *u, Node *v, int length, int elevation, direction direction);

    int addAllNodesFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig, bool stream);
    int addAllEdgesFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig, bool stream);
    int addAllRestrictionsFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig);
    int addAllEdgeGeometriesFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig);

    // save graph
//    void saveAllEdgesInCsv();


    void edgeBasedStrongConnectIterative(
            Node::Edge *edge,
            unordered_map<Node::Edge *, pair<int, int>> &memo,
            int &index,
            unordered_map<long, pair<Node::Edge *, Node::Edge *>> &newEdges,
            unordered_map<long, Node *> &newNodes
    );

    void maxGroup();


    // nearest neighbor search
    BallTreeNode *createBallTree();
    Node *nearestNode(pair<double, double> target);

    // other functions
    Node* getRandomNode();
    pair<Node*, Node*> getTwoNearNodes(int maxDistance);
    static bool isRestricted(Node::Edge *p, Node::Edge *c, bool reversed);


};


void saveAllEdgesInCsv(string connectionString);