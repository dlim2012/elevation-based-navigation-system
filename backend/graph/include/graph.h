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


enum direction {
    bidirectional,
    oneway,
    onewayReversed,
    none
};


struct BallTreeNode {
    Node *node;
    BallTreeNode *left = nullptr;
    BallTreeNode *right = nullptr;
    double radius = 0.0;

    BallTreeNode(Node *node) {
        this->node = node;
    }
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
    bool populateEdgesField;
    int maxElevation = 0;

private:

    // add restrictions
    void addRestriction(
            unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions,
            Node::Edge *fromEdge,
            Node::Edge *toEdge);
    void addRestrictions(Node::Edge *fromEdge, Node::Edge *toEdge, string &type_enum);

    BallTreeNode *createBallTree(vector<Node *> data);
    void nearestNode(BallTreeNode *ballTreeNode, pair<double, double> &target, Node *&node, double &minDistance);

public:

    // initialize
    Graph(bool populateEdgesField);
    Graph(string &connectionString, MapConfig mapConfig);
    // todo: make graph from csv

    // add to graph
    Node *addNode(long id, double lon, double lat, int elevation);
    void addNode(Node *node);
    void addEdge(long id, Node *u, Node *v, double length, int elevation, direction direction);

    int addAllNodesFromDB(pqxx::connection_base &C, bool stream);
    int addAllEdgesFromDB(pqxx::connection_base &C, bool stream);
    int addAllRestrictionsFromDB(pqxx::connection_base &C);

    // save graph
    void saveAllEdgesInCsv(string connectionString, MapConfig mapConfig);


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
    static bool isRestricted(Node::Edge *p, Node::Edge *c, bool reversed);


};


void saveAllEdgesInCsv(string connectionString);