#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>

using namespace std;

enum direction_no {
    bidirectional = 0,
    outgoing = 1,
    incoming = 2,
};

struct Node {

    struct Edge {
        long id;
        double length;
        Node* target;
        enum direction_no direction;

        Edge(long id, double length, Node* target, enum direction_no valid){
            this->id = id;
            this->length=length;
            this->target = target;
            this->direction = valid;
        }
    };

    long id;
    double x;
    double y;
    vector<Edge> edges;

    Node(long id, double x, double y) {
        this->id = id;
        this->edges = vector<Edge>();
        this->x = x;
        this->y = y;
    }

    bool operator==(const Node* node2) const{
        return this->id == node2->id;
    }
};

struct BallTreeNode {
    Node* node;
    BallTreeNode* left = nullptr;
    BallTreeNode* right = nullptr;
    double radius = 0.0;
    BallTreeNode(Node* node){
        this->node = node;
    }
};


static double distance(Node* node1, Node* node2);
static double distance(Node* node1, double x, double y);


class Graph{
public:
    string version;
    // todo: make nodes and ballTree private
    unordered_map<long, Node*> nodes;
    BallTreeNode* ballTree;
private:
    string nodesCsvHeader;
    string edgesCsvHeader;

    void addNodes(const string& nodesFile);
    void addEdges(const string& edgesFile);
    void getGroup(Node* node, unordered_set<Node*>& visited, unordered_set<long>& groupNodes);
    unordered_set<long> getMaxGroup();
    void maxGroup();
    void writeToCsv();
    void createBallTree();
    BallTreeNode* createBallTree(vector<Node*> data);
    void nearestNode(BallTreeNode *ballTreeNode, pair<double, double> &target, Node *&node, double &minDistance);


public:
    explicit Graph(const string& version, const bool preprocessed);
    pair<double, vector<Node*>> shortestPath(long start_id, long end_id);
    Node* nearestNode(pair<double, double> &target);
    pair<double, vector<Node*>> findRandomPath(Node* start, Node* end, double maxLength, int seed);
};
