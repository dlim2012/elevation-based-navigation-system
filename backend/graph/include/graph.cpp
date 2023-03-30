#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <float.h>
#include <random>

#include "utils.h"
#include "graph.h"

using namespace std;


double distance(Node* node1, Node* node2){
    return sqrt(pow(node1->x - node2->x, 2) + pow(node1->y - node2->y, 2));
}

static double distance(Node* node1, pair<double, double>& p){
    return sqrt(pow(node1->x - p.first, 2) + pow(node1->y - p.second, 2));
}

void Graph::addNodes(const string& nodesFile){
    ifstream fin;
    string line, word;
    fin.open(nodesFile);
    if (fin.is_open()) {
        getline(fin, nodesCsvHeader);
        while (fin) {
            getline(fin, line);
            vector<string> words = parseCsvLine(line);
            if (words.size() < 4) {
                break;
            }
            long id = stol(words[1]);
            nodes[id] = new Node(id, stod(words[2]), stod(words[3]));
        }
    } else
        cout << "File " + nodesFile + " not opened" << endl;
    fin.close();;
}

void Graph::addEdges(const string& edgesFile){
    ifstream fin;
    string line, word;
    fin.open(edgesFile);
    if (fin.is_open()) {
        getline(fin, edgesCsvHeader);
        while (fin) {
            getline(fin, line);
            vector<string> words = parseCsvLine(line);
            if (words.size() < 5)
                break;
            int edge_id = stoi(words[0]);

            long u_id(stol(words[1])), v_id(stol(words[2]));
            double length(stod(words[3]));

            if (this->nodes.find(u_id) == this->nodes.end()){
                cout << "Node " << u_id << " not found";
                continue;
            }

            if (this->nodes.find(v_id) == this->nodes.end()){
                cout << "Node " << v_id << " not found";
                continue;
            }
            bool oneway = words[4] == "True";

            this->nodes[u_id]->edges.emplace_back(edge_id, length, this->nodes[v_id], bidirectional);
            this->nodes[v_id]->edges.emplace_back(edge_id, length, this->nodes[u_id], bidirectional);
//            this->nodes[u_id]->edges.emplace_back(edge_id, length, this->nodes[v_id], oneway ? outgoing : bidirectional);
//            this->nodes[v_id]->edges.emplace_back(edge_id, length, this->nodes[u_id], oneway ? incoming : bidirectional);
        }
    } else
        cout << "File " + edgesFile + " not opened" << endl;
    fin.close();
}

void Graph::getGroup(Node* node, unordered_set<Node*>& visited, unordered_set<long>& groupNodes){
    for (auto edge: node->edges){
        Node* neighbor = edge.target;
        if (visited.find(neighbor) == visited.end()){
            visited.insert(neighbor);
            groupNodes.insert(neighbor->id);
            getGroup(neighbor, visited, groupNodes);
        }
    }
}

unordered_set<long> Graph::getMaxGroup(){
    unordered_set<Node*> visited;
    unordered_set<long> group;
    size_t tot = 0;
    for (auto it: this->nodes){
        Node* curNode = it.second;
        unordered_set<long> groupNodes;
        if (visited.find(curNode) == visited.end()){
            visited.insert(curNode);
            groupNodes.insert(curNode->id);
            getGroup(curNode, visited, groupNodes);
            tot += groupNodes.size();
            if (groupNodes.size() > group.size()){
                swap(groupNodes, group);
            }
        }
    }

    return group;
}

void Graph::maxGroup(){
    unordered_set<long> group = getMaxGroup();

    unordered_map<long, Node*> newNodes;
    for (long node_id: group){
        Node* node = this->nodes[node_id];
        vector<Node::Edge> newEdges;
        for (Node::Edge edge: node->edges){
            if (group.find(edge.target->id) != group.end()){
                newEdges.push_back(edge);
            }
        }
        node->edges = newEdges;
        newNodes[node->id] = node;
    }
    this->nodes = newNodes;
}

void Graph::writeToCsv(){
    ofstream nodesFile;
    ofstream edgesFile;

    nodesFile.open(this->version + ".maxGroup.nodes.csv");
    edgesFile.open(this->version + ".maxGroup.edges.csv");

    nodesFile << nodesCsvHeader << endl;
    edgesFile << edgesCsvHeader << endl;

    nodesFile << std::setprecision(7);
    edgesFile << std::setprecision(3);

    int ni(0), ei(0);
    unordered_set<int> visitedEdges;
    for (auto it: this->nodes){
        Node* node = it.second;
        nodesFile << ni++ << "," << node->id << "," << node->x << "," << node->y << endl;
        for (Node::Edge edge: node->edges){
            if (edge.direction == bidirectional && visitedEdges.find(edge.id) == visitedEdges.end()) {
                visitedEdges.insert(edge.id);
                edgesFile << ei++ << "," << node->id << "," << edge.target->id << "," << edge.length << ","
                          << "False" << endl;
            } else if (edge.direction == outgoing){
                edgesFile << ei++ << "," << node->id << "," << edge.target->id << "," << edge.length << ","
                          << "True" << endl;
            }
        }
    }
    nodesFile.close();
    edgesFile.close();
}

Graph::Graph(const string& version, const bool preprocessed) {
    // Prerequisite: maxGroup nodes, edges to .csv file
    this->version = version;

    // Add all nodes
    addNodes(version + ".nodes.csv");

    // Add all edges (drop edges if node doesn't exist)
    addEdges(version + ".edges.csv");

    cout << nodes.size() << endl;

    // Get the biggest connected graph
    if (!preprocessed) {
        maxGroup();
        writeToCsv();
    }

    cout << "Creating BallTree" << endl;
    createBallTree();
    cout << "Graph Ready" << endl;
}

pair<double, vector<Node*>> Graph::shortestPath(long start_id, long end_id){

    // a-star algorithm

    struct candidate {
        double length;
        double heuristic;
        Node* node;
        candidate(double length, Node* cur, Node* target){
            this->length = length;
            this->heuristic = length + distance(cur, target);
            this->node = cur;
        }
    };

    struct comparator {
        bool operator() (const candidate& a, candidate& b){
            return a.heuristic > b.heuristic;
        }
    };

    vector<Node*> path;
    if (start_id == end_id){
        cout << 0.0 << endl;
        path.push_back(nodes[start_id]);
        return {0.0, path};
    }
    Node *start(this->nodes[start_id]), *end(this->nodes[end_id]);

    priority_queue<candidate, vector<candidate>, comparator> pq;

    unordered_map<long, double> minLength; // node_id, length
    pq.emplace(0.0, start, end);
    minLength[start->id] = 0.0;

    while (!pq.empty()){
        candidate c = pq.top();
        pq.pop();

        if (minLength.find(c.node->id) != minLength.end() && c.length > minLength[c.node->id] + 1e-10){
            continue;
        }

        for (Node::Edge edge: c.node->edges){
            if (edge.target->id == end->id) {
                double length = c.length + edge.length;
                double curLength = length;
                cout << "found! (length: " << curLength << ")" << endl;

                Node* cur = end;
                while (cur != start){
                    path.push_back(cur);
                    for (Node::Edge edge2: cur->edges){
                        auto it = minLength.find(edge2.target->id);

                        if (it != minLength.end() && abs(it->second + edge2.length - curLength) < 1e-10){
                            curLength = it->second;
                            cur = edge2.target;
                            break;
                        }
                    }
                }
                path.push_back(cur);
                reverse(path.begin(), path.end());
                return {length, path};
            }

            double length = c.length + edge.length;
            auto it = minLength.find(edge.target->id);
            if (it == minLength.end() || it->second > length) {
                minLength[edge.target->id] = length;
                pq.emplace(length, edge.target, end);
            }
        }
    }
    return {-1.0, path};
}

void Graph::createBallTree(){
    vector<Node*> data;
    for (auto it: nodes){
        data.push_back(it.second);
    }
    if (data.size() == 0){
        cout << "Graph is empty." << endl;
        return;
    }
    ballTree = createBallTree(data);
};

BallTreeNode* Graph::createBallTree(vector<Node*> data) {
    BallTreeNode* root;
    cout << data.size() << endl;
    if (data.size() == 1){
        root = new BallTreeNode(data[0]);
    } else {
        long x_mean(0.0), y_mean(0.0), x_spread(0.0), y_spread(0.0);
        for (Node* node: data){
            x_mean += node->x;
            y_mean += node->y;
        }
        x_mean = x_mean / data.size();
        y_mean = y_mean / data.size();
        for (Node* node: data){
            x_spread += pow(node->x - x_mean, 2);
            y_spread += pow(node->y - y_mean, 2);
        }
        int pivotIndex = data.size() / 2;
        Node* pivotNode = data[pivotIndex];
        if (x_spread > y_spread){
            sort(data.begin(), data.end(), [](const Node* a, Node* b){
                return a->x < b->x;
            });
        } else {
            sort(data.begin(), data.end(), [](const Node* a, Node* b){
                return a->y < b->y;
            });
        }
        root = new BallTreeNode(data[pivotIndex]);
        root->left = createBallTree(vector<Node*>(data.begin(), data.begin() + pivotIndex));
        if (pivotIndex + 1 < data.size())
            root->right = createBallTree(vector<Node*>(data.begin() + pivotIndex + 1, data.end()));


        for (Node* node: data){
            root->radius = max(root->radius, distance(node, pivotNode));
        }
    }
    return root;
}

void Graph::nearestNode(BallTreeNode* ballTreeNode, pair<double, double>& target, Node* &node, double &minDistance) {
    if (ballTreeNode == nullptr)
        return;
    if (distance(ballTreeNode->node, target) - ballTreeNode->radius >= minDistance)
        return;

    double curDistance = distance(ballTreeNode->node, target);
    if (curDistance < minDistance){
        node = ballTreeNode->node;
        minDistance = curDistance;
    }
    nearestNode(ballTreeNode->left, target, node, minDistance);
    nearestNode(ballTreeNode->right, target, node, minDistance);
}

Node* Graph::nearestNode(pair<double, double> &target) {

    Node* res = ballTree->node;
    double minDistance = DBL_MAX;
    nearestNode(ballTree, target, res, minDistance);
    return res;
}

pair<double, vector<Node*>> Graph::findRandomPath(Node *start, Node *end, double maxLength, int seed) {
    srand(seed);
    vector<Node*> path;

    while (true) {
        path.clear();
        unordered_set<long> visited;
        double length = 0.0;

        Node* node = start;
        path.push_back(start);
        visited.insert(start->id);

        while (node->id != end->id && length + distance(node, end) < maxLength) {
            vector<Node::Edge> validEdges;
            for (Node::Edge edge: node->edges){
                if (visited.find(edge.target->id) == visited.end()){
                    validEdges.push_back(edge);
                }
            }
            if (validEdges.empty()){
                break;
            }

            int randomIndex = rand() % validEdges.size();
            Node::Edge edge = validEdges[randomIndex];

            node = edge.target;
            path.push_back(node);
            visited.insert(node->id);
            length += edge.length;
        }
        if (node->id == end->id)
            return {length, path};
    }
}