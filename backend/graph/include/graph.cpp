#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <float.h>
#include <random>
#include <algorithm>
#include <pqxx/pqxx>
#include <chrono>

#include "graph.h"



using namespace std;
using namespace pqxx;

BallTreeNode::BallTreeNode(Node *node, int count) {
    this->node = node;
    this->count = count;
}

Geometry::Geometry(long startId, long endId, const vector<pair<double, double>> &points){
    this->startId = startId;
    this->endId = endId;
    this->points = points;
}

Graph::Graph(bool populateEdgesField) {
    this->populateEdgesField = populateEdgesField;
    srand(time(nullptr));
    this->edgeCount = 0;
}

Graph::Graph(string& connectionString,
      HighwayConfig highwayConfig, LocationConfig locationConfig,
      unordered_map<long, Geometry*>* sharedEdgesGeometries){
    this->edgeGeometries = sharedEdgesGeometries;
    this->populateEdgesField = true;
    connection C(connectionString);
    addAllNodesFromDB(C, highwayConfig, locationConfig, false);
    addAllEdgesFromDB(C, highwayConfig, locationConfig, false);
    addAllRestrictionsFromDB(C, highwayConfig, locationConfig);
    maxGroup();
    addAllEdgeGeometriesFromDB(C, highwayConfig, locationConfig);
    this->ballTree = createBallTree();
    this->populateEdgesField = false;
    this->edges.clear();
}

void deleteBallTree(BallTreeNode* node){
    if (node->left != nullptr){
        deleteBallTree(node->left);
    }
    if (node->right != nullptr){
        deleteBallTree(node->right);
    }
    delete node->left;
    delete node->right;
}

Graph::~Graph(){
    cout << "delete Graph" << endl;
    if (this->ballTree != nullptr) {
        deleteBallTree(this->ballTree);
    }
    delete this->ballTree;
    for (pair<long, Node*> p: nodes){
        delete p.second;
    }
}

//void Graph::saveAllEdgesInCsv() {
//
////    Graph *graph = new Graph(C, highwayConfig);
//
//    ofstream edgesFile;
//    edgesFile.open("csv_files/mass-edges.csv");
//
//    edgesFile << fixed << setprecision(8) << "index,longitude1,latitude1,longitude2,latitude2" << endl;
//
//    int index(0);
//    for (auto pair1: this->getEdges()) {
//        vector<pair<double, double>> &v = this->edgeGeometries[pair1.first];
//        for (int i = 0; i < v.size() - 1; i++) {
//            edgesFile << index++ << "," << v[i].first << "," << v[i].second << "," << v[i + 1].first << ","
//                      << v[i + 1].second << endl;
//        }
//    }
//    edgesFile.close();
//}


Node *Graph::addNode(long id, double lon, double lat, int elevation) {
    return this->nodes[id] = new Node(id, lon, lat, elevation);
}

void Graph::addNode(Node *node) {
    this->nodes[node->getId()] = node;
}



void Graph::addEdge(long id, Node *u, Node *v, int length, int elevation, direction direction) {
    if (direction == none) {
        return;
    }

    Node::Edge *edge1, *edge2;

    if (direction == bidirectional || direction == oneway) {
        if (!u->hasEdge(id)){
            this->edgeCount++;
        }
        edge1 = new Node::Edge(id, u, v, length, elevation);
        u->setEdge(id, edge1);
        v->setReversedEdge(id, edge1);
    } else {
        edge1 = nullptr;
    }

    if (u != v && (direction == bidirectional || direction == onewayReversed)) {
        if (!v->hasEdge(id)) {
            this->edgeCount++;
        }
        edge2 = new Node::Edge(id, v, u, length, elevation);
        v->setEdge(id, edge2);
        u->setReversedEdge(id, edge2);
    } else {
        edge2 = nullptr;
    }

    if (this->populateEdgesField) {
        if (u->getId() < v->getId()) {
            this->edges[id] = {edge1, edge2};
        } else {
            this->edges[id] = {edge2, edge1};
        }
    }
}

int Graph::addAllNodesFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig, bool stream) {

    auto t0 = chrono::high_resolution_clock::now();
    cout << "Fetching nodes data from DB... ";
    prepare_nodes_query(C, highwayConfig, locationConfig);
    pqxx::work w(C);
    if (stream) {
        throw runtime_error("Not Implemented.");
    } else {
        pqxx::result r = w.exec_prepared("nodes");
        w.commit();
        for (pqxx::row row1: r) {
            long id = stol(row1[0].c_str());
            this->nodes[id] = new Node(
                    id,
                    stod(row1[1].c_str()),
                    stod(row1[2].c_str()),
                    stoi(row1[3].c_str()));
        }
        auto t1 = chrono::high_resolution_clock::now();
        auto t_s = (double) chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000;
        cout << "done (count: " << r.size() << ", " << t_s << " seconds)" << endl;
        return r.size();
    }
}

int Graph::addAllEdgesFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig, bool stream) {

    auto t0 = chrono::high_resolution_clock::now();
    cout << "Fetching edges data from DB... ";

    prepare_edges_query(C, highwayConfig, locationConfig);
    pqxx::work w(C);
    if (stream) {
        throw runtime_error("Not Implemented.");
    } else {
        pqxx::result r = w.exec_prepared("edges");
        w.commit();
        for (pqxx::row row1: r) {
            long id = stol(row1[0].c_str());
            int length = (int) (stod(row1[3].c_str()) * 100);
            int elevation = stoi(row1[4].c_str());

            long u_id = stol(row1[1].c_str());
            long v_id = stol(row1[2].c_str());

            // some edges have same u and v: decided to include them

            auto itU = this->nodes.find(u_id);
            if (itU == this->nodes.end())
                continue;
            auto itV = this->nodes.find(v_id);
            if (itV == this->nodes.end())
                continue;

            direction direction;
            if (stoi(row1[5].c_str()) == 0 || highwayConfig == cycling || highwayConfig == hiking) {
                direction = bidirectional;
            } else {
                direction = oneway;
            }
            this->addEdge(id, itU->second, itV->second, length, elevation, direction);
        }

        int count = 0;
        for (auto p: this->edges) {
            if (p.second.first != nullptr)
                count++;
            if (p.second.second != nullptr) {
                count++;
            }
        }
        this->edgeCount = count;

        auto t1 = chrono::high_resolution_clock::now();
        auto t_s = (double) chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000;
        cout << "done (count: " << r.size() << ", directed edges: " << count << ", " << t_s << " seconds)" << endl;

        return r.size();
    }
}

void Graph::addRestriction(unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions, Node::Edge *fromEdge,
                           Node::Edge *toEdge) {
    // make an element for the 'from' edge in restrictions if not exists
    unordered_set<Node::Edge *> *restricted;
    auto restrictionsIt = restrictions->find(fromEdge);
    if (restrictionsIt == restrictions->end()) {
        restrictions->insert({{fromEdge, unordered_set<Node::Edge *>()}});
        restricted = &(restrictions->at(fromEdge));
    } else {
        restricted = &(restrictionsIt->second);
    }

    restricted->insert(toEdge);
}

void Graph::addRestrictions(Node::Edge *fromEdge, Node::Edge *toEdge, string &type_enum) {
    Node *viaNode = fromEdge->getV();
    // add the restrictions map to 'via' if not exists
    if (viaNode->getRestrictions() == nullptr) {
        viaNode->setRestrictions(new unordered_map<Node::Edge *, unordered_set<Node::Edge *>>());
        viaNode->setReversedRestrictions(new unordered_map<Node::Edge *, unordered_set<Node::Edge *>>());
    }
    unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions = viaNode->getRestrictions();
    unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *reversedRestrictions = viaNode->getReversedRestrictions();

    if (type_enum[1] == 'n') {
        addRestriction(restrictions, fromEdge, toEdge);
        addRestriction(reversedRestrictions, toEdge, fromEdge);
    } else {
        for (pair<const long, Node::Edge *> &pair1: viaNode->getEdges()) {
            if (pair1.second != toEdge) {
                addRestriction(restrictions, fromEdge, pair1.second);
                addRestriction(reversedRestrictions, pair1.second, fromEdge);
            }
        }
    }
}

int Graph::addAllRestrictionsFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig) {
    // Note: the populateEdgesField has to be true to add restrictions to edges
    if (highwayConfig == cycling || highwayConfig == hiking || highwayConfig == motorway_without_restrictions){
        return 0;
    }

    auto t0 = chrono::high_resolution_clock::now();
    cout << "Fetching restrictions data from DB... ";

    prepare_restrictions_query(C, highwayConfig, locationConfig);
    if (!this->populateEdgesField) {
        throw runtime_error("The populateEdgesField has to be true to add restrictions to edges.");
    }

    pqxx::work w(C);
    pqxx::result r = w.exec_prepared("restrictions");
    w.commit();

    for (pqxx::row row1: r) {
        string type_enum = row1[0].c_str();
        long from_id = stol(row1[1].c_str());
        long via_id = stol(row1[2].c_str());
        long to_id = stol(row1[3].c_str());
        Node::Edge *from(nullptr), *to(nullptr);


        if (type_enum[0] == 'n') {

            // find the 'via' node (skip if not found)
            auto nodesIt = this->nodes.find(via_id);
            if (nodesIt == this->nodes.end())
                continue;
            Node *via = nodesIt->second;

            // find the 'from' edge (skip if not found or invalid)
            auto edgesIt = this->edges.find(from_id);
            if (edgesIt == this->edges.end())
                continue;
            pair<Node::Edge *, Node::Edge *> edgePair = edgesIt->second;
            if (edgePair.first != nullptr && edgePair.first->getV() == via) {
                from = edgePair.first;
            } else if (edgePair.second != nullptr && edgePair.second->getV() == via) {
                from = edgePair.second;
            }
            if (from == nullptr)
                continue;

            // find the 'to' edge (skip if not found)
            auto nodeEdgesIt = via->getEdges().find(to_id);
            if (nodeEdgesIt == via->getEdges().end())
                continue;
            to = nodeEdgesIt->second;

            addRestrictions(from, to, type_enum);
        } else if (type_enum[0] == 'e') {

            Node::Edge *from, *via(nullptr), *to;

            // find the two candidates for the 'via' edge that matches the edge id
            auto edgesIt = this->edges.find(via_id);
            if (edgesIt == this->edges.end())
                continue;
            pair<Node::Edge *, Node::Edge *> p = edgesIt->second;

            // find the 'via' edge and the 'to' edge
            if (p.first == nullptr) {
                via = p.second;
                auto nodeEdgesIt = via->getV()->getEdges().find(to_id);
                if (nodeEdgesIt == via->getV()->getEdges().end())
                    continue;
                to = nodeEdgesIt->second;
            } else if (p.second == nullptr) {
                via = p.first;
                auto nodeEdgesIt = via->getV()->getEdges().find(to_id);
                if (nodeEdgesIt == via->getV()->getEdges().end())
                    continue;
            } else {
                auto nodeEdgesIt = p.first->getV()->getEdges().find(to_id);
                if (nodeEdgesIt != p.first->getV()->getEdges().end()) {
                    to = nodeEdgesIt->second;
                    via = p.first;
                } else {
                    nodeEdgesIt = p.second->getV()->getEdges().find(to_id);
                    if (nodeEdgesIt != p.second->getV()->getEdges().end()) {
                        to = nodeEdgesIt->second;
                        via = p.second;
                    } else {
                        continue;
                    }
                }
            }

            // find the 'from' edge
            edgesIt = this->edges.find(from_id);
            if (edgesIt == this->edges.end())
                continue;
            p = edgesIt->second;
            if (p.first != nullptr && p.first->getV() == via->getU()) {
                from = p.first;
            } else if (p.second != nullptr && p.second->getV() == via->getU()) {
                from = p.second;
            } else
                continue;

            addRestrictions(from, via, type_enum);
            addRestrictions(via, to, type_enum);

        }

    }
    auto t1 = chrono::high_resolution_clock::now();
    auto t_s = (double) chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000;
    cout << "done (count: " << r.size() << ", " << t_s << " seconds)" << endl;

    return r.size();
}



struct StackElement {
    Node::Edge *parEdge;
    Node::Edge *childEdge;
    unordered_map<long, Node::Edge *>::iterator itCur;

    StackElement(Node::Edge *parEdge, Node::Edge *childEdge, unordered_map<long, Node::Edge *>::iterator itCur) {
        this->parEdge = parEdge;
        this->childEdge = childEdge;
        this->itCur = itCur;
    }
};

void Graph::edgeBasedStrongConnectIterative(
        Node::Edge *edge,
        unordered_map<Node::Edge *, pair<int, int>> &memo,
        int &index,
        unordered_map<long, pair<Node::Edge *, Node::Edge *>> &newEdges,
        unordered_map<long, Node *> &newNodes
) {

    stack<Node::Edge *> tarjanStack;
    unordered_set<Node::Edge *> onTarjanStack;
    stack<StackElement> stack;

    tarjanStack.push(edge);
    onTarjanStack.insert(edge);
    memo[edge] = {index, index};
    index++;
    stack.push({edge, nullptr, edge->getV()->getEdges().begin()});

    while (!stack.empty()) {
        Node::Edge *parEdge = stack.top().parEdge;
        Node::Edge *childEdge = stack.top().childEdge;
        unordered_map<long, Node::Edge *>::iterator itCur = stack.top().itCur;
        stack.pop();

        pair<int, int> &parMemo = memo[parEdge];
        if (childEdge != nullptr) {
            parMemo.second = min(parMemo.second, memo[childEdge].second);
        }

        while (itCur != parEdge->getV()->getEdges().end()) {
            childEdge = (itCur)->second;

            if (isRestricted(parEdge, childEdge, false)) {
                itCur++;
                continue;
            }

            auto memoIt = memo.find(childEdge);
            if (memoIt == memo.end()) {
                break;
            }
            if (onTarjanStack.find(childEdge) != onTarjanStack.end())
                parMemo.second = min(parMemo.second, memoIt->second.first);
            itCur++;
        }
        if (itCur != parEdge->getV()->getEdges().end()) {
            tarjanStack.push(childEdge);
            onTarjanStack.insert(childEdge);
            memo[childEdge] = {index, index};
            index++;
            stack.push({parEdge, childEdge, ++itCur});
            stack.push({childEdge, nullptr, childEdge->getV()->getEdges().begin()});
        } else if (parMemo.first == parMemo.second) {
            unordered_map<long, pair<Node::Edge *, Node::Edge *>> edgeGroup;
            unordered_map<long, Node *> nodeGroup;

            while (true) {
                edge = tarjanStack.top();
                tarjanStack.pop();
                onTarjanStack.erase(edge);

                auto edgePairIt = edgeGroup.find(edge->getId());
                if (edgePairIt == edgeGroup.end()) {

                    if (edge->getU()->getId() < edge->getV()->getId()) {
                        edgeGroup[edge->getId()] = {edge, nullptr};
                    } else {
                        edgeGroup[edge->getId()] = {nullptr, edge};
                    }
                } else {
                    if (edge->getU()->getId() < edge->getV()->getId()) {
                        edgePairIt->second.first = edge;
                    } else {
                        edgePairIt->second.second = edge;
                    }
                }

                if (nodeGroup.find(edge->getU()->getId()) == nodeGroup.end()) {
                    nodeGroup[edge->getU()->getId()] = edge->getU();
                }
                if (nodeGroup.find(edge->getV()->getId()) == nodeGroup.end()) {
                    nodeGroup[edge->getV()->getId()] = edge->getV();
                }

                if (edge == parEdge) {
                    break;
                }
            }

            if (nodeGroup.size() > newNodes.size()) {
                swap(nodeGroup, newNodes);
                swap(edgeGroup, newEdges);
            }

        }
    }


}


void Graph::maxGroup() {
    auto t0 = chrono::high_resolution_clock::now();
    cout << "Extracting the biggest strongly connected components... ";
    int index(0);
    stack<Node::Edge *> tarjanStack;
    unordered_set<Node::Edge *> onTarjanStack;
    unordered_map<Node::Edge *, pair<int, int>> memo;
    unordered_map<long, pair<Node::Edge *, Node::Edge *>> newEdges;
    unordered_map<long, Node *> newNodes;
    stack<Node::Edge *> stack;

    // use Tarjan's algorithm to get MaxGroup and save the new Nodes;
    for (pair<const long, pair<Node::Edge *, Node::Edge *>> &p: this->edges) {
        if (p.second.first != nullptr && memo.find(p.second.first) == memo.end())
            edgeBasedStrongConnectIterative(p.second.first, memo, index, newEdges, newNodes);
        if (p.second.second != nullptr && memo.find(p.second.second) == memo.end())
            edgeBasedStrongConnectIterative(p.second.second, memo, index, newEdges, newNodes);
    }

    for (auto it: this->nodes){
        Node* node = it.second;
        if (newNodes.find(it.first) == newNodes.end()){
            delete node;
        }
    }

    int edgeCount = 0;
    for (pair<const long, Node *> &pair1: newNodes) {
        Node *node = pair1.second;
        for (auto edgeIt = node->getEdges().begin(); edgeIt != node->getEdges().end();edgeIt++) {
            edgeCount++;
        }
    }

    int countNull(0);
    for (auto it: newEdges) {
        if (it.second.first == nullptr && it.second.second == nullptr)
            countNull++;
    }

    swap(this->nodes, newNodes);
    swap(this->edges, newEdges);
    this->edgeCount = edgeCount;

    auto t1 = chrono::high_resolution_clock::now();
    auto t_s = (double) chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000;
    cout << "done (node: " << this->nodes.size() << ", edge: "
         << edgeCount << ", " << t_s << " seconds)" << endl;


}

int Graph::addAllEdgeGeometriesFromDB(pqxx::connection_base &C, HighwayConfig highwayConfig, LocationConfig locationConfig){

    auto t0 = chrono::high_resolution_clock::now();
    cout << "Fetching edge geometries... ";

    prepare_edges_geom_query(C, highwayConfig, locationConfig);
    pqxx::work w(C);
    pqxx::result r = w.exec_prepared("edges_geom");
    long lineCount = 0;
    for (pqxx::row row: r) {
        long id = stol(row[0].c_str());
        string geom_text = row[1].c_str();
        long u = stol(row[2].c_str());
        long v = stol(row[3].c_str());
        if (this->edges.find(id) == this->edges.end()){
            continue;
        }
        (*this->edgeGeometries)[id] = new Geometry(u, v, parseLineString(geom_text));
        lineCount += this->edgeGeometries->at(id)->points.size();
    }

    auto t1 = chrono::high_resolution_clock::now();
    auto t_s = (double) chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000;
    cout << "done (row count: " << r.size() << ", line count: " << lineCount << ", " << t_s << " seconds)" << endl;
    return r.size();
}

bool Graph::isRestricted(Node::Edge *p, Node::Edge *c, bool reversed) {
    if (p == nullptr){
        return false;
    }
    unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions =
            reversed ? p->getU()->getReversedRestrictions() : p->getV()->getRestrictions();

    if (restrictions != nullptr) {
        auto restrictionsIt = restrictions->find(p);
        if (restrictionsIt != restrictions->end()) {
            if (restrictionsIt->second.find(c) != restrictionsIt->second.end()) {
                return true;
            }
        }
    }
    return false;
}

BallTreeNode *Graph::createBallTree() {
    auto t0 = chrono::high_resolution_clock::now();
    cout << "Creating ball tree... ";
    vector<Node *> data;
    for (auto it: nodes) {
        data.push_back(it.second);
    }
    if (data.empty()) {
        cout << "Graph is empty." << endl;
        return nullptr;
    }
    this->ballTree = createBallTree(data);
    auto t1 = chrono::high_resolution_clock::now();
    auto t_s = (double) chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000;
    cout << "done (" << t_s << " seconds)" << endl;

    return this->ballTree;
}

BallTreeNode *Graph::createBallTree(vector<Node *> data) {
    BallTreeNode *root;
    if (data.size() == 1) {
        root = new BallTreeNode(data[0], 1);
    } else {
        double x_mean(0.0), y_mean(0.0), x_spread(0.0), y_spread(0.0);
        for (Node *node: data) {
            x_mean += node->getLon();
            y_mean += node->getLat();
        }
        x_mean = x_mean / (double) data.size();
        y_mean = y_mean / (double) data.size();
        for (Node *node: data) {
            x_spread += pow(node->getLon() - x_mean, 2);
            y_spread += pow(node->getLat() - y_mean, 2);
        }
        int pivotIndex = (int) data.size() / 2;
        Node *pivotNode = data[pivotIndex];
        if (x_spread > y_spread) {
            sort(data.begin(), data.end(), [](const Node *a, Node *b) {
                return a->getLon() < b->getLon();
            });
        } else {
            sort(data.begin(), data.end(), [](const Node *a, Node *b) {
                return a->getLat() < b->getLat();
            });
        }
        root = new BallTreeNode(data[pivotIndex], data.size());
        root->left = createBallTree(vector<Node *>(data.begin(), data.begin() + pivotIndex));
        if (pivotIndex + 1 < data.size())
            root->right = createBallTree(vector<Node *>(data.begin() + pivotIndex + 1, data.end()));

        for (Node *node: data) {
            root->radius = max(root->radius, distance(node, pivotNode));
        }
    }
    return root;
}

void Graph::nearestNode(BallTreeNode *ballTreeNode, pair<double, double> &target, Node *&node, int &minDistance) {
    if (ballTreeNode == nullptr)
        return;
    int curDistance = distance(ballTreeNode->node, target);

    if (curDistance - ballTreeNode->radius >= minDistance)
        return;

    if (curDistance < minDistance) {
        node = ballTreeNode->node;
        minDistance = curDistance;
    }
    nearestNode(ballTreeNode->left, target, node, minDistance);
    nearestNode(ballTreeNode->right, target, node, minDistance);
}

Node *Graph::nearestNode(pair<double, double> target) {
    if (ballTree == nullptr) {
        throw runtime_error("BallTree not created.");
    }

    Node *res = ballTree->node;
    int minDistance = INT_MAX;
    nearestNode(ballTree, target, res, minDistance);
    return res;
}

Node* Graph::getRandomNode() {
    auto it = this->nodes.begin();
    int random = rand() % this->nodes.size();
    std::advance(it, random);
    return it->second;
}

Node* nearNode(BallTreeNode *ballTreeNode, pair<double, double> &target, int maxDistance) {

    if (ballTreeNode == nullptr)
        return nullptr;

    int curDistance = distance(ballTreeNode->node, target);
    if (curDistance - ballTreeNode->radius > maxDistance) {
        return nullptr;
    }

    if (curDistance < maxDistance) {
        return ballTreeNode->node;
    }

    Node* node;
    if (rand() % 2 == 1) {
        node = nearNode(ballTreeNode->left, target, maxDistance);
        if (node != nullptr){
            return node;
        }
        node = nearNode(ballTreeNode->right, target, maxDistance);
    } else {
        node = nearNode(ballTreeNode->left, target, maxDistance);
        if (node != nullptr){
            return node;
        }
        node = nearNode(ballTreeNode->right, target, maxDistance);
    }
    return node;
}
pair<Node*, Node*> Graph::getTwoNearNodes(int maxDistance){
    Node *node, *node2;

    // retry up to n times to search for two nodes in distance within good range
    // retry up to n2 times before giving same node
    int n = 20;
    int n2 = 100;
    for (int i=1; i<=n || (i <= n2 && (node2 == nullptr || node2 == node)); i++){
        node = this->getRandomNode();
        pair<double, double> target = {node->getLon(), node->getLat()};
        node2 = nearNode(this->ballTree, target, maxDistance);
        if (node2 != nullptr && distance(node, node2) > maxDistance * 2 / 3){
            break;
        }
    }
    return {node, node2};
}


bool Graph::hasNode(Node* node) {
    return this->nodes.find(node->getId()) != this->nodes.end();
}

int Graph::getEdgeCount(){
    return this->edgeCount;
}
