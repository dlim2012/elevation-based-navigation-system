#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <float.h>
#include <random>
#include <algorithm>
#include <pqxx/pqxx>

#include "graph.h"



using namespace std;
using namespace pqxx;

Graph::Graph(bool populateEdgesField) {
    this->populateEdgesField = populateEdgesField;
    srand(time(nullptr));
}

Graph::Graph(string &connectionString, MapConfig mapConfig) {
    this->populateEdgesField = true;
    connection C(connectionString);
    prepare_nodes_query(C, mapConfig);
    prepare_edges_query(C, mapConfig);
    prepare_restrictions_query(C, mapConfig);
    addAllNodesFromDB(C, false);
    addAllEdgesFromDB(C, false);
    addAllRestrictionsFromDB(C);
    maxGroup();
    this->ballTree = createBallTree();
}


void saveAllEdgesInCsv(string connectionString, MapConfig mapConfig) {
    unordered_map<long, vector<pair<double, double>>> umap = allEdges(connectionString);

    Graph *graph = new Graph(connectionString, mapConfig);

    ofstream edgesFile;
    edgesFile.open("csv_files/mass-edges.csv");

    edgesFile << fixed << setprecision(8) << "index,longitude1,latitude1,longitude2,latitude2" << endl;

    int index(0);
    for (auto pair1: graph->edges) {
        vector<pair<double, double>> &v = umap[pair1.first];
        for (int i = 0; i < v.size() - 1; i++) {
            edgesFile << index++ << "," << v[i].first << "," << v[i].second << "," << v[i + 1].first << ","
                      << v[i + 1].second << endl;
        }
    }
    edgesFile.close();
}


Node *Graph::addNode(long id, double lon, double lat, int elevation) {
    return this->nodes[id] = new Node(id, lon, lat, elevation);
}

void Graph::addNode(Node *node) {
    this->nodes[node->id] = node;
}



void Graph::addEdge(long id, Node *u, Node *v, double length, int elevation, direction direction) {
    if (direction == none) {
        return;
    }

    Node::Edge *edge1, *edge2;


    if (direction == bidirectional || direction == oneway) {
        edge1 = new Node::Edge(id, u, v, length, elevation);
        u->edges[id] = edge1;
        v->reversedEdges[id] = edge1;
    } else {
        edge1 = nullptr;
    }

    if (direction == bidirectional || direction == onewayReversed) {
        edge2 = new Node::Edge(id, v, u, length, elevation);
        v->edges[id] = edge2;
        u->reversedEdges[id] = edge2;
    } else {
        edge2 = nullptr;
    }

    if (this->populateEdgesField) {
        if (u->id < v->id) {
            this->edges[id] = {edge1, edge2};
        } else {
            this->edges[id] = {edge2, edge1};
        }
    }

    this->maxElevation = max(this->maxElevation, elevation);
}

int Graph::addAllNodesFromDB(pqxx::connection_base &C, bool stream) {
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
        cout << "Fetched all nodes data from DB (count: " << r.size() << ")" << endl;
        return r.size();
    }
}

int Graph::addAllEdgesFromDB(pqxx::connection_base &C, bool stream) {
    pqxx::work w(C);
    if (stream) {
        throw runtime_error("Not Implemented.");
    } else {
        pqxx::result r = w.exec_prepared("edges");
        w.commit();
        for (pqxx::row row1: r) {
            long id = stol(row1[0].c_str());
            double length = stod(row1[3].c_str());
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
            if (stoi(row1[5].c_str()) == 0) {
                direction = bidirectional;
            } else {
                direction = oneway;
            }
            this->addEdge(id, itU->second, itV->second, length, elevation, direction);
        }
        cout << "Fetched all edges data from DB (count: " << r.size() << ")" << endl;

        int count = 0;
        for (auto p: this->edges) {
            if (p.second.first != nullptr)
                count++;
            if (p.second.second != nullptr) {
                count++;
            }
        }
        cout << count << endl;
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
    Node *viaNode = fromEdge->v;
    // add the restrictions map to 'via' if not exists
    if (viaNode->restrictions == nullptr) {
        viaNode->restrictions = new unordered_map<Node::Edge *, unordered_set<Node::Edge *>>();
        viaNode->reversedRestrictions = new unordered_map<Node::Edge *, unordered_set<Node::Edge *>>();
    }
    unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions = viaNode->restrictions;
    unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *reversedRestrictions = viaNode->reversedRestrictions;

    if (type_enum[1] == 'n') {
        addRestriction(restrictions, fromEdge, toEdge);
        addRestriction(reversedRestrictions, toEdge, fromEdge);
    } else {
        for (pair<const long, Node::Edge *> &pair1: viaNode->edges) {
            if (pair1.second != toEdge) {
                addRestriction(restrictions, fromEdge, pair1.second);
                addRestriction(reversedRestrictions, pair1.second, fromEdge);
            }
        }
    }
}

int Graph::addAllRestrictionsFromDB(pqxx::connection_base &C) {
    // Note: the populateEdgesField has to be true to add restrictions to edges
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
            if (edgePair.first != nullptr && edgePair.first->v == via) {
                from = edgePair.first;
            } else if (edgePair.second != nullptr && edgePair.second->v == via) {
                from = edgePair.second;
            }
            if (from == nullptr)
                continue;

            // find the 'to' edge (skip if not found)
            auto nodeEdgesIt = via->edges.find(to_id);
            if (nodeEdgesIt == via->edges.end())
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
                auto nodeEdgesIt = via->v->edges.find(to_id);
                if (nodeEdgesIt == via->v->edges.end())
                    continue;
                to = nodeEdgesIt->second;
            } else if (p.second == nullptr) {
                via = p.first;
                auto nodeEdgesIt = via->v->edges.find(to_id);
                if (nodeEdgesIt == via->v->edges.end())
                    continue;
            } else {
                auto nodeEdgesIt = p.first->v->edges.find(to_id);
                if (nodeEdgesIt != p.first->v->edges.end()) {
                    to = nodeEdgesIt->second;
                    via = p.first;
                } else {
                    nodeEdgesIt = p.second->v->edges.find(to_id);
                    if (nodeEdgesIt != p.second->v->edges.end()) {
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
            if (p.first != nullptr && p.first->v == via->u) {
                from = p.first;
            } else if (p.second != nullptr && p.second->v == via->u) {
                from = p.second;
            } else
                continue;

            addRestrictions(from, via, type_enum);
            addRestrictions(via, to, type_enum);

        }

    }
    cout << "Fetched all restrictions data from DB (count: " << r.size() << ")" << endl;
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
    stack.push({edge, nullptr, edge->v->edges.begin()});

    while (!stack.empty()) {
        Node::Edge *parEdge = stack.top().parEdge;
        Node::Edge *childEdge = stack.top().childEdge;
        unordered_map<long, Node::Edge *>::iterator itCur = stack.top().itCur;
        stack.pop();

        pair<int, int> &parMemo = memo[parEdge];
        if (childEdge != nullptr) {
            parMemo.second = min(parMemo.second, memo[childEdge].second);
        }

        while (itCur != parEdge->v->edges.end()) {
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
        if (itCur != parEdge->v->edges.end()) {
            tarjanStack.push(childEdge);
            onTarjanStack.insert(childEdge);
            memo[childEdge] = {index, index};
            index++;
            stack.push({parEdge, childEdge, ++itCur});
            stack.push({childEdge, nullptr, childEdge->v->edges.begin()});
        } else if (parMemo.first == parMemo.second) {
            unordered_map<long, pair<Node::Edge *, Node::Edge *>> edgeGroup;
            unordered_map<long, Node *> nodeGroup;

            while (true) {
                edge = tarjanStack.top();
                tarjanStack.pop();
                onTarjanStack.erase(edge);

                auto edgePairIt = edgeGroup.find(edge->id);
                if (edgePairIt == edgeGroup.end()) {

                    if (edge->u->id < edge->v->id) {
                        edgeGroup[edge->id] = {edge, nullptr};
                    } else {
                        edgeGroup[edge->id] = {nullptr, edge};
                    }
                } else {
                    if (edge->u->id < edge->v->id) {
                        edgePairIt->second.first = edge;
                    } else {
                        edgePairIt->second.second = edge;
                    }
                }

                if (nodeGroup.find(edge->u->id) == nodeGroup.end()) {
                    nodeGroup[edge->u->id] = edge->u;
                }
                if (nodeGroup.find(edge->v->id) == nodeGroup.end()) {
                    nodeGroup[edge->v->id] = edge->v;
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


    int edgeCount1 = 0;
    int edgeCount2 = 0;
    for (pair<const long, Node *> &pair1: newNodes) {
        Node *node = pair1.second;
        // remove unused edges
        for (auto edgeIt = node->edges.begin(); edgeIt != node->edges.end();) {
            Node::Edge *edge = edgeIt->second;
            auto nodeEdgesIt = newEdges.find(edge->id);
            if (nodeEdgesIt == newEdges.end()) {
                edgeIt = node->edges.erase(edgeIt);
            } else if (edge != nodeEdgesIt->second.first && edge != nodeEdgesIt->second.second) {
                edgeIt = node->edges.erase(edgeIt);
            } else {
                edgeIt++;
                edgeCount1++;
            }
        }
        // remove unused reversed edges
        for (auto edgeIt = node->reversedEdges.begin(); edgeIt != node->reversedEdges.end();) {
            Node::Edge *edge = edgeIt->second;
            auto nodeEdgesIt = newEdges.find(edge->id);
            if (nodeEdgesIt == newEdges.end()) {
                edgeIt = node->reversedEdges.erase(edgeIt);
            } else if (edge != nodeEdgesIt->second.first && edge != nodeEdgesIt->second.second) {
                edgeIt = node->reversedEdges.erase(edgeIt);
            } else {
                edgeIt++;
                edgeCount1++;
            }
        }

    }


    int countNull(0);
    for (auto it: newEdges) {
        if (it.second.first == nullptr && it.second.second == nullptr)
            countNull++;
    }
    cout << "null " << countNull << endl;

    swap(this->nodes, newNodes);
    swap(this->edges, newEdges);

    cout << "Extracted the biggest strongly connected components (node: " << this->nodes.size() << ", edge: "
         << edgeCount1 << ")" << endl;

}

bool Graph::isRestricted(Node::Edge *p, Node::Edge *c, bool reversed) {
    if (p == nullptr){
        return false;
    }
    unordered_map<Node::Edge *, unordered_set<Node::Edge *>> *restrictions =
            reversed ? p->u->reversedRestrictions : p->v->restrictions;

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

Node* Graph::getRandomNode() {
    auto it = this->nodes.cbegin();
    int random = rand() % this->nodes.size();
    std::advance(it, random);
    return it->second;
}

BallTreeNode *Graph::createBallTree() {
    vector<Node *> data;
    for (auto it: nodes) {
        data.push_back(it.second);
    }
    if (data.empty()) {
        cout << "Graph is empty." << endl;
        return nullptr;
    }
    this->ballTree = createBallTree(data);
    cout << "BallTree created." << endl;
    return this->ballTree;
}

BallTreeNode *Graph::createBallTree(vector<Node *> data) {
    BallTreeNode *root;
    if (data.size() == 1) {
        root = new BallTreeNode(data[0]);
    } else {
        double x_mean(0.0), y_mean(0.0), x_spread(0.0), y_spread(0.0);
        for (Node *node: data) {
            x_mean += node->lon;
            y_mean += node->lat;
        }
        x_mean = x_mean / (double) data.size();
        y_mean = y_mean / (double) data.size();
        for (Node *node: data) {
            x_spread += pow(node->lon - x_mean, 2);
            y_spread += pow(node->lat - y_mean, 2);
        }
        int pivotIndex = (int) data.size() / 2;
        Node *pivotNode = data[pivotIndex];
        if (x_spread > y_spread) {
            sort(data.begin(), data.end(), [](const Node *a, Node *b) {
                return a->lon < b->lon;
            });
        } else {
            sort(data.begin(), data.end(), [](const Node *a, Node *b) {
                return a->lat < b->lat;
            });
        }
        root = new BallTreeNode(data[pivotIndex]);
        root->left = createBallTree(vector<Node *>(data.begin(), data.begin() + pivotIndex));
        if (pivotIndex + 1 < data.size())
            root->right = createBallTree(vector<Node *>(data.begin() + pivotIndex + 1, data.end()));

        for (Node *node: data) {
            root->radius = max(root->radius, distance(node, pivotNode));
        }
    }
    return root;
}

void Graph::nearestNode(BallTreeNode *ballTreeNode, pair<double, double> &target, Node *&node, double &minDistance) {
    if (ballTreeNode == nullptr)
        return;
    if (distance(ballTreeNode->node, target) - ballTreeNode->radius >= minDistance)
        return;

    double curDistance = distance(ballTreeNode->node, target);
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
    double minDistance = DBL_MAX;
    nearestNode(ballTree, target, res, minDistance);
    return res;
}
