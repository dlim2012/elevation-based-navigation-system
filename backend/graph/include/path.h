
#include "node.h"


const uint HASHMULT = 31; // prime number
const uint HASHMOD = 134217613; // prime number, HASHMOD * (HASHMULT + 1) < UINT_MAX


class Path {
public:


    vector<Node::Edge*> edges;
    double length;
    int elevation;

    Path();
    Path(vector<Node::Edge*>& edges, double length, int elevation);
    void addEdge(Node::Edge *edge);
    void revert();
    Node* getStart();
    Node* getEnd();
};

class PathEdges {
public:

    struct PathEdge {
        Node::Edge *edge;
        double length;
        int elevation;

        PathEdge(Node::Edge *edge, double length, int elevation);
    };


    struct biggerElevation {
        bool operator() (const PathEdges* pathEdges1, const PathEdges* pathEdges2){
            return pathEdges1->elevation > pathEdges2->elevation;
        }
    };

    Node* start;
    vector<PathEdge *> pathEdges;
    double length;
    int elevation;

    PathEdges(Node* start);

    PathEdges(Node* start, vector<Node::Edge *>& edges);

    PathEdges(Node* start, vector<PathEdge *> &pathEdges, size_t index);

    int size() const;

    double getLength() const;

    int getElevation() const;

    Node* getStart();

    Node* getEnd();

    PathEdge* at(size_t index);

    Node::Edge* lastEdge();

    void addEdge(Node::Edge *edge);

    PathEdges* cutAfter(size_t index);

    PathEdges* randomCut();

    Path* toPath();

};

// todo: is this supposed to be declared in a header file?
struct pathEdgesHash{
public:
    size_t operator() (const PathEdges* pathEdges) const{
        size_t hash = 0;
        for (PathEdges::PathEdge* pathEdge: pathEdges->pathEdges){
            hash = (hash * HASHMULT + ((int) (pathEdge->edge->v->id % HASHMOD))) % HASHMOD;
        }
        return hash;
    }
};

struct pathEdgesEqual{
public:
    bool operator() (const PathEdges* pathEdges1, const PathEdges* pathEdges2) const {
        if (pathEdges1->size() != pathEdges2->size()){
            return false;
        }
        size_t size = pathEdges1->size();
        for (size_t i=0; i<size; i++){
            if (pathEdges1->pathEdges[i]->edge != pathEdges2->pathEdges[i]->edge){
                return false;
            }
        }
        return true;
    }
};