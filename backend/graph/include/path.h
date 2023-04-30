#include <vector>

#include "node.h"


const uint HASHMULT = 31; // prime number
const uint HASHMOD = 134217613; // prime number, HASHMOD * (HASHMULT + 1) < UINT_MAX


class Path {
private:

    vector<Node::Edge*> edges;
    int length;
    int elevation;
    bool confirmedlyOptimal;

public:
    Path();
    Path(vector<Node::Edge*>& edges, int length, int elevation);
    ~Path();
    void addEdge(Node::Edge *edge);
    void revert();
    Node* getStart();
    Node* getEnd();
    vector<Node::Edge*>& getEdges();
    int getLength();
    int getElevation();
    bool isConfirmedOptimal();
    void setLength(int length);
    void setElevation(int elevation);
    void confirmOptimal();
    bool empty() const;
};


class PathEdges {
public:

    class PathEdge {
    private:
        Node::Edge *edge;
        long length;
        int elevation;
        int referCount;
    public:
        PathEdge(Node::Edge *edge, int length, int elevation);
        ~PathEdge();
        Node::Edge* getEdge();
        int getLength();
        int getElevation();
        void increaseReferCount();
        bool decreaseReferCount(); // returns if refer count = 0
        PathEdge* clone();
    };


    struct biggerElevation {
        bool operator() (const PathEdges* pathEdges1, const PathEdges* pathEdges2){
            return pathEdges1->elevation > pathEdges2->elevation;
        }
    };

    struct smallerElevation {
        bool operator() (const PathEdges* pathEdges1, const PathEdges* pathEdges2){
            return pathEdges1->elevation < pathEdges2->elevation;
        }
    };

private:
    Node* start;
    Node* end;


    vector<PathEdge *> pathEdges;
    int length; // in centimeters
    int elevation;
    bool confirmedlyOptimal;

public:

    explicit PathEdges(Node* start);

    PathEdges(Node* start, vector<Node::Edge *>& edges);

    PathEdges(Node* start, vector<PathEdge *> &pathEdges, int index);

    explicit PathEdges(Path* path);

    ~PathEdges();

    bool empty() const;

    int size() const;

    int getLength() const;

    int getElevation() const;

    vector<PathEdge *>& getPathEdges();

    Node* getStart();

    Node* getEnd();

    PathEdge* at(size_t index);
    
    Node::Edge* getFirstEdge();

    Node::Edge* getLastEdge();

    bool isConfirmedOptimal();

    void addEdge(Node::Edge *edge);

    void confirmOptimal();

    PathEdges* cutBefore(int index);

    PathEdges* cutAfter(int index);

    PathEdges* randomCutReturnFront();

    Path* toPath();

};

struct pathEdgesHash{
public:
    size_t operator() (PathEdges *pathEdges) const;
};

struct pathEdgesEqual{
public:
    bool operator() (PathEdges *pathEdges1, PathEdges* pathEdges2) const;
};
