#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;

enum direction {
    bidirectional,
    oneway,
    onewayReversed,
    none
};

class Node {

public:
    class Edge {
        long id;
        int length; // saving in the unit of centimeters
        int elevation;
        Node *u;
        Node *v;

    public:
        Edge(Node *v);

        Edge(Node *node, bool reversed);

        Edge(long id, Node *u, Node *v,
             int length,
             int elevation);
        ~Edge();



        long getId() const;

        int getLength();

        int getElevation();

        Node* getU();

        Node* getV();

        static int getLength(Edge *edge);

        static int getElevation(Edge *edge);


    };

private:
    long id;
    double lon;
    double lat;
    int elevation;
    unordered_map<long, Edge *> edges;
    unordered_map<long, Edge *> reversedEdges;
    unordered_map<Edge *, unordered_set<Edge *>> *restrictions; // {previous edge, restricted next edge}
    unordered_map<Edge *, unordered_set<Edge *>> *reversedRestrictions;

public:

    Node(long id, double lon, double lat, int elevation);

    Node(long id, double lon, double lat, int elevation, unordered_map<long, Edge *> &edges);
    ~Node();

    Node *cloneWithoutEdges();

    bool operator==(const Node *node2) const;

    long getId();

    double getLon() const;

    double getLat() const;

    int getElevation() const;

    unordered_map<long, Edge *>& getEdges();

    unordered_map<long, Edge *>& getReversedEdges();

    bool hasEdge(long id);

    bool hasReversedEdge(long id);

    void setEdge(long id, Edge*);

    void setReversedEdge(long id, Edge*);

    unordered_map<Edge *, unordered_set<Edge *>> *getRestrictions();

    unordered_map<Edge *, unordered_set<Edge *>> *getReversedRestrictions();

    void setRestrictions(unordered_map<Edge *, unordered_set<Edge *>> * restrictions);
    void setReversedRestrictions(unordered_map<Edge *, unordered_set<Edge *>> * reversedRestrictions);

};

int distance(Node *node1, Node *node2);

int distance(Node *node1, pair<double, double> &p);

int distance(double lon1, double lon2, double lat1, double lat2);
