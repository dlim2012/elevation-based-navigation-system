#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;

class Node {

public:
    class Edge {
    public:
        long id;
        int length; // saving in the unit of centimeters
        int elevation;
        Node *u;
        Node *v;

        Edge(Node *v);

        Edge(Node *node, bool reversed);

        Edge(long id, Node *u, Node *v,
             int length,
             int elevation);

        int getLength();

        int getElevation();

        static int getLength(Edge *edge);

        static int getElevation(Edge *edge);

    };

    long id;
    double lon;
    double lat;
    int elevation;
    unordered_map<long, Edge *> edges;
    unordered_map<long, Edge *> reversedEdges;
    unordered_map<Edge *, unordered_set<Edge *>> *restrictions; // {previous edge, restricted next edge}
    unordered_map<Edge *, unordered_set<Edge *>> *reversedRestrictions;


    Node(long id, double lon, double lat, int elevation);

    Node(long id, double lon, double lat, int elevation, unordered_map<long, Edge *> &edges);

    Node *cloneWithoutEdges();

    bool operator==(const Node *node2) const;
};

int distance(Node *node1, Node *node2);

int distance(Node *node1, pair<double, double> &p);

int distance(double lon1, double lon2, double lat1, double lat2);
