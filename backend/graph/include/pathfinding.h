#ifndef GRAPH_H
#define GRAPH_H
#include "graph.h"
#endif

enum DuplicateEdge {
    ignoreDuplicateEdges=0,
    minimizeDuplicateDirectedEdges=1,
    minimizeDuplicateUndirectedEdges=2,
    duplicateEdgeEnumCount = 3
};

// unordered_map<Node*, int> dijkstraAlgorithm
// Path* dijkstraAlgorithm
// unordered_set<Node*> getReachableNodes
// Path* elenaPathFindMinUsingDijkstra
// PathEdges* elenaPathSearchMaxUsingGeneticAlgorithm

int getMaxWeight(int curMinWeight, double maxWeightRatio);

unordered_map<Node*, int> dijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int &curMinWeight,
        double maxWeightRatio,
        int getWeight(Node::Edge *edge),
        bool reversed,
        unordered_map<Node*, Node::Edge*> *prevEdgeMap,
        unordered_set<Node*> *validNodes
);

unordered_map<Node::Edge *, int> edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int &curMinWeight,
        double maxWeightRatio,
        int getWeight(Node::Edge *edge),
        bool reversed,
        unordered_map<Node::Edge *, Node::Edge *> *prevEdgeMap,
        Node::Edge *&lastEdge,
        unordered_set<Node::Edge *> *validEdges
);

Path* pathFromPrevEdgeMap(
        Node* start,
        Node* end,
        unordered_map<Node*, Node::Edge*>& prevEdgeMap
);

Path* edgeBasedPathFromPrevEdgeMap(
        Node::Edge* lastEdge,
        unordered_map<Node::Edge *, Node::Edge *>& prevEdgeMap
);

Path *dijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int getWeight(Node::Edge *edge),
        unordered_set<Node *> *validNodes
);

Path *edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int getWeight(Node::Edge *edge),
        unordered_set<Node::Edge *> *validEdges
);

unordered_set<Node*> getReachableNodes(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
);

unordered_set<Node::Edge *> getReachableEdges(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
);

Path *findMinElevationUsingDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
);

Path *findMinElevationUsingEdgeBasedDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
);

Path* geneticAlgorithm(
        Graph *graph, Node *start, Node *end, int maxLength,
        unordered_map<Node*, int>& minLengthStart, unordered_map<Node*, int>& minLengthEnd,
        Path* shortestPath,
        size_t numProduce, size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch, bool maximize);

Path* edgeBasedGeneticAlgorithm(
        Graph *graph, Node *start, Node *end, int maxLength,
        unordered_map<Node::Edge*, int>& minLengthStart, unordered_map<Node::Edge*, int>& minLengthEnd,
        Path* shortestPath,
        size_t numProduce, size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch, bool maximize
        );



