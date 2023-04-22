#ifndef GRAPH_H
#define GRAPH_H
#include "graph.h"
#endif

enum DuplicateEdge {
    ignoreDuplicateEdges=0,
    avoidDuplicateDirectedEdges=1,
    unallowDuplicateDirectedEdges=2,
    avoidDuplicateUndirectedEdges=3,
    unallowDuplicateUndirectedEdges=4,
    duplicateEdgeEnumCount = 5
};

// unordered_map<Node*, int> dijkstraAlgorithm
// Path* dijkstraAlgorithm
// unordered_set<Node*> getReachableNodes
// Path* elenaPathFindMinUsingDijkstra
// PathEdges* elenaPathSearchMaxUsingGeneticAlgorithm

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
        double maxWeightRatio,
        int maxLength
);

unordered_set<Node::Edge *> getReachableEdges(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio,
        int maxLength
);

Path *elenaPathFindMinUsingDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio,
        int maxLength
);

Path *elenaPathFindMinUsingEdgeBasedDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio,
        int maxLength
);

PathEdges* elenaPathSearchMaxUsingGeneticAlgorithm(
        Graph *graph, Node *start, Node *end, double maxLengthRatio, size_t numProduce,
        size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge, int maxLength,
        int maxMilliseconds);

PathEdges* elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        double maxLengthRatio,
        size_t numProduce,
        size_t numMaxSelect,
        int numEpoch,
        DuplicateEdge duplicateEdge,
        int maxLength,
        int maxMilliseconds
        );



