#ifndef GRAPH_H
#define GRAPH_H
#include "graph.h"
#endif

unordered_map<Node::Edge *, double> edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        double &curMinWeight,
        double maxWeightRatio,
        double getWeight(Node::Edge *edge),
        bool reversed,
        unordered_map<Node::Edge *, Node::Edge *> *prevEdgeMap,
        Node::Edge *&lastEdge,
        unordered_set<Node::Edge *> *validEdges
);

Path *edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        double getWeight(Node::Edge *edge),
        bool reversed,
        unordered_set<Node::Edge *> *validEdges
);


unordered_set<Node::Edge *> getReachableEdges(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
);


Path *elenaPathFindMinUsingDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
);

PathEdges* elenaPathSearchMaxUsingGeneticAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        double maxLengthRatio,
        size_t numProduce,
        size_t numMaxSelect,
        int numEpoch
        );



