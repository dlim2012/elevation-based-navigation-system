#include <float.h>
#include <unordered_map>
#include <random>
#include <deque>
#include <chrono>

#include "pathfinding.h"


#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

using namespace std;

int getMaxWeight(int curMinWeight, double maxWeightRatio){
    if (maxWeightRatio <= 0){
        return curMinWeight;
    }
    return curMinWeight > (INT_MAX - 1) / maxWeightRatio ? INT_MAX : (int) curMinWeight * maxWeightRatio;
}

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
        ){

    if (graph->nodes.find(start->id) == graph->nodes.end()
        || graph->nodes.find(end->id) == graph->nodes.end()
            ) {
        throw invalid_argument("nodes do not belong to the same graph");
    }


    struct Candidate {
        int weight;
        Node *prevNode;

        Candidate(int weight, Node *prevNode) {
            this->weight = weight;
            this->prevNode = prevNode;
        }
    };

    struct comparator {
        bool operator()(const Candidate a, const Candidate b) {
            return a.weight > b.weight;
        }
    };

    unordered_map<Node *, int> minWeightNode;
    unordered_set<Node *> visited;

    priority_queue<Candidate, vector<Candidate>, comparator> pq;

    Node* curNode = reversed ? end : start;
    minWeightNode[curNode] = 0;
    pq.emplace(0, curNode);
    int curMaxWeight = getMaxWeight(curMinWeight, maxWeightRatio);

    while (!pq.empty()) { // measured in a modern computer: 3~4 microseconds per iteration
        Candidate candidate = pq.top();
        pq.pop();

        if (visited.find(candidate.prevNode) != visited.end())
            continue;
        visited.insert(candidate.prevNode);

        for (pair<const long, Node::Edge *> &pair1: reversed ? candidate.prevNode->reversedEdges
                                                             : candidate.prevNode->edges) {
            curNode = reversed ? pair1.second->u : pair1.second->v;
            if (validNodes != nullptr && validNodes->find(curNode) == validNodes->end())
                continue;

            if (visited.find(curNode) != visited.end()){
                continue;
            }

            int curWeight = candidate.weight + getWeight(pair1.second);
            if (curWeight > curMaxWeight) {
                continue;
            }

            auto it = minWeightNode.find(curNode);
            if (it == minWeightNode.end() || it->second > curWeight) {

                pq.emplace(curWeight, curNode);

                // update min weight
                minWeightNode[curNode] = curWeight;

                // update previous edge
                if (prevEdgeMap != nullptr) {
                    (*prevEdgeMap)[curNode] = pair1.second;
                }

                // update when arrived at the target
                if ((reversed && curNode->id == start->id)
                    || (!reversed && curNode->id == end->id)) {
                    if (curMinWeight >= curWeight) {
                        curMinWeight = curWeight;
                        curMaxWeight = getMaxWeight(curMinWeight, maxWeightRatio);
                    }
                }
            }
        }
    }
    return minWeightNode;
}

unordered_map<Node::Edge *, int> edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int &curMinWeight,
        double maxWeightRatio,
        int getWeight(Node::Edge *edge),
        bool reversed, // if true: use incomingEdges instead
        unordered_map<Node::Edge *, Node::Edge *> *prevEdgeMap,
        Node::Edge *&lastEdge,
        unordered_set<Node::Edge *> *validEdges
) {

    if (graph->nodes.find(start->id) == graph->nodes.end()
        || graph->nodes.find(end->id) == graph->nodes.end()
            ) {
        throw invalid_argument("nodes do not belong to the same graph");
    }

    struct Candidate {
        int weight;
        Node::Edge *prevEdge;

        Candidate(int weight, Node::Edge *prevEdge) {
            this->weight = weight;
            this->prevEdge = prevEdge;
        }
    };

    struct comparator {
        bool operator()(const Candidate a, const Candidate b) {
            return a.weight > b.weight;
        }
    };

    unordered_map<Node::Edge *, int> minWeightEdge;
    unordered_set<Node::Edge*> visited;

    priority_queue<Candidate, vector<Candidate>, comparator> pq;

    Node::Edge *curEdge = new Node::Edge(reversed ? end : start, reversed);
    minWeightEdge[curEdge] = 0;
    pq.emplace(0, curEdge);
    int curMaxWeight = getMaxWeight(curMinWeight, maxWeightRatio);

    while (!pq.empty()) {
        Candidate candidate = pq.top();
        pq.pop();

        if (visited.find(candidate.prevEdge) != visited.end())
            continue;
        visited.insert(candidate.prevEdge);

        for (pair<const long, Node::Edge *> &pair1: reversed ? candidate.prevEdge->u->reversedEdges
                                                             : candidate.prevEdge->v->edges) {
            curEdge = pair1.second;
            if (validEdges != nullptr && validEdges->find(curEdge) == validEdges->end())
                continue;

            if (visited.find(curEdge) != visited.end()){
                continue;
            }

            if (Graph::isRestricted(candidate.prevEdge, curEdge, reversed)) {
                continue;
            }

            int curWeight = candidate.weight + getWeight(curEdge);
            if (curWeight > curMaxWeight) {
                continue;
            }

            auto it = minWeightEdge.find(curEdge);
            if (it == minWeightEdge.end() || it->second > curWeight) {

                pq.emplace(curWeight, curEdge);

                // update min weight
                minWeightEdge[curEdge] = curWeight;

                // update previous edge
                if (prevEdgeMap != nullptr) {
                    (*prevEdgeMap)[curEdge] = candidate.prevEdge;
                }

                // update when arrived at the target
                if ((reversed && curEdge->u->id == start->id)
                    || (!reversed && curEdge->v->id == end->id)) {
                        if (curMinWeight >= curWeight) {
                        curMinWeight = curWeight;
                        curMaxWeight = getMaxWeight(curMinWeight, maxWeightRatio);
                        lastEdge = curEdge;
                    }
                }
            }
        }
    }
    return minWeightEdge;
}

Path *dijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int getWeight(Node::Edge *edge),
        unordered_set<Node *> *validNodes
        ){
    Path* path = new Path();
    if (start == end) {
        return path;
    }
    unordered_map<Node *, Node::Edge *> prevEdgeMap;
    int curMinWeight = INT_MAX;


    auto t0 = chrono::high_resolution_clock::now();
    // find the path reversed so that the path can be reconstructed non-reversed
    unordered_map<Node *, int> minWeightEdge = dijkstraAlgorithm(
            graph, start, end, curMinWeight, -1.0, getWeight, true, &prevEdgeMap, validNodes);
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;

    Node::Edge* prevEdge;
    while (start != end) {
        prevEdge = prevEdgeMap[start];
        path->addEdge(prevEdge);
        start = prevEdge->v;
    }

    return path;
}

Path *edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int getWeight(Node::Edge *edge),
        unordered_set<Node::Edge *> *validEdges
) {
    Path* path = new Path();
    if (start == end) {
        return path;
    }
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMap;
    Node::Edge *edge = nullptr;
    int curMinWeight = INT_MAX;
    // find the path reversed so that the path can be reconstructed non-reversed
    unordered_map<Node::Edge *, int> minWeightEdge = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, -1.0, getWeight, true, &prevEdgeMap, edge, validEdges);


    if (edge == nullptr){
        throw runtime_error("path not found.");
    }

    while (edge->id != -1) {
        path->addEdge(edge);
        edge = prevEdgeMap[edge];
    }

    return path;
}

unordered_set<Node*> getReachableNodes(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
        ){


    if (graph->nodes.find(start->id) == graph->nodes.end()
        || graph->nodes.find(end->id) == graph->nodes.end()
            ) {
        throw invalid_argument("start and/or end nodes do not belong to the graph");
    }
    unordered_set<Node *> validNodes;

    int curMinWeight = INT_MAX;

    auto t0 = chrono::high_resolution_clock::now();

    unordered_map<Node *, int> minWeightStart = dijkstraAlgorithm(
            graph, start, end, curMinWeight, maxWeightRatio, Node::Edge::getLength, false, nullptr, nullptr
    );

    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;

    // If all edges are reachable, return empty set to indicate this fact
    if (minWeightStart.size() == graph->edgeCount + 1){
        return validNodes;
    }

    unordered_map<Node *, int> minWeightEnd = dijkstraAlgorithm(
            graph, start, end, curMinWeight, maxWeightRatio, Node::Edge::getLength, true, nullptr, nullptr
    );
    t1 = chrono::high_resolution_clock::now();
    t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;

    int maxWeight = getMaxWeight(curMinWeight, maxWeightRatio);
    for (pair<Node *const, int> &pair1: minWeightEnd) {
        auto it = minWeightStart.find(pair1.first);
        if (it == minWeightStart.end())
            continue;

        if (pair1.second + it->second > maxWeight)
            continue;

        validNodes.insert(pair1.first);
    }
    return validNodes;
}

unordered_set<Node::Edge *> getReachableEdges(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    if (graph->nodes.find(start->id) == graph->nodes.end()
        || graph->nodes.find(end->id) == graph->nodes.end()
            ) {
        throw invalid_argument("start and/or end nodes do not belong to the graph");
    }
    unordered_set<Node::Edge *> validEdges;

    int curMinWeight = INT_MAX;
    Node::Edge *lastEdge = nullptr;

    unordered_map<Node::Edge *, int> minWeightStart = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, maxWeightRatio, Node::Edge::getLength, false, nullptr, lastEdge, nullptr
    );
    if (lastEdge == nullptr) {
        throw runtime_error("runtime error");
    }

    // If all edges are reachable, return empty set to indicate this fact
    if (minWeightStart.size() == graph->edgeCount + 1){
        return validEdges;
    }

    lastEdge = nullptr;
    unordered_map<Node::Edge *, int> minWeightEnd = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, maxWeightRatio, Node::Edge::getLength, true, nullptr, lastEdge, nullptr
    );

    if (lastEdge == nullptr) {
        throw runtime_error("runtime error");
    }

    int maxWeight = getMaxWeight(curMinWeight, maxWeightRatio);
    for (pair<Node::Edge *const, int> &pair1: minWeightEnd) {
        auto it = minWeightStart.find(pair1.first);
        if (it == minWeightStart.end())
            continue;

        if (pair1.second + it->second - Node::Edge::getLength(pair1.first) > maxWeight)
            continue;

        validEdges.insert(pair1.first);
    }
    return validEdges;
}


Path *elenaPathFindMinUsingDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    unordered_set<Node *> validNodes = getReachableNodes(graph, start, end, maxWeightRatio);
    if (validNodes.empty()){
        return new Path();
    }
    return dijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, validNodes.empty() ? nullptr : &validNodes);
}

Path *elenaPathFindMinUsingEdgeBasedDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    unordered_set<Node::Edge *> validEdges = getReachableEdges(graph, start, end, maxWeightRatio);
    if (validEdges.empty()){
        return new Path();
    }
    return edgeBasedDijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, validEdges.empty() ? nullptr : &validEdges);
}

PathEdges* removeDuplicateDirectedEdge(PathEdges* pathEdges){
    if (pathEdges->pathEdges.empty()){
        return pathEdges;
    }

    // find the last index of each edge
    unordered_map<Node::Edge*, int> umap;
    int i = 0;
    vector<PathEdges::PathEdge*> &pathEdgesVector = pathEdges->pathEdges;
    for (PathEdges::PathEdge* pathEdge: pathEdgesVector){
        umap[pathEdge->edge] = i++;
    }

    vector<pair<int,int>> intervals;
    int l(0);
    for (int r=0; r<pathEdgesVector.size(); r++){
        i = umap[pathEdgesVector[r]->edge];
        if (r != i){
            intervals.emplace_back(l+1, r+1);
            l = r = i;
        }
    }

    if (intervals.empty()){
        return pathEdges;
    }
    intervals.emplace_back(l+1, pathEdgesVector.size());
    PathEdges *newPathEdges = new PathEdges(pathEdges->start);
    newPathEdges->addEdge(pathEdgesVector[0]->edge);
    for (pair<int, int>& p: intervals){
        for (int i=p.first; i<p.second; i++){
            newPathEdges->addEdge(pathEdgesVector[i]->edge);
        }
    }

    return newPathEdges;
}


PathEdges* removeDuplicateUndirectedEdge(PathEdges* pathEdges){
    // Note: no restrictions applied

    if (pathEdges->pathEdges.empty()){
        return pathEdges;
    }

    // find the last index of each edge
    unordered_map<long, int> umap;
    int i = 0;
    vector<PathEdges::PathEdge*> &pathEdgesVector = pathEdges->pathEdges;
    for (PathEdges::PathEdge* pathEdge: pathEdgesVector){
        umap[pathEdge->edge->id] = i++;
    }

    vector<pair<int,int>> intervals;
    int l(0);
    for (int r=0; r<pathEdgesVector.size(); r++){
        i = umap[pathEdgesVector[r]->edge->id];
        if (r != i){
            if (pathEdgesVector[i]->edge == pathEdgesVector[r]->edge){
                intervals.emplace_back(l+1, r+1);
                l = r = i;
            } else {
                if (l+1 < r) {
                    intervals.emplace_back(l + 1, r);
                    Node::Edge* e1 = pathEdgesVector[r]->edge;
                    Node::Edge* e2 = pathEdgesVector[i]->edge;
                    l = r = i;
                }
            }
        }
    }

    if (intervals.empty()){
        return pathEdges;
    }

    intervals.emplace_back(l+1, pathEdgesVector.size());

    PathEdges *newPathEdges = new PathEdges(pathEdges->start);
    newPathEdges->addEdge(pathEdgesVector[0]->edge);
    for (pair<int, int>& p: intervals){
        for (int i=p.first; i<p.second; i++){
            newPathEdges->addEdge(pathEdgesVector[i]->edge);
        }
    }

    return newPathEdges;
}

PathEdges* mutatePathEdges(
        PathEdges* pathEdges,
        unordered_map<Node *, int>& minLengthStart,
        unordered_map<Node *, int>& minLengthEnd,
        int maxLength,
        DuplicateEdge duplicateEdge){
    unordered_set<long> visitedEdgeId;
    unordered_set<Node::Edge*> visitedEdge;
    bool useVisitedEdgeId = duplicateEdge == avoidDuplicateUndirectedEdges || duplicateEdge == unallowDuplicateUndirectedEdges;
    bool useVisitedEdge = duplicateEdge == avoidDuplicateDirectedEdges || duplicateEdge == unallowDuplicateDirectedEdges;
    bool avoid = useVisitedEdgeId || useVisitedEdge;

    PathEdges* newPathEdges;
    if (rand() & 1) {
        int targetCutIndex = rand() % (pathEdges->size() + 1);
        int index = 0;
        if (avoid){
            for (PathEdges::PathEdge* pathEdge: pathEdges->pathEdges) {
                if (index == targetCutIndex
                    || (duplicateEdge == unallowDuplicateUndirectedEdges
                        && visitedEdgeId.find(pathEdge->edge->id) != visitedEdgeId.end())
                    || (duplicateEdge == unallowDuplicateDirectedEdges
                        && visitedEdge.find(pathEdge->edge) != visitedEdge.end())) {
                    break;
                }
                index++;
                if (useVisitedEdgeId){
                    visitedEdgeId.insert(pathEdge->edge->id);
                } else if (useVisitedEdge){
                    visitedEdge.insert(pathEdge->edge);
                }
            }
        }
        newPathEdges = pathEdges->cutBefore(index);

        Node::Edge *prevEdge = newPathEdges->lastEdge();
        while (true) {
            vector<Node::Edge*> nextEdgeCandidates;
            vector<Node::Edge*> nextEdgeCandidatesVisited;
            for (pair<const long, Node::Edge *> pair1: newPathEdges->getEnd()->edges) {
                Node::Edge *curEdge = pair1.second;

                auto it = minLengthEnd.find(curEdge->v);

                if (it == minLengthEnd.end() || newPathEdges->length + curEdge->length + it->second > maxLength) {
                    continue;
                }

                if (duplicateEdge == ignoreDuplicateEdges
                    || (useVisitedEdgeId && visitedEdgeId.find(curEdge->id) == visitedEdgeId.end())
                    || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())){
                    nextEdgeCandidates.push_back(curEdge);
                } else if (avoid){
                    nextEdgeCandidatesVisited.push_back(curEdge);
                }
            }

            if (newPathEdges->getEnd() == pathEdges->getEnd()){
                break;
            }

            // even if current location is the target location, make moves if possible to maximize elevation gain
            if (nextEdgeCandidates.empty() && nextEdgeCandidatesVisited.empty()) {
                if (newPathEdges->getEnd() != pathEdges->getEnd()){
                    if (duplicateEdge != unallowDuplicateUndirectedEdges
                            && duplicateEdge != unallowDuplicateDirectedEdges){
                        throw runtime_error("error"); // temporal error statement
                    }
                    return pathEdges;
                }
                break;
            }

            // visit unvisited node first if not ignoring duplicate edges
            if (!nextEdgeCandidates.empty()) {
                prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
            } else {
                prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
            }


            newPathEdges->addEdge(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->id);
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }

        }
    } else {
        int targetCutRIndex = rand() % (pathEdges->size() + 1);
        int rIndex = 0;
        vector<Node::Edge*> newEdges;

        if (avoid) {
            for (auto it=pathEdges->pathEdges.rbegin(); rIndex < targetCutRIndex ; it++){
                if (duplicateEdge == unallowDuplicateUndirectedEdges
                    && visitedEdgeId.find((*it)->edge->id) == visitedEdgeId.end()
                    || (duplicateEdge == unallowDuplicateDirectedEdges
                        && visitedEdge.find((*it)->edge) != visitedEdge.end())){
                    break;
                }
                rIndex++;
                if (useVisitedEdgeId){
                    visitedEdgeId.insert((*it)->edge->id);
                } else if (useVisitedEdge){
                    visitedEdge.insert((*it)->edge);
                }
            }
        }
        PathEdges::PathEdge *prevPathEdge = (rIndex == pathEdges->size()) ?
                                            nullptr :  pathEdges->at(pathEdges->size()-1-rIndex);
        Node::Edge *prevEdge = prevPathEdge == nullptr ? nullptr : prevPathEdge->edge;
        int length = prevPathEdge == nullptr ? 0 : pathEdges->length - prevPathEdge->length + prevEdge->length;
        Node *prevNode = prevPathEdge == nullptr ? pathEdges->end : prevEdge->u;
        while (true) {
            vector<Node::Edge*> nextEdgeCandidates;
            vector<Node::Edge*> nextEdgeCandidatesVisited;
            for (pair<const long, Node::Edge *> pair1: prevNode->reversedEdges) {
                Node::Edge *curEdge = pair1.second;

                if (Graph::isRestricted(prevEdge, curEdge, true)) {
                    continue;
                }

                auto it = minLengthStart.find(curEdge->u);
                if (it == minLengthStart.end() || length + + curEdge->length + it->second > maxLength) {
                    continue;
                }


                if (duplicateEdge == ignoreDuplicateEdges
                    || (useVisitedEdgeId && visitedEdgeId.find(curEdge->id) == visitedEdgeId.end())
                    || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())   ){
                    nextEdgeCandidates.push_back(curEdge);
                } else if (avoid){
                    nextEdgeCandidatesVisited.push_back(curEdge);
                }

            }

            if (!newEdges.empty() && newEdges.back()->u == pathEdges->getStart()){
                break;
            }

            // even if current location is the target location, make moves if possible to maximize elevation gain
            if (nextEdgeCandidates.empty() && nextEdgeCandidatesVisited.empty()) {
                if (!newEdges.empty() && newEdges.back()->u != pathEdges->getStart()) {
                    return pathEdges;
                }
                break;
            }

            // visit unvisited node first if minimizeRevisit == true
            if (!nextEdgeCandidates.empty()) {
                prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
            } else {
                prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
            }

            newEdges.push_back(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->id);
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
            length += prevEdge->length;
            prevNode = prevEdge->u;
        }


        newPathEdges = new PathEdges(pathEdges->start);
        for (auto it=newEdges.rbegin(); it!=newEdges.rend(); it++){
            newPathEdges->addEdge(*it);
        }
        if (rIndex < pathEdges->size()) {
            auto it = pathEdges->pathEdges.rbegin() + rIndex;
            while (true) {
                newPathEdges->addEdge((*it)->edge);
                if (it-- == pathEdges->pathEdges.rbegin()) {
                    break;
                }
            }
        }
    }

    if (duplicateEdge == unallowDuplicateUndirectedEdges) {
        newPathEdges = removeDuplicateUndirectedEdge(newPathEdges);
    } else if (duplicateEdge == unallowDuplicateDirectedEdges){
        newPathEdges = removeDuplicateDirectedEdge(newPathEdges);
    }
    return newPathEdges;
}

PathEdges* edgeBasedMutatePathEdges(
        PathEdges* pathEdges,
        unordered_map<Node::Edge *, int>& minLengthStart,
        unordered_map<Node::Edge *, int>& minLengthEnd,
        int maxLength,
        DuplicateEdge duplicateEdge){
    unordered_set<long> visitedEdgeId;
    unordered_set<Node::Edge*> visitedEdge;
    bool useVisitedEdgeId = duplicateEdge == avoidDuplicateUndirectedEdges || duplicateEdge == unallowDuplicateUndirectedEdges;
    bool useVisitedEdge = duplicateEdge == avoidDuplicateDirectedEdges || duplicateEdge == unallowDuplicateDirectedEdges;
    bool avoid = useVisitedEdgeId || useVisitedEdge;
    PathEdges* newPathEdges;
    if (rand() & 1) {
        int targetCutIndex = rand() % (pathEdges->size() + 1);
        int index = 0;
        if (avoid){
            for (PathEdges::PathEdge* pathEdge: pathEdges->pathEdges){
                if (index == targetCutIndex
                    || (duplicateEdge == unallowDuplicateUndirectedEdges
                        && visitedEdgeId.find(pathEdge->edge->id) != visitedEdgeId.end())
                    || (duplicateEdge == unallowDuplicateDirectedEdges
                        && visitedEdge.find(pathEdge->edge) != visitedEdge.end())) {
                    break;
                }
                index++;
                if (useVisitedEdgeId){
                    visitedEdgeId.insert(pathEdge->edge->id);
                } else if (useVisitedEdge){
                    visitedEdge.insert(pathEdge->edge);
                }
            }
        }
        newPathEdges = pathEdges->cutAfter(index);

        Node::Edge *prevEdge = newPathEdges->lastEdge();
        while (true) {
            vector<Node::Edge*> nextEdgeCandidates;
            vector<Node::Edge*> nextEdgeCandidatesVisited;
            for (pair<const long, Node::Edge *> pair1: newPathEdges->getEnd()->edges) {
                Node::Edge *curEdge = pair1.second;

                if (Graph::isRestricted(prevEdge, curEdge, false)) {
                    continue;
                }

                auto it = minLengthEnd.find(curEdge);


                if (it == minLengthEnd.end() || newPathEdges->length + it->second > maxLength) {
                    continue;
                }

                if (duplicateEdge == ignoreDuplicateEdges
                    || (useVisitedEdgeId && visitedEdgeId.find(curEdge->id) == visitedEdgeId.end())
                    || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())){
                    nextEdgeCandidates.push_back(curEdge);
                } else if (avoid){
                    nextEdgeCandidatesVisited.push_back(curEdge);
                }
            }

            if (newPathEdges->getEnd() == pathEdges->getEnd()){
                break;
            }

            // even if current location is the target location, make moves if possible to maximize elevation gain
            if (nextEdgeCandidates.empty() && nextEdgeCandidatesVisited.empty()) {
                if (newPathEdges->getEnd() != pathEdges->getEnd()){
                    if (duplicateEdge != unallowDuplicateUndirectedEdges
                        && duplicateEdge != unallowDuplicateDirectedEdges){
                        throw runtime_error("error"); // temporal error statement
                    }
                    return pathEdges;
                }
                break;
            }

            // visit unvisited node first if not ignoring duplicate edges
            if (!nextEdgeCandidates.empty()) {
                prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
            } else {
                prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
            }

            newPathEdges->addEdge(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->id);
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
        }
    } else {
        int targetCutRIndex = rand() % (pathEdges->size() + 1);
        int rIndex = 0;
        vector<Node::Edge*> newEdges;

        if (avoid) {
            for (auto it=pathEdges->pathEdges.rbegin(); rIndex < targetCutRIndex ; it++){
                if (duplicateEdge == unallowDuplicateUndirectedEdges
                    && visitedEdgeId.find((*it)->edge->id) == visitedEdgeId.end()
                    || (duplicateEdge == unallowDuplicateDirectedEdges
                        && visitedEdge.find((*it)->edge) != visitedEdge.end())){
                    break;
                }
                rIndex++;
                if (useVisitedEdgeId){
                    visitedEdgeId.insert((*it)->edge->id);
                } else if (useVisitedEdge){
                    visitedEdge.insert((*it)->edge);
                }
            }
        }
        PathEdges::PathEdge *prevPathEdge = (rIndex == pathEdges->size()) ?
                                            nullptr :  pathEdges->at(pathEdges->size()-1-rIndex);
        Node::Edge *prevEdge = prevPathEdge == nullptr ? nullptr : prevPathEdge->edge;
        int length = prevPathEdge == nullptr ? 0 : pathEdges->length - prevPathEdge->length + prevEdge->length;
        Node *prevNode = prevPathEdge == nullptr ? pathEdges->end : prevEdge->u;
        while (true) {
            vector<Node::Edge*> nextEdgeCandidates;
            vector<Node::Edge*> nextEdgeCandidatesVisited;
            for (pair<const long, Node::Edge *> pair1: prevNode->reversedEdges) {
                Node::Edge *curEdge = pair1.second;

                if (Graph::isRestricted(prevEdge, curEdge, true)) {
                    continue;
                }

                auto it = minLengthStart.find(curEdge);
                if (it == minLengthStart.end() || length + it->second > maxLength) {
                    continue;
                }

                if (duplicateEdge == ignoreDuplicateEdges
                    || (useVisitedEdgeId && visitedEdgeId.find(curEdge->id) == visitedEdgeId.end())
                    || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())   ){
                    nextEdgeCandidates.push_back(curEdge);
                } else if (avoid){
                    nextEdgeCandidatesVisited.push_back(curEdge);
                }
            }

            if (!newEdges.empty() && newEdges.back()->u == pathEdges->getStart()){
                break;
            }

            // even if current location is the target location, make moves if possible to maximize elevation gain
            if (nextEdgeCandidates.empty() && nextEdgeCandidatesVisited.empty()) {
                if (!newEdges.empty() && newEdges.back()->u != pathEdges->getStart()) {
                    return pathEdges;
                }
                break;
            }

            // visit unvisited node first if minimizeRevisit == true
            if (!nextEdgeCandidates.empty()) {
                prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
            } else {
                prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
            }

            newEdges.push_back(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->id);
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
            length += prevEdge->length;
            prevNode = prevEdge->u;
        }


        newPathEdges = new PathEdges(pathEdges->start);
        for (auto it=newEdges.rbegin(); it!=newEdges.rend(); it++){
            newPathEdges->addEdge(*it);
        }
        if (rIndex < pathEdges->size()) {
            auto it = pathEdges->pathEdges.rbegin() + rIndex;
            while (true) {
                newPathEdges->addEdge((*it)->edge);
                if (it-- == pathEdges->pathEdges.rbegin()) {
                    break;
                }
            }
        }
    }

    if (duplicateEdge == unallowDuplicateUndirectedEdges) {
        newPathEdges = removeDuplicateUndirectedEdge(newPathEdges);
    } else if (duplicateEdge == unallowDuplicateDirectedEdges){
        newPathEdges = removeDuplicateDirectedEdge(newPathEdges);
    }
    return newPathEdges;
}

unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> naturalSelection(
        unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual>& pathEdgesSet, size_t numMaxSelect
        ){



    if (pathEdgesSet.size() < numMaxSelect){
        return pathEdgesSet;
    }
    priority_queue<PathEdges*, vector<PathEdges*>, PathEdges::biggerElevation> pq;
    for (auto it: pathEdgesSet){
        if (pq.size() >= numMaxSelect){
            if (it->elevation > pq.top()->elevation){
                pq.pop();
                pq.push(it);
            }
        } else {
            pq.push(it);
        }
    }
    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> res;
    while (!pq.empty()){
        res.insert(pq.top());
        pq.pop();
    }
    return res;
}

int weightedSelection(int numChoices, default_random_engine rng){
    int totalWeights = 0;
    for (int i=0; i<numChoices; i++){
        totalWeights += (i+1) * (i+1);
    }
    uniform_int_distribution<int> distribution(0, totalWeights-1);
    int r = distribution(rng);
    int i = 1;
    while (r >= 0){
        r -= i * i;
        i++;
    }
    return numChoices - i + 1;

}

PathEdges* elenaPathSearchMaxUsingGeneticAlgorithm(
        Graph *graph, Node *start, Node *end, double maxLengthRatio, size_t numProduce,
        size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch){

    auto t0 = chrono::high_resolution_clock::now();

    // 1) minWeightEnd and the shortest path

    int curMinWeight = INT_MAX;

    unordered_map<Node*, int> minWeightStart = dijkstraAlgorithm(
            graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
            false, nullptr, nullptr
    );

    unordered_map<Node *, Node::Edge *> prevEdgeMap;
    unordered_map<Node *, int> minWeightEnd = dijkstraAlgorithm(
            graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
            true, &prevEdgeMap, nullptr
    );

    // Find the shortest path
    PathEdges* pathEdges = new PathEdges(start);
    Node* lastNode = start;
    Node::Edge* lastEdge;
    while (lastNode != end) {
        lastEdge = prevEdgeMap[lastNode];
        pathEdges->addEdge(lastEdge);
        lastNode = lastEdge->v;
    }

    int maxWeight = getMaxWeight(curMinWeight, maxLengthRatio);

    // First generation: mutations from the shortest path
    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    for (size_t i=0; i<numProduce; i++){
        pathEdgesSet.insert(
                mutatePathEdges(pathEdges, minWeightStart, minWeightEnd, maxWeight, duplicateEdge));
    }

    for (int epoch=0; epoch < numEpoch; epoch++){
//
//        cout << "====================================" <<  endl;
//        for (PathEdges* pathEdges: pathEdgesSet){
//            cout << pathEdges->elevation << " " << pathEdges->length << " " << maxWeight << endl;
//        }
//        cout << "====================================;" <<  endl;

        // natural selection
        pathEdgesSet = naturalSelection(pathEdgesSet, numMaxSelect);


        // prepare for crossovers: search for intersections
        unordered_map<Node::Edge*, unordered_map<PathEdges*, vector<size_t>>> intersections;
        for (auto it: pathEdgesSet){
            size_t index = 0;
            for (PathEdges::PathEdge* pathEdge: it->pathEdges){
                intersections[pathEdge->edge][it].push_back(index++);
            }
        }

        for (auto it=intersections.begin(); it!=intersections.end();){
            if (it->second.size() < 2){
                it = intersections.erase(it);
            } else {
                it++;
            }
        }

        // 4) crossovers: based on intersections find valid paths, generate mutants of each new crossovers with m% probability
        // Choose an intersected edge
        for (int produce=0; produce<numProduce; produce++) {
            auto intersectionsIt = intersections.begin();
            if (intersections.size() == 0){
                auto it = pathEdgesSet.begin();
                advance(it, rand() % pathEdgesSet.size());
                pathEdgesSet.insert(
                        mutatePathEdges(*it, minWeightStart, minWeightEnd, maxWeight, duplicateEdge));
                continue;
            }

            advance(intersectionsIt, (rand() % intersections.size()));
            unordered_map<PathEdges *, vector<size_t>> &choices = intersectionsIt->second;

            // choose two paths
            int index1(rand() % intersectionsIt->second.size()), index2(rand() % (intersectionsIt->second.size() - 1));
            if (index1 == index2)
                index2++;
            auto choicesIt1(choices.begin()), choicesIt2(choices.begin());
            advance(choicesIt1, index1);
            advance(choicesIt2, index2);
            PathEdges *pathEdges1 = choicesIt1->first;
            PathEdges *pathEdges2 = choicesIt2->first;

            // use the second choice as the prefix if the first choice does not give any path within the length limit
            if (pathEdges1->at(choicesIt1->second.front())->length
                + pathEdges2->length - pathEdges2->at(choicesIt2->second.back())->length > maxWeight) {
                swap(pathEdges1, pathEdges2);
                swap(choicesIt1, choicesIt2);
            }
            vector<size_t>& indices1 = choicesIt1->second;
            vector<size_t>& indices2 = choicesIt2->second;

//            // find the maximum value for index1 and randomly choose one
//            int maxIndex1 = 0;
//            for (int i=1; i<indices1.size(); i++){
//                if (pathEdges1->at(indices1[i])->length + pathEdges2->length - pathEdges2->at(indices2.back())->length > maxWeight)
//                    break;
//                maxIndex1 = i;
//            }
//            index1 = rand() % (maxIndex1 + 1);
////            index1 = weightedSelection(maxIndex1 + 1, graph->rng);
//
//            // find the minimum value for index2 and randomly choose one
//            int minIndex2 = indices2.size()-1;
//            for (int i=indices2.size()-2; i>=0; i--){
//                if (pathEdges1->at(indices1[index1])->length + pathEdges2->length - pathEdges2->at(indices2[i])->length > maxWeight)
//                    break;
//                minIndex2 = i;
//            }
//            index2 = (rand() % (int) (indices2.size() - minIndex2)) + minIndex2;

//            index2 = indices2.size()-1 - weightedSelection(indices2.size()-minIndex2, graph->rng);

            index1 = indices1.front();
            index2 = indices2.back();

            // concatenate two paths
            // duplicate edges will be handled in the mutation step
            pathEdges = pathEdges1->cutAfter(index1);
            if (index2 + 1 < pathEdges2->size()) {

                // check if restricted
                Node::Edge* edge = pathEdges1->at(index1)->edge;
                if (edge->v->restrictions != nullptr) {
                    auto restrictionsIt = edge->v->restrictions->find(edge);
                    if (restrictionsIt != edge->v->restrictions->end()) {
                        if (restrictionsIt->second.find(pathEdges2->at(index2 + 1)->edge)
                            != restrictionsIt->second.end()) {
                            continue;
                        }
                    }
                }

                auto pathEdgesIt = pathEdges2->pathEdges.begin();

                advance(pathEdgesIt, index2 + 1);
                while (pathEdgesIt != pathEdges2->pathEdges.end()) {
                    pathEdges->addEdge((*pathEdgesIt)->edge);
                    pathEdgesIt++;
                }
            }



            pathEdgesSet.insert(
                    mutatePathEdges(pathEdges, minWeightStart, minWeightEnd, maxWeight, duplicateEdge));
        }
        if (epoch >= minEpoch && chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - t0).count() > maxMilliseconds){
            break;
        }
    }

    int maxElevation = -1.0;
    for (auto it: pathEdgesSet){
        if (it->getElevation() > maxElevation){
            maxElevation = it->getElevation();
            pathEdges = it;
        }
    }
    return pathEdges;

}

PathEdges* elenaPathSearchMaxUsingEdgeBasedGeneticAlgorithm(
        Graph *graph, Node *start, Node *end, double maxLengthRatio, size_t numProduce,
        size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch){

    auto t0 = chrono::high_resolution_clock::now();
    // 1) minWeightEnd and the shortest path

    int curMinWeight = INT_MAX;
    Node::Edge *lastEdge = nullptr;

    unordered_map<Node::Edge *, int> minWeightStart = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
            false, nullptr, lastEdge, nullptr
    );

    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMap;
    unordered_map<Node::Edge *, int> minWeightEnd = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength,
            true, &prevEdgeMap, lastEdge, nullptr
    );
    if (lastEdge == nullptr){
        throw runtime_error("path not found.");
    }

    // Find the shortest path
    PathEdges* pathEdges = new PathEdges(start);
    while (lastEdge->id != -1) {
        pathEdges->addEdge(lastEdge);
        lastEdge = prevEdgeMap[lastEdge];
    }

    int maxWeight = getMaxWeight(curMinWeight, maxLengthRatio);

    // First generation: mutations from the shortest path
    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    for (size_t i=0; i<numProduce; i++){
        pathEdgesSet.insert(
                edgeBasedMutatePathEdges(pathEdges, minWeightStart, minWeightEnd, maxWeight, duplicateEdge));
    }

    for (int epoch=0; epoch < numEpoch; epoch++){

        // natural selection
        pathEdgesSet = naturalSelection(pathEdgesSet, numMaxSelect);


        // prepare for crossovers: search for intersections
        unordered_map<Node::Edge*, unordered_map<PathEdges*, vector<size_t>>> intersections;
        for (auto it: pathEdgesSet){
            size_t index = 0;
            for (PathEdges::PathEdge* pathEdge: it->pathEdges){
                intersections[pathEdge->edge][it].push_back(index++);
            }
        }

        for (auto it=intersections.begin(); it!=intersections.end();){
            if (it->second.size() < 2){
                it = intersections.erase(it);
            } else {
                it++;
            }
        }

        // 4) crossovers: based on intersections find valid paths, generate mutants of each new crossovers
        // Choose an intersected edge
        for (int produce=0; produce<numProduce; produce++) {
            auto intersectionsIt = intersections.begin();
            if (intersections.size() == 0 || rand() % 4 < 2){
                auto it = pathEdgesSet.begin();
                advance(it, rand() % pathEdgesSet.size());
                pathEdgesSet.insert(
                        edgeBasedMutatePathEdges(*it, minWeightStart, minWeightEnd, maxWeight, duplicateEdge));
                continue;
            }

            advance(intersectionsIt, (rand() % intersections.size()));
            unordered_map<PathEdges *, vector<size_t>> &choices = intersectionsIt->second;

            // choose two paths
            int index1(rand() % intersectionsIt->second.size()), index2(rand() % (intersectionsIt->second.size() - 1));
            if (index1 == index2)
                index2++;
            auto choicesIt1(choices.begin()), choicesIt2(choices.begin());
            advance(choicesIt1, index1);
            advance(choicesIt2, index2);
            PathEdges *pathEdges1 = choicesIt1->first;
            PathEdges *pathEdges2 = choicesIt2->first;

            // use the second choice as the prefix if the first choice does not give any path within the length limit
            if (pathEdges1->at(choicesIt1->second.front())->length
                + pathEdges2->length - pathEdges2->at(choicesIt2->second.back())->length > maxWeight) {
                swap(pathEdges1, pathEdges2);
                swap(choicesIt1, choicesIt2);
            }
            vector<size_t>& indices1 = choicesIt1->second;
            vector<size_t>& indices2 = choicesIt2->second;

//            // find the maximum value for index1 and randomly choose one
//            int maxIndex1 = 0;
//            for (int i=1; i<indices1.size(); i++){
//                if (pathEdges1->at(indices1[i])->length + pathEdges2->length - pathEdges2->at(indices2.back())->length > maxWeight)
//                    break;
//                maxIndex1 = i;
//            }
//            index1 = rand() % (maxIndex1 + 1);
////            index1 = weightedSelection(maxIndex1 + 1, graph->rng);
//
//            // find the minimum value for index2 and randomly choose one
//            int minIndex2 = indices2.size()-1;
//            for (int i=indices2.size()-2; i>=0; i--){
//                if (pathEdges1->at(indices1[index1])->length + pathEdges2->length - pathEdges2->at(indices2[i])->length > maxWeight)
//                    break;
//                minIndex2 = i;
//            }
//            index2 = (rand() % (int) (indices2.size() - minIndex2)) + minIndex2;

//            index2 = indices2.size()-1 - weightedSelection(indices2.size()-minIndex2, graph->rng);

            index1 = indices1.front();
            index2 = indices2.back();

            // concatenate two paths
            // duplicate edges will be handled in the mutation step
            pathEdges = pathEdges1->cutAfter(index1);
            if (index2 + 1 < pathEdges2->size()) {

                // check if restricted
                Node::Edge* edge = pathEdges1->at(index1)->edge;
                if (edge->v->restrictions != nullptr) {
                    auto restrictionsIt = edge->v->restrictions->find(edge);
                    if (restrictionsIt != edge->v->restrictions->end()) {
                        if (restrictionsIt->second.find(pathEdges2->at(index2 + 1)->edge)
                            != restrictionsIt->second.end()) {
                            continue;
                        }
                    }
                }

                auto pathEdgesIt = pathEdges2->pathEdges.begin();

                advance(pathEdgesIt, index2 + 1);
                while (pathEdgesIt != pathEdges2->pathEdges.end()) {
                    pathEdges->addEdge((*pathEdgesIt)->edge);
                    pathEdgesIt++;
                }
            }

            pathEdgesSet.insert(
                    edgeBasedMutatePathEdges(pathEdges, minWeightStart, minWeightEnd, maxWeight, duplicateEdge));
        }

        if (epoch >= minEpoch && chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - t0).count() > maxMilliseconds){
            break;
        }
    }

    int maxElevation = -1.0;
    for (auto it: pathEdgesSet){
        if (it->getElevation() > maxElevation){
            maxElevation = it->getElevation();
            pathEdges = it;
        }
    }
    return pathEdges;

}
