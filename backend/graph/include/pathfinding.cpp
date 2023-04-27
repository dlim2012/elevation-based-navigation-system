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
        throw invalid_argument("At least one node does not belong to the graph");
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
            if (validNodes != nullptr && validNodes->find(curNode) == validNodes->end()) {
                continue;
            }

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
                if ((reversed && curNode == start)
                    || (!reversed && curNode == end)) {
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
        throw invalid_argument("At least one node does not belong to the graph");
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
                if ((reversed && curEdge->u == start)
                    || (!reversed && curEdge->v == end)) {
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

Path* pathFromPrevEdgeMap(
        Node* start,
        Node* end,
        unordered_map<Node*, Node::Edge*>& prevEdgeMap
        ){

    Path* path = new Path();
    Node::Edge* prevEdge;
    while (start != end) {
        prevEdge = prevEdgeMap[start];
        path->addEdge(prevEdge);
        start = prevEdge->v;
    }
    return path;
}

Path *dijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        int getWeight(Node::Edge *edge),
        unordered_set<Node *> *validNodes
        ){
    if (start == end) {
        return new Path();
    }
    unordered_map<Node *, Node::Edge *> prevEdgeMap;
    int curMinWeight = INT_MAX;

    auto t0 = chrono::high_resolution_clock::now();
    // find the path reversed so that the path can be reconstructed non-reversed
    unordered_map<Node *, int> minWeightEdge = dijkstraAlgorithm(
            graph, start, end, curMinWeight, -1.0, getWeight, true, &prevEdgeMap, validNodes);
    auto t1 = chrono::high_resolution_clock::now();
    auto t = (double) chrono::duration_cast<chrono::microseconds>(t1 - t0).count() / 1000000;

    return pathFromPrevEdgeMap(start, end, prevEdgeMap);
}

Path* edgeBasedPathFromPrevEdgeMap(
        Node::Edge* lastEdge,
        unordered_map<Node::Edge *, Node::Edge *>& prevEdgeMap
        ){
    Path* path = new Path();

    while (lastEdge->id != -1) {
        path->addEdge(lastEdge);
        lastEdge = prevEdgeMap[lastEdge];
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
    if (start == end) {
        return new Path();
    }
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMap;
    Node::Edge *lastEdge = nullptr;
    int curMinWeight = INT_MAX;
    // find the path reversed so that the path can be reconstructed non-reversed
    unordered_map<Node::Edge *, int> minWeightEdge = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, -1.0, getWeight, true, &prevEdgeMap, lastEdge, validEdges);


    if (lastEdge == nullptr){
        throw runtime_error("path not found.");
    }
    return edgeBasedPathFromPrevEdgeMap(lastEdge, prevEdgeMap);
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
    if (minWeightStart.size() == graph->edgeCount + 1
        && max_element(minWeightStart.begin(), minWeightStart.end(), [] (
                const pair<Node*, int> p1, const pair<Node*, int> p2){
            return p1.second < p2.second;
        })->second <= curMinWeight * maxWeightRatio){
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
    if (minWeightStart.size() == graph->edgeCount + 1
        && max_element(minWeightStart.begin(), minWeightStart.end(), [] (
            const pair<Node::Edge*, int> p1, const pair<Node::Edge*, int> p2){
        return p1.second < p2.second;
    })->second <= curMinWeight * maxWeightRatio){
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


Path *findMinElevationUsingDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    unordered_set<Node *> validNodes = getReachableNodes(graph, start, end, maxWeightRatio);
    return dijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, validNodes.empty() ? nullptr : &validNodes);
}

Path *findMinElevationUsingEdgeBasedDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    unordered_set<Node::Edge *> validEdges = getReachableEdges(graph, start, end, maxWeightRatio);
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
        Node* start,
        Node* end,
        unordered_map<Node *, int>& minLengthStart,
        unordered_map<Node *, int>& minLengthEnd,
        unordered_map<Node *, Node::Edge* >& prevEdgeMapMinElevationStart,
        unordered_map<Node *, Node::Edge *>& prevEdgeMapMinElevationEnd,
        int maxLength,
        DuplicateEdge duplicateEdge,
        bool maximize
){
    unordered_set<long> visitedEdgeId;
    unordered_set<Node::Edge*> visitedEdge;
    bool useVisitedEdgeId = duplicateEdge == avoidDuplicateUndirectedEdges || duplicateEdge == unallowDuplicateUndirectedEdges;
    bool useVisitedEdge = duplicateEdge == avoidDuplicateDirectedEdges || duplicateEdge == unallowDuplicateDirectedEdges;
    bool useVisited = useVisitedEdgeId || useVisitedEdge;

    // for minimization: choose first edge randomly, then, follow the min elevation path if possible
    bool firstEdge = true;
    PathEdges* newPathEdges;
    if (rand() & 1) {
        if (pathEdges != nullptr) {
            int targetCutIndex = rand() % (pathEdges->size() + 1);
            int index = 0;
            if (useVisited) {
                for (PathEdges::PathEdge *pathEdge: pathEdges->pathEdges) {
                    if (index == targetCutIndex
                        || (duplicateEdge == unallowDuplicateUndirectedEdges
                            && visitedEdgeId.find(pathEdge->edge->id) != visitedEdgeId.end())
                        || (duplicateEdge == unallowDuplicateDirectedEdges
                            && visitedEdge.find(pathEdge->edge) != visitedEdge.end())) {
                        break;
                    }
                    index++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert(pathEdge->edge->id);
                    } else if (useVisitedEdge) {
                        visitedEdge.insert(pathEdge->edge);
                    }
                }
            }
            newPathEdges = pathEdges->cutBefore(index);
        } else {
            newPathEdges = new PathEdges(start);
        }

        Node::Edge *prevEdge = newPathEdges->lastEdge();
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge){
                    firstEdge = false;
                } else {
                    Node::Edge *curEdge = prevEdgeMapMinElevationEnd[newPathEdges->getEnd()];
                    if (curEdge != nullptr) {
                        auto it = minLengthEnd.find(curEdge->v);
                        if (it != minLengthEnd.end() &&
                            newPathEdges->length + curEdge->length + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge){
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
                    } else if (useVisited){
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }
                }

                if (newPathEdges->getEnd() == end){
                    break;
                }

                // visit unvisited node first if not ignoring duplicate edges
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
            }

            newPathEdges->addEdge(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->id);
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }

        }
    } else {
        vector<Node::Edge *> newEdges;
        PathEdges::PathEdge *prevPathEdge;
        int rIndex = 0;
        if (pathEdges != nullptr) {
            int targetCutRIndex = rand() % (pathEdges->size() + 1);

            if (useVisited) {
                for (auto it = pathEdges->pathEdges.rbegin(); rIndex < targetCutRIndex; it++) {
                    if (duplicateEdge == unallowDuplicateUndirectedEdges
                        && visitedEdgeId.find((*it)->edge->id) == visitedEdgeId.end()
                        || (duplicateEdge == unallowDuplicateDirectedEdges
                            && visitedEdge.find((*it)->edge) != visitedEdge.end())) {
                        break;
                    }
                    rIndex++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert((*it)->edge->id);
                    } else if (useVisitedEdge) {
                        visitedEdge.insert((*it)->edge);
                    }
                }
            }
            prevPathEdge = (rIndex == pathEdges->size()) ?
                                                nullptr : pathEdges->at(pathEdges->size() - 1 - rIndex);
        } else {
            prevPathEdge = nullptr;
        }

        Node::Edge *prevEdge = prevPathEdge == nullptr ? nullptr : prevPathEdge->edge;
        int length = prevPathEdge == nullptr ? 0 : pathEdges->length - prevPathEdge->length + prevEdge->length;
        Node *prevNode = prevPathEdge == nullptr ? end : prevEdge->u;
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge){
                    firstEdge = false;
                } else {
                    Node::Edge *curEdge = prevEdgeMapMinElevationStart[prevNode];
                    if (curEdge != nullptr) {
                        auto it = minLengthStart.find(curEdge->u);
                        if (it != minLengthStart.end() && length + curEdge->length + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge) {
                vector<Node::Edge *> nextEdgeCandidates;
                vector<Node::Edge *> nextEdgeCandidatesVisited;
                for (pair<const long, Node::Edge *> pair1: prevNode->reversedEdges) {
                    Node::Edge *curEdge = pair1.second;

                    auto it = minLengthStart.find(curEdge->u);
                    if (it == minLengthStart.end() || length + +curEdge->length + it->second > maxLength) {
                        continue;
                    }

                    if (duplicateEdge == ignoreDuplicateEdges
                        || (useVisitedEdgeId && visitedEdgeId.find(curEdge->id) == visitedEdgeId.end())
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())) {
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited) {
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }

                }

                if (!newEdges.empty() && newEdges.back()->u == start) {
                    break;
                }

                // visit unvisited node first if minimizeRevisit == true
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
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


        newPathEdges = new PathEdges(start);
        for (auto it=newEdges.rbegin(); it!=newEdges.rend(); it++){
            newPathEdges->addEdge(*it);
        }
        if (pathEdges != nullptr) {
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
        Node* start,
        Node* end,
        unordered_map<Node::Edge *, int>& minLengthStart,
        unordered_map<Node::Edge *, int>& minLengthEnd,
        unordered_map<Node::Edge *, Node::Edge* >& prevEdgeMapMinElevationStart,
        unordered_map<Node::Edge *, Node::Edge *>& prevEdgeMapMinElevationEnd,
        int maxLength,
        DuplicateEdge duplicateEdge,
        bool maximize
        ){
    unordered_set<long> visitedEdgeId;
    unordered_set<Node::Edge*> visitedEdge;
    bool useVisitedEdgeId = duplicateEdge == avoidDuplicateUndirectedEdges || duplicateEdge == unallowDuplicateUndirectedEdges;
    bool useVisitedEdge = duplicateEdge == avoidDuplicateDirectedEdges || duplicateEdge == unallowDuplicateDirectedEdges;
    bool useVisited = useVisitedEdgeId || useVisitedEdge;

    PathEdges* newPathEdges;
    bool firstEdge = true;
    if (rand() & 1) {
        if (pathEdges != nullptr) {
            int targetCutIndex = rand() % (pathEdges->size() + 1);
            int index = 0;
            if (useVisited) {
                for (PathEdges::PathEdge *pathEdge: pathEdges->pathEdges) {
                    if (index == targetCutIndex
                        || (duplicateEdge == unallowDuplicateUndirectedEdges
                            && visitedEdgeId.find(pathEdge->edge->id) != visitedEdgeId.end())
                        || (duplicateEdge == unallowDuplicateDirectedEdges
                            && visitedEdge.find(pathEdge->edge) != visitedEdge.end())) {
                        break;
                    }
                    index++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert(pathEdge->edge->id);
                    } else if (useVisitedEdge) {
                        visitedEdge.insert(pathEdge->edge);
                    }
                }
            }
            newPathEdges = pathEdges->cutBefore(index);
        } else {
            newPathEdges = new PathEdges(start);
        }

        Node::Edge *prevEdge = newPathEdges->lastEdge();
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge) {
                    firstEdge = false;
                } else {
                    Node::Edge *curEdge = prevEdgeMapMinElevationEnd[newPathEdges->lastEdge()];
                    if (curEdge != nullptr){
                        auto it = minLengthEnd.find(curEdge);
                        if (it != minLengthEnd.end() && newPathEdges->length + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge) {
                vector<Node::Edge *> nextEdgeCandidates;
                vector<Node::Edge *> nextEdgeCandidatesVisited;
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
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())) {
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited) {
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }
                }

                if (newPathEdges->getEnd() == end) {
                    break;
                }

                // visit unvisited node first if not ignoring duplicate edges
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
            }

            newPathEdges->addEdge(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->id);
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
        }
    } else {
        int rIndex = 0;
        vector<Node::Edge*> newEdges;
        PathEdges::PathEdge *prevPathEdge;
        if (pathEdges != nullptr) {
            int targetCutRIndex = rand() % (pathEdges->size() + 1);

            if (useVisited) {
                for (auto it = pathEdges->pathEdges.rbegin(); rIndex < targetCutRIndex; it++) {
                    if (duplicateEdge == unallowDuplicateUndirectedEdges
                        && visitedEdgeId.find((*it)->edge->id) == visitedEdgeId.end()
                        || (duplicateEdge == unallowDuplicateDirectedEdges
                            && visitedEdge.find((*it)->edge) != visitedEdge.end())) {
                        break;
                    }
                    rIndex++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert((*it)->edge->id);
                    } else if (useVisitedEdge) {
                        visitedEdge.insert((*it)->edge);
                    }
                }
            }
            prevPathEdge = (rIndex == pathEdges->size()) ?
                            nullptr : pathEdges->at(pathEdges->size() - 1 - rIndex);
        } else {
            prevPathEdge = nullptr;
        }
        Node::Edge *prevEdge = prevPathEdge == nullptr ? nullptr : prevPathEdge->edge;
        int length = prevPathEdge == nullptr ? 0 : pathEdges->length - prevPathEdge->length + prevEdge->length;
        Node *prevNode = prevPathEdge == nullptr ? end : prevEdge->u;
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge){
                    firstEdge = false;
                } else {
                    Node::Edge *curEdge = prevEdgeMapMinElevationStart[prevEdge];
                    if (curEdge != nullptr) {
                        auto it = minLengthStart.find(curEdge);
                        if (it != minLengthStart.end() && length + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge) {
                vector<Node::Edge *> nextEdgeCandidates;
                vector<Node::Edge *> nextEdgeCandidatesVisited;
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
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())) {
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited) {
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }
                }

                if (!newEdges.empty() && newEdges.back()->u == start) {
                    break;
                }

                // visit unvisited node first if minimizeRevisit == true
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
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


        newPathEdges = new PathEdges(start);
        for (auto it=newEdges.rbegin(); it!=newEdges.rend(); it++){
            newPathEdges->addEdge(*it);
        }
        if (pathEdges != nullptr) {
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
    }

    if (duplicateEdge == unallowDuplicateUndirectedEdges) {
        newPathEdges = removeDuplicateUndirectedEdge(newPathEdges);
    } else if (duplicateEdge == unallowDuplicateDirectedEdges){
        newPathEdges = removeDuplicateDirectedEdge(newPathEdges);
    }
    return newPathEdges;
}



unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> naturalSelection(
        unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual>& pathEdgesSet, size_t numMaxSelect, bool maximize
        ){

    if (pathEdgesSet.size() < numMaxSelect){
        return pathEdgesSet;
    }

    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> res;
    if (maximize) {
        priority_queue<PathEdges *, vector<PathEdges *>, PathEdges::biggerElevation> pq;
        for (auto it: pathEdgesSet) {
            if (pq.size() >= numMaxSelect) {
                if (it->elevation > pq.top()->elevation) {
                    pq.pop();
                    pq.push(it);
                }
            } else {
                pq.push(it);
            }
        }
        while (!pq.empty()){
            res.insert(pq.top());
            pq.pop();
        }
    } else {
        priority_queue<PathEdges *, vector<PathEdges *>, PathEdges::smallerElevation> pq;
        for (auto it: pathEdgesSet) {
            if (pq.size() >= numMaxSelect) {
                if (it->elevation < pq.top()->elevation) {
                    pq.pop();
                    pq.push(it);
                }
            } else {
                pq.push(it);
            }
        }
        while (!pq.empty()){
            res.insert(pq.top());
            pq.pop();
        }
    }

    return res;
}


Path* geneticAlgorithm(
        Graph *graph, Node *start, Node *end, int maxLength,
        unordered_map<Node*, int>& minLengthStart, unordered_map<Node*, int>& minLengthEnd,
        PathEdges* shortestPathEdges,
        size_t numProduce,
        size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch, bool maximize) {

    auto t0 = chrono::high_resolution_clock::now();

    PathEdges *pathEdges;
    unordered_map<Node *, Node::Edge* > prevEdgeMapMinElevationStart;
    unordered_map<Node *, Node::Edge *> prevEdgeMapMinElevationEnd;
    if (!maximize){
        unordered_set<Node*> validNodes;
        for (pair<Node *const, int> &pair1: minLengthEnd) {
            auto it = minLengthStart.find(pair1.first);
            if (it == minLengthStart.end())
                continue;
            if (pair1.second + it->second > maxLength)
                continue;
            validNodes.insert(pair1.first);
        }
        int curMinWeight = INT_MAX;

        auto t0 = chrono::high_resolution_clock::now();
        // find the path reversed so that the path can be reconstructed non-reversed
        dijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0, Node::Edge::getElevation, true, &prevEdgeMapMinElevationEnd, &validNodes);
        Path* dijkstraPath = pathFromPrevEdgeMap(start, end, prevEdgeMapMinElevationEnd);
        if (dijkstraPath->length <= maxLength){
            return pathFromPrevEdgeMap(start, end, prevEdgeMapMinElevationEnd);
        }
        dijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0, Node::Edge::getElevation, false, &prevEdgeMapMinElevationStart, &validNodes);
    }

    // first generation: optional shortestPathEdges and random paths
    unordered_set<PathEdges *, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    if (shortestPathEdges != nullptr){
        pathEdgesSet.insert(shortestPathEdges);
    }
    for (size_t i = 0; i < numProduce; i++) {
        pathEdgesSet.insert(
                mutatePathEdges(nullptr, start, end, minLengthStart, minLengthEnd,
                                prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd, maxLength, duplicateEdge, maximize));
    }

    for (int epoch = 0; epoch < numEpoch; epoch++) {

        // natural selection
        pathEdgesSet = naturalSelection(pathEdgesSet, numMaxSelect, maximize);


        // prepare for crossovers: search for intersections
        unordered_map<Node::Edge *, unordered_map<PathEdges *, vector<size_t>>> intersections;
        for (auto it: pathEdgesSet) {
            size_t index = 0;
            for (PathEdges::PathEdge *pathEdge: it->pathEdges) {
                intersections[pathEdge->edge][it].push_back(index++);
            }
        }

        for (auto it = intersections.begin(); it != intersections.end();) {
            if (it->second.size() < 2) {
                it = intersections.erase(it);
            } else {
                it++;
            }
        }

        // 4) crossovers: based on intersections find valid paths, generate mutants of each new crossovers with m% probability
        // Choose an intersected edge
        for (int produce = 0; produce < numProduce; produce++) {
            auto intersectionsIt = intersections.begin();
            if (intersections.size() == 0) {
                auto it = pathEdgesSet.begin();
                advance(it, rand() % pathEdgesSet.size());
                pathEdgesSet.insert(
                        mutatePathEdges(*it, start, end, minLengthStart, minLengthEnd,
                                        prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                        maxLength, duplicateEdge, maximize));
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
                + pathEdges2->length - pathEdges2->at(choicesIt2->second.back())->length > maxLength) {
                swap(pathEdges1, pathEdges2);
                swap(choicesIt1, choicesIt2);
            }
            vector<size_t> &indices1 = choicesIt1->second;
            vector<size_t> &indices2 = choicesIt2->second;

            index1 = indices1.front();
            index2 = indices2.back();

            // concatenate two paths
            // duplicate edges will be handled in the mutation step
            pathEdges = pathEdges1->cutAfter(index1);
            if (index2 + 1 < pathEdges2->size()) {

                // check if restricted
                Node::Edge *edge = pathEdges1->at(index1)->edge;
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

            if (maximize || rand() & 1) {
                pathEdgesSet.insert(
                        mutatePathEdges(pathEdges, start, end, minLengthStart, minLengthEnd,
                                        prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                        maxLength, duplicateEdge, maximize));
            } else {
                pathEdgesSet.insert(
                        pathEdges);
            }

        }

        if (epoch >= minEpoch && chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - t0).count() > maxMilliseconds) {
            break;
        }
    }

    if (maximize) {
        int maxElevation = -1;
        for (auto it: pathEdgesSet) {
            if (it->getElevation() > maxElevation) {
                maxElevation = it->getElevation();
                pathEdges = it;
            }
        }
    } else {
        int minElevation = INT_MAX;
        for (auto it: pathEdgesSet) {
            if (it->getElevation() < minElevation) {
                minElevation = it->getElevation();
                pathEdges = it;
            }
        }
    }
    return pathEdges->toPath();

}

Path* edgeBasedGeneticAlgorithm(
        Graph *graph, Node *start, Node *end, int maxLength,
        unordered_map<Node::Edge*, int>& minLengthStart, unordered_map<Node::Edge*, int>& minLengthEnd,
        PathEdges* shortestPathEdges,
        size_t numProduce, size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch, bool maximize
){

    auto t0 = chrono::high_resolution_clock::now();

    PathEdges* pathEdges;
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMapMinElevationStart;
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMapMinElevationEnd;

    if (!maximize){
        unordered_set<Node::Edge*> validEdges;
        for (pair<Node::Edge *const, int> &pair1: minLengthEnd) {
            auto it = minLengthStart.find(pair1.first);
            if (it == minLengthStart.end())
                continue;
            if (pair1.second + it->second > maxLength)
                continue;
            validEdges.insert(pair1.first);
        }

        Node::Edge *lastEdge = nullptr;
        int curMinWeight = INT_MAX;
        // find the path reversed so that the path can be reconstructed non-reversed
        edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0, Node::Edge::getElevation, true, &prevEdgeMapMinElevationEnd, lastEdge, &validEdges);
        Path* dijkstraPath = edgeBasedPathFromPrevEdgeMap(lastEdge, prevEdgeMapMinElevationEnd);
        if (dijkstraPath->length < maxLength){
            return dijkstraPath;
        }
        edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0, Node::Edge::getElevation, false, &prevEdgeMapMinElevationStart, lastEdge, &validEdges);
    }

    // First generation
    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    if (shortestPathEdges != nullptr) {
        pathEdgesSet.insert(shortestPathEdges);
    }
    for (size_t i=0; i<numProduce; i++){
        pathEdgesSet.insert(
                edgeBasedMutatePathEdges(nullptr, start, end, minLengthStart, minLengthEnd,
                                         prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                         maxLength, duplicateEdge, maximize)
                );
    }

    for (int epoch=0; epoch < numEpoch; epoch++){

        // natural selection
        pathEdgesSet = naturalSelection(pathEdgesSet, numMaxSelect, maximize);


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
            if (intersections.size() == 0){
                auto it = pathEdgesSet.begin();
                advance(it, rand() % pathEdgesSet.size());
                pathEdgesSet.insert(
                        edgeBasedMutatePathEdges(*it, start, end, minLengthStart, minLengthEnd,
                                                 prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                 maxLength, duplicateEdge, maximize));
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
                + pathEdges2->length - pathEdges2->at(choicesIt2->second.back())->length > maxLength) {
                swap(pathEdges1, pathEdges2);
                swap(choicesIt1, choicesIt2);
            }
            vector<size_t>& indices1 = choicesIt1->second;
            vector<size_t>& indices2 = choicesIt2->second;

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

            if (maximize || rand() & 1) {
                pathEdgesSet.insert(
                        edgeBasedMutatePathEdges(pathEdges, start, end, minLengthStart, minLengthEnd,
                                                 prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                 maxLength, duplicateEdge, maximize));
           } else {
                pathEdgesSet.insert(pathEdges);
            }
        }

        if (epoch >= minEpoch && chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - t0).count() > maxMilliseconds){
            break;
        }
    }


    if (maximize) {
        int maxElevation = -1;
        for (auto it: pathEdgesSet) {
            if (it->getElevation() > maxElevation) {
                maxElevation = it->getElevation();
                pathEdges = it;
            }
        }
    } else {
        int minElevation = INT_MAX;
        for (auto it: pathEdgesSet) {
            if (it->getElevation() < minElevation) {
                minElevation = it->getElevation();
                pathEdges = it;
            }
        }
    }
    return pathEdges->toPath();

}