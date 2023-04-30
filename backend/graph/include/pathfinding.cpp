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

    if (!graph->hasNode(start) || !graph->hasNode(end)) {
        throw invalid_argument("start and/or end nodes do not belong to the graph");
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

        if (visited.find(candidate.prevNode) != visited.end()) {
            continue;
        }
        visited.insert(candidate.prevNode);

        for (pair<const long, Node::Edge *> &pair1: reversed ? candidate.prevNode->getReversedEdges()
                                                             : candidate.prevNode->getEdges()) {
            curNode = reversed ? pair1.second->getU() : pair1.second->getV();
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


    if (!graph->hasNode(start) || !graph->hasNode(end)) {
        throw invalid_argument("start and/or end nodes do not belong to the graph");
    }

    if (start == end){
        return unordered_map<Node::Edge *, int>();
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

    Node::Edge *startEdge = new Node::Edge(reversed ? end : start, reversed);
    Node::Edge *curEdge = startEdge;
    minWeightEdge[curEdge] = 0;
    pq.emplace(0, curEdge);
    int curMaxWeight = getMaxWeight(curMinWeight, maxWeightRatio);

    while (!pq.empty()) {
        Candidate candidate = pq.top();
        pq.pop();

        if (visited.find(candidate.prevEdge) != visited.end()) {
            continue;
        }
        visited.insert(candidate.prevEdge);

        for (pair<const long, Node::Edge *> &pair1: reversed ? candidate.prevEdge->getU()->getReversedEdges()
                                                             : candidate.prevEdge->getV()->getEdges()) {
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
                if ((reversed && curEdge->getU() == start)
                    || (!reversed && curEdge->getV() == end)) {
                        if (curMinWeight >= curWeight) {
                        curMinWeight = curWeight;
                        curMaxWeight = getMaxWeight(curMinWeight, maxWeightRatio);
                        lastEdge = curEdge;
                    }
                }
            }
        }
    }

    minWeightEdge.erase(startEdge);
    delete startEdge;
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
        start = prevEdge->getV();
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

    Path* path = pathFromPrevEdgeMap(start, end, prevEdgeMap);

    return path;
}

Path* edgeBasedPathFromPrevEdgeMap(
        Node::Edge* lastEdge,
        unordered_map<Node::Edge *, Node::Edge *>& prevEdgeMap
        ){
    Path* path = new Path();

    if (lastEdge == nullptr){
        return path;
    }

    while (true) {
        auto it = prevEdgeMap.find(lastEdge);
        if (it == prevEdgeMap.end()){
            break;
        }
        path->addEdge(lastEdge);
        lastEdge = it->second;
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

    Path *path = edgeBasedPathFromPrevEdgeMap(lastEdge, prevEdgeMap);

    // clear memory
    unordered_map<Node::Edge *, Node::Edge *>().swap(prevEdgeMap);

    return path;
}

unordered_set<Node*> getReachableNodes(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
        ){


    if (!graph->hasNode(start) || !graph->hasNode(end)) {
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
    if (minWeightStart.size() == graph->getEdgeCount() + 1
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

    if (!graph->hasNode(start) || !graph->hasNode(end)) {
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
    if (minWeightStart.size() == graph->getEdgeCount() + 1
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
    Path* path = dijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, validNodes.empty() ? nullptr : &validNodes);

    // clear memory
    unordered_set<Node *>().swap(validNodes);

    return path;
}

Path *findMinElevationUsingEdgeBasedDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    unordered_set<Node::Edge *> validEdges = getReachableEdges(graph, start, end, maxWeightRatio);

    Path* path = edgeBasedDijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, validEdges.empty() ? nullptr : &validEdges);

    return path;
}

PathEdges* removeDuplicateDirectedEdge(PathEdges* pathEdges, bool prioritizeRestrictions){
    if (pathEdges->getPathEdges().empty()){
        return pathEdges;
    }

    // find the last index of each edge
    unordered_map<Node::Edge*, int> umap;
    int i = 0;
    vector<PathEdges::PathEdge*> &pathEdgesVector = pathEdges->getPathEdges();
    for (PathEdges::PathEdge* pathEdge: pathEdgesVector){
        umap[pathEdge->getEdge()] = i++;
    }

    vector<pair<int,int>> intervals;
    int l(0);
    for (int r=0; r<pathEdgesVector.size(); r++){
        // todo
        i = umap[pathEdgesVector[r]->getEdge()];
        if (r != i){
            if (!prioritizeRestrictions || r == 0 || i + 1 == pathEdgesVector.size()
                || !Graph::isRestricted(pathEdgesVector[r-1]->getEdge(), pathEdgesVector[i]->getEdge(), false)){
                intervals.emplace_back(l, r);
                l = r = i;
            }
        }
    }

    if (intervals.empty()){
        return pathEdges;
    }
    intervals.emplace_back(l, pathEdgesVector.size());
    PathEdges *_pathEdges = new PathEdges(pathEdges->getStart());
    for (pair<int, int>& p: intervals){
        for (int i=p.first; i<p.second; i++){
            _pathEdges->addEdge(pathEdgesVector[i]->getEdge());
        }
    }

    delete pathEdges;

    return _pathEdges;
}


PathEdges* removeDuplicateUndirectedEdge(PathEdges* pathEdges, bool prioritizeRestrictions){

    if (pathEdges->getPathEdges().empty()){
        return pathEdges;
    }

    // find the last index of each edge
    unordered_map<long, int> umap;
    int i = 0;
    vector<PathEdges::PathEdge*> &pathEdgesVector = pathEdges->getPathEdges();
    for (PathEdges::PathEdge* pathEdge: pathEdgesVector){
        umap[pathEdge->getEdge()->getId()] = i++;
    }

    vector<pair<int,int>> intervals;
    int l(0);
    bool remove = true;
    for (int r=0; r<pathEdgesVector.size();){
        i = umap[pathEdgesVector[r]->getEdge()->getId()];
        if (r != i){
            if (pathEdgesVector[i]->getEdge() == pathEdgesVector[r]->getEdge()){
                if (r == 0 || i+1 == pathEdgesVector.size()
                    || !Graph::isRestricted(pathEdgesVector[r-1]->getEdge(), pathEdgesVector[i]->getEdge(), false)){
                    intervals.emplace_back(l, r);
                    l = r = i;
                    remove = false;
                } else {
                    r++;
                }
            } else {
                if (r == 0 || i+1 == pathEdgesVector.size()
                    || !Graph::isRestricted(pathEdgesVector[r-1]->getEdge(), pathEdgesVector[i+1]->getEdge(), false)){
                    if (l < r) {
                        intervals.emplace_back(l, r);
                    }
                    l = r = i + 1;
                    remove = false;
                } else {
                    r++;
                }
            }
        } else {
            r++;
        }
    }

    if (remove){
        return pathEdges;
    }

    intervals.emplace_back(l, pathEdgesVector.size());

    PathEdges *_pathEdges = new PathEdges(pathEdges->getStart());
    for (pair<int, int>& p: intervals){

        for (int i=p.first; i<p.second; i++){
            _pathEdges->addEdge(pathEdgesVector[i]->getEdge());
        }
    }

    delete pathEdges;

    return _pathEdges;
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
    bool useVisitedEdgeId = duplicateEdge == minimizeDuplicateUndirectedEdges;
    bool useVisitedEdge = duplicateEdge == minimizeDuplicateDirectedEdges;
    bool useVisited = useVisitedEdgeId || useVisitedEdge;

    // for minimization: choose first edge randomly, then, follow the min elevation path if possible
    bool firstEdge = true;
    PathEdges* _pathEdges;
    Node::Edge *prevEdge;
    Node::Edge *curEdge;
    if (rand() & 1) {
        if (pathEdges != nullptr) {
            int index = 0;
            if (!pathEdges->empty() && useVisited) {
                int targetCutIndex = rand() % (pathEdges->size());
                for (PathEdges::PathEdge *pathEdge: pathEdges->getPathEdges()) {
                    if (index == targetCutIndex
                        || (duplicateEdge == minimizeDuplicateUndirectedEdges
                            && visitedEdgeId.find(pathEdge->getEdge()->getId()) != visitedEdgeId.end())
                        || (duplicateEdge == minimizeDuplicateDirectedEdges
                            && visitedEdge.find(pathEdge->getEdge()) != visitedEdge.end())) {
                        break;
                    }
                    index++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert(pathEdge->getEdge()->getId());
                    } else if (useVisitedEdge) {
                        visitedEdge.insert(pathEdge->getEdge());
                    }
                }
            }
            _pathEdges = pathEdges->cutBefore(index);
        } else {
            _pathEdges = new PathEdges(start);
        }

        prevEdge = _pathEdges->getLastEdge();
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge){
                    firstEdge = false;
                } else {
                    curEdge = prevEdgeMapMinElevationEnd[_pathEdges->getEnd()];
                    if (curEdge != nullptr) {
                        auto it = minLengthEnd.find(curEdge->getV());
                        if (it != minLengthEnd.end() &&
                            _pathEdges->getLength() + curEdge->getLength() + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge){
                vector<Node::Edge*> nextEdgeCandidates;
                vector<Node::Edge*> nextEdgeCandidatesVisited;
                for (pair<const long, Node::Edge *> pair1: _pathEdges->getEnd()->getEdges()) {
                    curEdge = pair1.second;

                    auto it = minLengthEnd.find(curEdge->getV());

                    if (it == minLengthEnd.end() || _pathEdges->getLength() + curEdge->getLength() + it->second > maxLength) {
                        continue;
                    }

                    if (duplicateEdge == ignoreDuplicateEdges
                        || (useVisitedEdgeId && visitedEdgeId.find(curEdge->getId()) == visitedEdgeId.end())
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())){
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited){
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }
                }

                if (_pathEdges->getEnd() == end){
                    break;
                }

                // visit unvisited node first if not ignoring duplicate edges
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
            }

            _pathEdges->addEdge(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->getId());
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }

        }
    } else {
        vector<Node::Edge *> _edges;
        PathEdges::PathEdge *prevPathEdge;
        Node *prevNode;
        int rIndex = 0;
        if (pathEdges != nullptr) {

            if (!pathEdges->empty() && useVisited) {
                int targetCutRIndex = rand() % (pathEdges->size());
                for (auto it = pathEdges->getPathEdges().rbegin(); rIndex < targetCutRIndex; it++) {
                    if (duplicateEdge == minimizeDuplicateUndirectedEdges
                        && visitedEdgeId.find((*it)->getEdge()->getId()) != visitedEdgeId.end()
                        || (duplicateEdge == minimizeDuplicateDirectedEdges
                            && visitedEdge.find((*it)->getEdge()) != visitedEdge.end())) {
                        break;
                    }
                    rIndex++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert((*it)->getEdge()->getId());
                    } else if (useVisitedEdge) {
                        visitedEdge.insert((*it)->getEdge());
                    }
                }
            }
            prevPathEdge = (rIndex == 0) ? nullptr : pathEdges->at(pathEdges->size() - rIndex);
        } else {
            prevPathEdge = nullptr;
        }

        prevEdge = prevPathEdge == nullptr ? nullptr : prevPathEdge->getEdge();
        prevNode = prevEdge == nullptr ? end : prevEdge->getU();
        int length = prevPathEdge == nullptr ? 0 : pathEdges->getLength() - prevPathEdge->getLength() + prevEdge->getLength();

        while (true) {
            bool foundEdge = false;
            if (!maximize) {
                if (firstEdge) {
                    firstEdge = false;
                } else {
                    curEdge = prevEdgeMapMinElevationStart[prevNode];
                    if (curEdge != nullptr) {
                        auto it = minLengthStart.find(curEdge->getU());
                        if (it != minLengthStart.end() && length + curEdge->getLength() + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge) {
                vector<Node::Edge *> nextEdgeCandidates;
                vector<Node::Edge *> nextEdgeCandidatesVisited;
                for (pair<const long, Node::Edge *> pair1: prevNode->getReversedEdges()) {
                    curEdge = pair1.second;

                    auto it = minLengthStart.find(curEdge->getU());
                    if (it == minLengthStart.end() || length + curEdge->getLength() + it->second > maxLength) {
                        continue;
                    }

                    if (duplicateEdge == ignoreDuplicateEdges
                        || (useVisitedEdgeId && visitedEdgeId.find(curEdge->getId()) == visitedEdgeId.end())
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())) {
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited) {
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }

                }

                if (!_edges.empty() && _edges.back()->getU() == start) {
                    break;
                }

                // visit unvisited node first if minimizeRevisit == true
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
            }

            _edges.push_back(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->getId());
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
            length += prevEdge->getLength();
            prevNode = prevEdge->getU();
        }


        _pathEdges = new PathEdges(start);
        for (auto it=_edges.rbegin(); it!=_edges.rend(); it++){
            _pathEdges->addEdge(*it);
        }
        if (pathEdges != nullptr) {
            if (rIndex > 0) {
                auto it = pathEdges->getPathEdges().rbegin() + (rIndex - 1);
                while (true) {
                    _pathEdges->addEdge((*it)->getEdge());
                    if (it-- == pathEdges->getPathEdges().rbegin()) {
                        break;
                    }
                }
            }
        }
    }

    if (duplicateEdge == minimizeDuplicateUndirectedEdges) {
        _pathEdges = removeDuplicateUndirectedEdge(_pathEdges, false);
    } else if (duplicateEdge == minimizeDuplicateDirectedEdges){
        _pathEdges = removeDuplicateDirectedEdge(_pathEdges, false);
    }


    return _pathEdges;
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
    bool useVisitedEdgeId = duplicateEdge == minimizeDuplicateUndirectedEdges;
    bool useVisitedEdge = duplicateEdge == minimizeDuplicateDirectedEdges;
    bool useVisited = useVisitedEdgeId || useVisitedEdge;

    Node::Edge *prevEdge;
    Node::Edge *curEdge;
    Node *prevNode = nullptr;
    PathEdges* _pathEdges;
    bool firstEdge = true;
    if (rand() & 1) {
        if (pathEdges != nullptr) {
            int index = 0;
            if (!pathEdges->empty() && useVisited) {
                int targetCutIndex = rand() % (pathEdges->size());
                for (PathEdges::PathEdge *pathEdge: pathEdges->getPathEdges()) {
                    if (index == targetCutIndex
                        || (duplicateEdge == minimizeDuplicateUndirectedEdges
                            && visitedEdgeId.find(pathEdge->getEdge()->getId()) != visitedEdgeId.end())
                        || (duplicateEdge == minimizeDuplicateDirectedEdges
                            && visitedEdge.find(pathEdge->getEdge()) != visitedEdge.end())) {
                        break;
                    }
                    index++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert(pathEdge->getEdge()->getId());
                    } else if (useVisitedEdge) {
                        visitedEdge.insert(pathEdge->getEdge());
                    }
                }
            }
            _pathEdges = pathEdges->cutBefore(index);
        } else {
            _pathEdges = new PathEdges(start);
        }

        prevEdge = _pathEdges->getLastEdge();
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge) {
                    firstEdge = false;
                } else {
                    curEdge = prevEdgeMapMinElevationEnd[_pathEdges->getLastEdge()];
                    if (curEdge != nullptr){
                        auto it = minLengthEnd.find(curEdge);
                        if (it != minLengthEnd.end() && _pathEdges->getLength() + it->second <= maxLength) {
                            prevEdge = curEdge;
                            foundEdge = true;
                        }
                    }
                }
            }
            if (!foundEdge) {
                vector<Node::Edge *> nextEdgeCandidates;
                vector<Node::Edge *> nextEdgeCandidatesVisited;
                for (pair<const long, Node::Edge *> pair1: _pathEdges->getEnd()->getEdges()) {
                    curEdge = pair1.second;

                    if (Graph::isRestricted(prevEdge, curEdge, false)) {
                        continue;
                    }

                    auto it = minLengthEnd.find(curEdge);

                    if (it == minLengthEnd.end() || _pathEdges->getLength() + it->second > maxLength) {
                        continue;
                    }

                    if (duplicateEdge == ignoreDuplicateEdges
                        || (useVisitedEdgeId && visitedEdgeId.find(curEdge->getId()) == visitedEdgeId.end())
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())) {
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited) {
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }
                }

                if (_pathEdges->getEnd() == end) {
                    break;
                }

                // visit unvisited node first if not ignoring duplicate edges
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
            }

            _pathEdges->addEdge(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->getId());
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
        }
    } else {
        int rIndex = 0;
        vector<Node::Edge*> _edges;
        PathEdges::PathEdge *prevPathEdge;
        if (pathEdges != nullptr) {

            if (!pathEdges->empty() && useVisited) {
                int targetCutRIndex = rand() % (pathEdges->size());
                for (auto it = pathEdges->getPathEdges().rbegin(); rIndex < targetCutRIndex; it++) {
                    if (duplicateEdge == minimizeDuplicateUndirectedEdges
                        && visitedEdgeId.find((*it)->getEdge()->getId()) != visitedEdgeId.end()
                        || (duplicateEdge == minimizeDuplicateDirectedEdges
                            && visitedEdge.find((*it)->getEdge()) != visitedEdge.end())) {
                        break;
                    }
                    rIndex++;
                    if (useVisitedEdgeId) {
                        visitedEdgeId.insert((*it)->getEdge()->getId());
                    } else if (useVisitedEdge) {
                        visitedEdge.insert((*it)->getEdge());
                    }
                }
            }
            prevPathEdge = (rIndex == 0) ?
                            nullptr : pathEdges->at(pathEdges->size() - rIndex );
        } else {
            prevPathEdge = nullptr;
        }
        prevEdge = prevPathEdge == nullptr ? nullptr : prevPathEdge->getEdge();
        int length = prevPathEdge == nullptr ? 0 : pathEdges->getLength() - prevPathEdge->getLength() + prevEdge->getLength();
        prevNode = prevPathEdge == nullptr ? end : prevEdge->getU();
        while (true) {
            bool foundEdge = false;
            if (!maximize){
                if (firstEdge){
                    firstEdge = false;
                } else {
                    curEdge = prevEdgeMapMinElevationStart[prevEdge];
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
                for (pair<const long, Node::Edge *> pair1: prevNode->getReversedEdges()) {
                    curEdge = pair1.second;

                    if (Graph::isRestricted(prevEdge, curEdge, true)) {
                        continue;
                    }

                    auto it = minLengthStart.find(curEdge);
                    if (it == minLengthStart.end() || length + it->second > maxLength) {
                        continue;
                    }

                    if (duplicateEdge == ignoreDuplicateEdges
                        || (useVisitedEdgeId && visitedEdgeId.find(curEdge->getId()) == visitedEdgeId.end())
                        || (useVisitedEdge && visitedEdge.find(curEdge) == visitedEdge.end())) {
                        nextEdgeCandidates.push_back(curEdge);
                    } else if (useVisited) {
                        nextEdgeCandidatesVisited.push_back(curEdge);
                    }
                }

                if (!_edges.empty() && _edges.back()->getU() == start) {
                    break;
                }

                // visit unvisited node first if minimizeRevisit == true
                if (!nextEdgeCandidates.empty()) {
                    prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
                } else {
                    prevEdge = nextEdgeCandidatesVisited[rand() % nextEdgeCandidatesVisited.size()];
                }
            }

            _edges.push_back(prevEdge);
            if (useVisitedEdgeId) {
                visitedEdgeId.insert(prevEdge->getId());
            } else if (useVisitedEdge){
                visitedEdge.insert(prevEdge);
            }
            length += prevEdge->getLength();
            prevNode = prevEdge->getU();
        }


        _pathEdges = new PathEdges(start);
        for (auto it=_edges.rbegin(); it!=_edges.rend(); it++){
            _pathEdges->addEdge(*it);
        }
        if (pathEdges != nullptr) {
            if (rIndex > 0) {
                auto it = pathEdges->getPathEdges().rbegin() + (rIndex - 1);
                while (true) {
                    _pathEdges->addEdge((*it)->getEdge());
                    if (it-- == pathEdges->getPathEdges().rbegin()) {
                        break;
                    }
                }
            }
        }
    }

    if (duplicateEdge == minimizeDuplicateUndirectedEdges) {
        _pathEdges = removeDuplicateUndirectedEdge(_pathEdges, true);
    } else if (duplicateEdge == minimizeDuplicateDirectedEdges){
        _pathEdges = removeDuplicateDirectedEdge(_pathEdges, true);
    }


    return _pathEdges;
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
        for (auto it=pathEdgesSet.begin(); it!=pathEdgesSet.end();) {
            if (pq.size() >= numMaxSelect) {
                if ((*it)->getElevation() > pq.top()->getElevation()) {
                    PathEdges* _pathEdges = pq.top();
                    pq.pop();
                    delete _pathEdges;

                    pq.push(*it);
                    it++;
                } else {
                    PathEdges* _pathEdges = *it;
                    it = pathEdgesSet.erase(it);
                    delete _pathEdges;
                }
            } else {
                pq.push(*it);
                it++;
            }
        }
        while (!pq.empty()){
            res.insert(pq.top());
            pq.pop();
        }

    } else {
        priority_queue<PathEdges *, vector<PathEdges *>, PathEdges::smallerElevation> pq;
        for (auto it=pathEdgesSet.begin(); it!=pathEdgesSet.end();) {
            if (pq.size() >= numMaxSelect) {
                if ((*it)->getElevation() < pq.top()->getElevation()) {
                    PathEdges* _pathEdges = pq.top();
                    pq.pop();
                    delete _pathEdges;

                    pq.push(*it);
                    it++;
                } else {
                    PathEdges* _pathEdges = *it;
                    it = pathEdgesSet.erase(it);
                    delete _pathEdges;
                }
            } else {
                pq.push(*it);
                it++;
            }
        }
        while (!pq.empty()){
            res.insert(pq.top());
            pq.pop();
        }
    }

    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual>().swap(pathEdgesSet);

    return res;
}


Path* geneticAlgorithm(
        Graph *graph, Node *start, Node *end, int maxLength,
        unordered_map<Node*, int>& minLengthStart, unordered_map<Node*, int>& minLengthEnd,
        Path* shortestPath,
        size_t numProduce,
        size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch, bool maximize) {

    if (start == end){
        return new Path();
    }

    auto t0 = chrono::high_resolution_clock::now();

    PathEdges *pathEdges;
    PathEdges *pathEdges1;
    PathEdges *pathEdges2;

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
        if (dijkstraPath->getLength() <= maxLength){
            PathEdges* dijkstraPathEdges = new PathEdges(dijkstraPath);
            delete dijkstraPath;
            if (duplicateEdge == minimizeDuplicateUndirectedEdges){
                dijkstraPathEdges = removeDuplicateUndirectedEdge(dijkstraPathEdges, false);
            } else if (duplicateEdge == minimizeDuplicateDirectedEdges){
                dijkstraPathEdges = removeDuplicateDirectedEdge(dijkstraPathEdges, false);
            }
            dijkstraPath = dijkstraPathEdges->toPath();
            delete dijkstraPathEdges;
            dijkstraPath->confirmOptimal();
            return dijkstraPath;
        }
        delete dijkstraPath;
        dijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0, Node::Edge::getElevation, false, &prevEdgeMapMinElevationStart, &validNodes);
    }

    // first generation: optional shortestPathEdges and random paths
    unordered_set<PathEdges *, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    if (shortestPath != nullptr){
        pathEdgesSet.insert(new PathEdges(shortestPath));
    }
    for (size_t i = 0; i < numProduce; i++) {
        PathEdges* _pathEdges = mutatePathEdges(
                nullptr, start, end, minLengthStart, minLengthEnd,
                prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd, maxLength, duplicateEdge, maximize);

        bool inserted = pathEdgesSet.insert(_pathEdges).second;
        if (!inserted){
            delete _pathEdges;
        }
    }

    for (int epoch = 0; epoch < numEpoch; epoch++) {

        // natural selection
        pathEdgesSet = naturalSelection(pathEdgesSet, numMaxSelect, maximize);

        // prepare for crossovers: search for intersections
        unordered_map<Node::Edge *, unordered_map<PathEdges *, vector<size_t>>> intersections;
        for (auto it: pathEdgesSet) {
            size_t index = 0;
            for (PathEdges::PathEdge *pathEdge: it->getPathEdges()) {
                intersections[pathEdge->getEdge()][it].push_back(index++);
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
                PathEdges* _pathEdges = mutatePathEdges(*it, start, end, minLengthStart, minLengthEnd,
                                                          prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                          maxLength, duplicateEdge, maximize);
                bool inserted = pathEdgesSet.insert(_pathEdges).second;
                if (!inserted){
                    delete _pathEdges;
                }
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
            pathEdges1 = choicesIt1->first;
            pathEdges2 = choicesIt2->first;

            // use the second choice as the prefix if the first choice does not give any path within the length limit
            if (pathEdges1->at(choicesIt1->second.front())->getLength()
                + pathEdges2->getLength() - pathEdges2->at(choicesIt2->second.back())->getLength() > maxLength) {
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
                Node::Edge *edge = pathEdges1->at(index1)->getEdge();
                if (edge->getV()->getRestrictions() != nullptr) {
                    auto restrictionsIt = edge->getV()->getRestrictions()->find(edge);
                    if (restrictionsIt != edge->getV()->getRestrictions()->end()) {
                        if (restrictionsIt->second.find(pathEdges2->at(index2 + 1)->getEdge())
                            != restrictionsIt->second.end()) {
                            continue;
                        }
                    }
                }

                auto pathEdgesIt = pathEdges2->getPathEdges().begin();

                advance(pathEdgesIt, index2 + 1);
                while (pathEdgesIt != pathEdges2->getPathEdges().end()) {
                    pathEdges->addEdge((*pathEdgesIt)->getEdge());
                    pathEdgesIt++;
                }
            }

            if (maximize || rand() & 1) {
                PathEdges* _pathEdges = mutatePathEdges(pathEdges, start, end, minLengthStart, minLengthEnd,
                                                          prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                          maxLength, duplicateEdge, maximize);
                bool inserted = pathEdgesSet.insert(_pathEdges).second;
                delete pathEdges;
                if (!inserted){
                    delete _pathEdges;
                }
            } else {
                bool inserted = pathEdgesSet.insert(pathEdges).second;
                if (!inserted){
                    delete pathEdges;
                }
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

    Path* path = pathEdges->toPath();

    for (auto it=pathEdgesSet.begin();it!=pathEdgesSet.end();){
        PathEdges* _pathEdgesToDelete = *it;
        it = pathEdgesSet.erase(it);
        delete _pathEdgesToDelete;
    }

    return path;

}

Path* edgeBasedGeneticAlgorithm(
        Graph *graph, Node *start, Node *end, int maxLength,
        unordered_map<Node::Edge*, int>& minLengthStart, unordered_map<Node::Edge*, int>& minLengthEnd,
        Path* shortestPath,
        size_t numProduce, size_t numMaxSelect, int numEpoch, DuplicateEdge duplicateEdge,
        int maxMilliseconds, int minEpoch, bool maximize
){

    if (start == end){
        return new Path();
    }

    auto t0 = chrono::high_resolution_clock::now();
    PathEdges* pathEdges;
    PathEdges *pathEdges1;
    PathEdges *pathEdges2;
    Node::Edge* edge;
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMapMinElevationStart;
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMapMinElevationEnd;

    if (!maximize){
        unordered_set<Node::Edge*> validEdges;
        for (pair<Node::Edge *const, int> &pair1: minLengthEnd) {
            auto it = minLengthStart.find(pair1.first);
            if (it == minLengthStart.end())
                continue;
            if (pair1.second + it->second - pair1.first->getLength() > maxLength)
                continue;
            validEdges.insert(pair1.first);
        }

        Node::Edge *lastEdge = nullptr;
        int curMinWeight = INT_MAX;

        // find the path reversed so that the path can be reconstructed non-reversed
        edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0,
                Node::Edge::getElevation, true,
                &prevEdgeMapMinElevationEnd, lastEdge, &validEdges);
        Path* dijkstraPath = edgeBasedPathFromPrevEdgeMap(lastEdge, prevEdgeMapMinElevationEnd);
        if (dijkstraPath->getLength() < maxLength){
            PathEdges* dijkstraPathEdges = new PathEdges(dijkstraPath);
            delete dijkstraPath;
            if (duplicateEdge == minimizeDuplicateUndirectedEdges){
                dijkstraPathEdges = removeDuplicateUndirectedEdge(dijkstraPathEdges, true);
            } else if (duplicateEdge == minimizeDuplicateDirectedEdges){
                dijkstraPathEdges = removeDuplicateDirectedEdge(dijkstraPathEdges, true);
            }

            // dijkstraPathEdges may be empty if there is a duplicate edge that cannot be removed due to road restrictions
            if (!dijkstraPathEdges->empty()) {
                dijkstraPath = dijkstraPathEdges->toPath();
                delete dijkstraPathEdges;
                dijkstraPath->confirmOptimal();
                return dijkstraPath;
            }
        }
        delete dijkstraPath;
        edgeBasedDijkstraAlgorithm(
                graph, start, end, curMinWeight, -1.0, Node::Edge::getElevation, false, &prevEdgeMapMinElevationStart, lastEdge, &validEdges);
    }

    // First generation
    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    if (shortestPath != nullptr) {
        pathEdgesSet.insert(new PathEdges(shortestPath));
    }
    for (size_t i=0; i<numProduce; i++){
        PathEdges* _pathEdges = edgeBasedMutatePathEdges(nullptr, start, end, minLengthStart, minLengthEnd,
                                                           prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                           maxLength, duplicateEdge, maximize);

        bool inserted = pathEdgesSet.insert(_pathEdges).second;
        if (!inserted){
            delete _pathEdges;
        }
    }

    for (int epoch=0; epoch < numEpoch; epoch++){

        // natural selection
        pathEdgesSet = naturalSelection(pathEdgesSet, numMaxSelect, maximize);

        // prepare for crossovers: search for intersections
        unordered_map<Node::Edge*, unordered_map<PathEdges*, vector<size_t>>> intersections;
        for (auto it: pathEdgesSet){
            size_t index = 0;
            for (PathEdges::PathEdge* pathEdge: it->getPathEdges()){
                intersections[pathEdge->getEdge()][it].push_back(index++);
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
                PathEdges* _pathEdges = edgeBasedMutatePathEdges(*it, start, end, minLengthStart, minLengthEnd,
                                                                   prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                                   maxLength, duplicateEdge, maximize);

                bool inserted = pathEdgesSet.insert(_pathEdges).second;
                if (!inserted){
                    delete _pathEdges;
                }
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
            pathEdges1 = choicesIt1->first;
            pathEdges2 = choicesIt2->first;

            // use the second choice as the prefix if the first choice does not give any path within the length limit
            if (pathEdges1->at(choicesIt1->second.front())->getLength()
                + pathEdges2->getLength() - pathEdges2->at(choicesIt2->second.back())->getLength() > maxLength) {
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
                edge = pathEdges1->at(index1)->getEdge();
                if (edge->getV()->getRestrictions() != nullptr) {
                    auto restrictionsIt = edge->getV()->getRestrictions()->find(edge);
                    if (restrictionsIt != edge->getV()->getRestrictions()->end()) {
                        if (restrictionsIt->second.find(pathEdges2->at(index2 + 1)->getEdge())
                            != restrictionsIt->second.end()) {
                            continue;
                        }
                    }
                }

                auto pathEdgesIt = pathEdges2->getPathEdges().begin();

                advance(pathEdgesIt, index2 + 1);
                while (pathEdgesIt != pathEdges2->getPathEdges().end()) {
                    pathEdges->addEdge((*pathEdgesIt)->getEdge());
                    pathEdgesIt++;
                }
            }

            if (maximize || rand() & 1) {
                PathEdges* _pathEdges = edgeBasedMutatePathEdges(pathEdges, start, end, minLengthStart, minLengthEnd,
                                                                prevEdgeMapMinElevationStart, prevEdgeMapMinElevationEnd,
                                                                maxLength, duplicateEdge, maximize);
                delete pathEdges;
                bool inserted = pathEdgesSet.insert(_pathEdges).second;
                if (!inserted){
                    delete _pathEdges;
                }

           } else {
                bool inserted = pathEdgesSet.insert(pathEdges).second;
                if (!inserted){
                    delete pathEdges;
                }
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

    Path* path = pathEdges->toPath();

    for (auto it=pathEdgesSet.begin();it!=pathEdgesSet.end();){
        PathEdges* _pathEdgesToDelete = *it;
        it = pathEdgesSet.erase(it);
        delete _pathEdgesToDelete;
    }

    return path;

}