#include <float.h>
#include <unordered_map>
#include <random>

#include "pathfinding.h"



using namespace std;



unordered_map<Node::Edge *, double> edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        double &curMinWeight,
        double maxWeightRatio,
        double getWeight(Node::Edge *edge),
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
        double weight;
        Node::Edge *prevEdge;

        Candidate(double weight, Node::Edge *prevEdge) {
            this->weight = weight;
            this->prevEdge = prevEdge;
        }
    };

    struct comparator {
        bool operator()(const Candidate a, const Candidate b) {
            return a.weight > b.weight;
        }
    };

    unordered_map<Node::Edge *, double> minWeightEdge;
    unordered_set<Node::Edge*> visited;

    priority_queue<Candidate, vector<Candidate>, comparator> pq;

    Node::Edge *curEdge = new Node::Edge(reversed ? end : start, reversed);
    minWeightEdge[curEdge] = 0.0;
    pq.emplace(0.0, curEdge);

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

            double curWeight = candidate.weight + getWeight(curEdge);
            if (maxWeightRatio > 0.0 && curWeight > curMinWeight * maxWeightRatio) {
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
                        lastEdge = curEdge;
                    }
                }
            }
        }
    }
    return minWeightEdge;
}

Path *edgeBasedDijkstraAlgorithm(
        Graph *graph,
        Node *start,
        Node *end,
        double getWeight(Node::Edge *edge),
        bool reversed,
        unordered_set<Node::Edge *> *validEdges
) {
    Path* path = new Path();
    if (start == end) {
        return path;
    }
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMap;
    Node::Edge *edge = nullptr;
    double curMinWeight = DBL_MAX;
    // find the path reversed so that the path can be reconstructed non-reversed
    unordered_map<Node::Edge *, double> minWeightEdge = edgeBasedDijkstraAlgorithm(graph, start, end, curMinWeight,
                                                                                   -1.0, getWeight, !reversed,
                                                                                   &prevEdgeMap, edge, validEdges);


    if (edge == nullptr){
        throw runtime_error("path not found.");
    }

    while (edge->id != -1) {
        path->addEdge(edge);
        edge = prevEdgeMap[edge];
    }

    return path;
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
        throw invalid_argument("nodes do not belong to the same graph");
    }

    double curMinWeight = DBL_MAX;
    Node::Edge *lastEdge = nullptr;
    double minWeight1 = DBL_MAX;
    unordered_map<Node::Edge *, double> minWeightStart = edgeBasedDijkstraAlgorithm(
            graph, start, end, minWeight1, maxWeightRatio, Node::Edge::getLength, false, nullptr, lastEdge, nullptr
    );
    if (lastEdge == nullptr) {
        throw runtime_error("runtime error");
    }

    lastEdge = nullptr;
    double minWeight2 = DBL_MAX;
    unordered_map<Node::Edge *, double> minWeightEnd = edgeBasedDijkstraAlgorithm(
            graph, start, end, minWeight2, maxWeightRatio, Node::Edge::getLength, true, nullptr, lastEdge, nullptr
    );

    if (lastEdge == nullptr || abs(minWeight1 - minWeight2) > 0.001) {
        throw runtime_error("runtime error");
    }


    unordered_set<Node::Edge *> validEdges;
    double maxWeight = curMinWeight * maxWeightRatio;
    for (pair<Node::Edge *const, double> &pair1: minWeightStart) {
        auto it = minWeightEnd.find(pair1.first);
        if (it == minWeightEnd.end())
            continue;

        if (pair1.second + it->second - Node::Edge::getLength(pair1.first) > maxWeight)
            continue;

        validEdges.insert(pair1.first);
    }
    cout << "minWeightStart: " << minWeightStart.size() << ", minWeightEnd: " << minWeightEnd.size() << endl;
    cout << "validEdges: " << validEdges.size() << endl;

    return validEdges;
}



Path *elenaPathFindMinUsingDijkstra(
        Graph *graph,
        Node *start,
        Node *end,
        double maxWeightRatio
) {

    unordered_set<Node::Edge *> validEdges = getReachableEdges(graph, start, end, maxWeightRatio);

    return edgeBasedDijkstraAlgorithm(graph, start, end, Node::Edge::getElevation, false, &validEdges);
}

PathEdges* mutatePathEdges(
        Graph *graph, PathEdges* pathEdges, Node* end, unordered_map<Node::Edge *, double>& minLengthEnd, double maxLength){
    pathEdges = pathEdges->randomCut();

    while (true) {
        vector<Node::Edge *> nextEdgeCandidates;
        Node::Edge* prevEdge = pathEdges->lastEdge();
        for (pair<const long, Node::Edge *> pair1: pathEdges->getEnd()->edges) {
            Node::Edge* curEdge = pair1.second;

            if (Graph::isRestricted(prevEdge, curEdge, false)) {
                continue;
            }

            auto it = minLengthEnd.find(curEdge);
            if (it == minLengthEnd.end() || pathEdges->length + it->second > maxLength) {
                continue;
            }
            nextEdgeCandidates.push_back(curEdge);
        }

        // even if current location is the target location, make moves if possible to maximize elevation gain
        if (nextEdgeCandidates.empty())
            break;

        prevEdge = nextEdgeCandidates[rand() % nextEdgeCandidates.size()];
        pathEdges->addEdge(prevEdge);
    }
    return pathEdges;
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

PathEdges* elenaPathSearchMaxUsingGeneticAlgorithm(Graph *graph, Node *start, Node *end, double maxLengthRatio, size_t numProduce, size_t numMaxSelect, int numEpoch){

    // 1) minWeightEnd and the shortest path

    double curMinWeight = DBL_MAX;
    Node::Edge *lastEdge = nullptr;
    unordered_map<Node::Edge *, Node::Edge *> prevEdgeMap;
    unordered_map<Node::Edge *, double> minWeightEnd = edgeBasedDijkstraAlgorithm(
            graph, start, end, curMinWeight, maxLengthRatio, Node::Edge::getLength, true, &prevEdgeMap, lastEdge, nullptr
    );
    if (lastEdge == nullptr){
        throw runtime_error("path not found.");
    }
    double maxWeight = curMinWeight * maxLengthRatio;

    // Find the shortest path
    PathEdges* pathEdges = new PathEdges(start);
    while (lastEdge->id != -1) {
        pathEdges->addEdge(lastEdge);
        lastEdge = prevEdgeMap[lastEdge];
    }

    // First generation: mutations from the shortest path
    unordered_set<PathEdges*, pathEdgesHash, pathEdgesEqual> pathEdgesSet;
    for (size_t i=0; i<numProduce; i++){
        pathEdgesSet.insert(mutatePathEdges(graph, pathEdges, end, minWeightEnd, maxWeight));
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

        // 4) crossovers: based on intersections find valid paths, generate mutants of each new crossovers with m% probability
        // Choose an intersected edge
        for (int produce=0; produce<numProduce; produce++) {
            auto intersectionsIt = intersections.begin();
            if (intersections.size() == 0){
                auto it = pathEdgesSet.begin();
                advance(it, rand() % pathEdgesSet.size());
                pathEdgesSet.insert(mutatePathEdges(graph, *it, end, minWeightEnd, maxWeight));
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
            bool flag = false;
            if (pathEdges1->at(choicesIt1->second.front())->length
                + pathEdges2->length - pathEdges2->at(choicesIt2->second.front())->length > maxWeight) {
                swap(pathEdges1, pathEdges2);
                flag = true;
            }
            vector<size_t>& indices1 = flag ? choicesIt2->second : choicesIt1->second;
            vector<size_t>& indices2 = flag ? choicesIt1->second : choicesIt2->second;

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

            index1 = 0;
            index2 = indices2.size()-1;

            // concatenate two paths
            pathEdges = pathEdges1->cutAfter(indices1[index1]);
            if (indices2[index2] + 1 < pathEdges2->size()) {

                // check if restricted
                Node::Edge* edge = pathEdges1->at(indices1[index1])->edge;
                if (edge->v->restrictions != nullptr) {
                    auto restrictionsIt = edge->v->restrictions->find(edge);
                    if (restrictionsIt != edge->v->restrictions->end()) {
                        if (restrictionsIt->second.find(pathEdges2->at(indices2[index2] + 1)->edge)
                            != restrictionsIt->second.end()) {
                            continue;
                        }
                    }
                }


                auto pathEdgesIt = pathEdges2->pathEdges.begin();
                advance(pathEdgesIt, indices2[index2] + 1);
                while (pathEdgesIt != pathEdges2->pathEdges.end()) {
                    pathEdges->addEdge((*pathEdgesIt)->edge);
                    pathEdgesIt++;
                }
            } else if (pathEdges->getEnd() != end)
                throw runtime_error("This is not supposed to happen");

            if (pathEdges->getLength() > maxWeight + 1e-10)
                throw runtime_error("This is not supposed to happen");


            pathEdgesSet.insert(mutatePathEdges(graph, pathEdges, end, minWeightEnd, maxWeight));
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
